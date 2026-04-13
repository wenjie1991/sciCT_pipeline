#include <algorithm>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int GZ_BUFFER_SIZE = 1 << 20;

struct Args {
    std::string input;
    std::string output;
    std::string matrix;
    std::string suffix;
};

struct FastqRecord {
    std::string header;
    std::string seq;
    std::string plus;
    std::string qual;
};

std::string trim_newline(std::string line) {
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
        line.pop_back();
    }
    return line;
}

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string current;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        char ch = line[i];
        if (ch == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                current.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (ch == ',' && !in_quotes) {
            fields.push_back(current);
            current.clear();
        } else {
            current.push_back(ch);
        }
    }

    fields.push_back(current);
    if (!fields.empty() && fields[0].size() >= 3 &&
        static_cast<unsigned char>(fields[0][0]) == 0xEF &&
        static_cast<unsigned char>(fields[0][1]) == 0xBB &&
        static_cast<unsigned char>(fields[0][2]) == 0xBF) {
        fields[0].erase(0, 3);
    }

    for (auto& field : fields) {
        field = trim_newline(field);
    }

    return fields;
}

int index_of(const std::vector<std::string>& header, std::string_view key) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == key) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

using BarcodeMap = std::unordered_map<std::string, std::string>;

BarcodeMap load_well_id_matrix(const std::string& matrix_path) {
    std::ifstream input(matrix_path);
    if (!input) {
        throw std::runtime_error("Failed to open matrix CSV: " + matrix_path);
    }

    std::string header_line;
    if (!std::getline(input, header_line)) {
        throw std::runtime_error("Matrix CSV is empty: " + matrix_path);
    }

    const auto header = split_csv_line(header_line);
    const int p1s7 = index_of(header, "PAGE-1-s7");
    const int p1s5 = index_of(header, "PAGE-1-s5");
    const int p2s7 = index_of(header, "PAGE-2-s7");
    const int p2s5 = index_of(header, "PAGE-2-s5");
    const int well = index_of(header, "Well-ID");

    if (p1s7 < 0 || p1s5 < 0 || p2s7 < 0 || p2s5 < 0 || well < 0) {
        throw std::runtime_error("Matrix CSV is missing required columns");
    }

    BarcodeMap barcode_to_well;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        const auto row = split_csv_line(line);
        const size_t required = static_cast<size_t>(std::max({p1s7, p1s5, p2s7, p2s5, well}) + 1);
        if (row.size() < required) {
            continue;
        }

        const std::string& page1s7 = row[p1s7];
        const std::string& page1s5 = row[p1s5];
        const std::string& page2s7 = row[p2s7];
        const std::string& page2s5 = row[p2s5];
        const std::string& well_id = row[well];

        barcode_to_well[page1s7 + '\t' + page1s5] = well_id;
        barcode_to_well[page1s7 + '\t' + page2s5] = well_id;
        barcode_to_well[page2s7 + '\t' + page2s5] = well_id;
        barcode_to_well[page2s7 + '\t' + page1s5] = well_id;
    }

    return barcode_to_well;
}

bool gz_readline(gzFile file, std::string& out) {
    out.clear();
    char buffer[GZ_BUFFER_SIZE];

    while (true) {
        char* result = gzgets(file, buffer, GZ_BUFFER_SIZE);
        if (result == nullptr) {
            return !out.empty();
        }

        out.append(buffer);
        const size_t len = out.size();
        if (len > 0 && out[len - 1] == '\n') {
            return true;
        }

        if (gzeof(file)) {
            return true;
        }
    }
}

bool read_record(gzFile input, FastqRecord& record) {
    if (!gz_readline(input, record.header)) {
        return false;
    }

    if (!gz_readline(input, record.seq) ||
        !gz_readline(input, record.plus) ||
        !gz_readline(input, record.qual)) {
        throw std::runtime_error("Malformed FASTQ: incomplete record");
    }

    return true;
}

void gz_write_all(gzFile output, const std::string& text) {
    const int written = gzwrite(output, text.data(), static_cast<unsigned int>(text.size()));
    if (written == 0) {
        int err_no = Z_OK;
        const char* err_msg = gzerror(output, &err_no);
        throw std::runtime_error("Failed to write gzip output: " + std::string(err_msg ? err_msg : "unknown error"));
    }
}

std::string rewrite_header(const std::string& original, const BarcodeMap& barcode_map, const std::string& suffix) {
    const std::string header = trim_newline(original);
    const size_t first_space = header.find(' ');
    const std::string first_part = header.substr(0, first_space);
    const std::string illumina_part = first_space == std::string::npos ? "" : header.substr(first_space + 1);

    std::vector<size_t> underscores;
    underscores.reserve(8);
    for (size_t i = 0; i < first_part.size(); ++i) {
        if (first_part[i] == '_') {
            underscores.push_back(i);
        }
    }

    if (underscores.size() < 4) {
        return original;
    }

    const size_t u0 = underscores[underscores.size() - 4];
    const size_t u1 = underscores[underscores.size() - 3];
    const size_t u2 = underscores[underscores.size() - 2];
    const size_t u3 = underscores[underscores.size() - 1];

    const std::string i7 = first_part.substr(u0 + 1, u1 - u0 - 1);
    const std::string i5 = first_part.substr(u1 + 1, u2 - u1 - 1);
    const std::string s7 = first_part.substr(u2 + 1, u3 - u2 - 1);
    const std::string s5_raw = first_part.substr(u3 + 1);
    const size_t read_suffix_pos = s5_raw.find('/');
    const std::string s5 = read_suffix_pos == std::string::npos ? s5_raw : s5_raw.substr(0, read_suffix_pos);
    const std::string read_suffix = read_suffix_pos == std::string::npos ? "" : s5_raw.substr(read_suffix_pos);

    const auto it = barcode_map.find(s7 + '\t' + s5);
    if (it == barcode_map.end()) {
        return original;
    }

    std::string rewritten = first_part.substr(0, u0);
    rewritten.push_back(':');
    rewritten += i7;
    rewritten.push_back('_');
    rewritten += i5;
    rewritten.push_back('_');
    rewritten += it->second;
    rewritten += suffix;
    rewritten += read_suffix;

    if (!illumina_part.empty()) {
        rewritten.push_back(' ');
        rewritten += illumina_part;
    }

    rewritten.push_back('\n');
    return rewritten;
}

Args parse_args(int argc, char* argv[]) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string current = argv[i];
        if (current.rfind("--input=", 0) == 0) {
            args.input = current.substr(8);
        } else if (current == "--input" && i + 1 < argc) {
            args.input = argv[++i];
        } else if (current.rfind("--output=", 0) == 0) {
            args.output = current.substr(9);
        } else if (current == "--output" && i + 1 < argc) {
            args.output = argv[++i];
        } else if (current.rfind("--matrix=", 0) == 0) {
            args.matrix = current.substr(9);
        } else if (current == "--matrix" && i + 1 < argc) {
            args.matrix = argv[++i];
        } else if (current.rfind("--suffix=", 0) == 0) {
            args.suffix = current.substr(9);
        } else if (current == "--suffix" && i + 1 < argc) {
            args.suffix = argv[++i];
        } else if (current == "--help" || current == "-h") {
            std::cout << "Usage: rewrite_fastq_barcodes_cpp --input IN.fq.gz --output OUT.fq.gz --matrix matrix.csv [--suffix SUFFIX]\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown or incomplete argument: " + current);
        }
    }

    if (args.input.empty() || args.output.empty() || args.matrix.empty()) {
        throw std::runtime_error("Missing required arguments: --input, --output, --matrix");
    }

    return args;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Args args = parse_args(argc, argv);
        const BarcodeMap barcode_map = load_well_id_matrix(args.matrix);

        gzFile input = gzopen(args.input.c_str(), "rb");
        if (input == nullptr) {
            throw std::runtime_error("Failed to open input FASTQ: " + args.input);
        }

        gzFile output = gzopen(args.output.c_str(), "wb");
        if (output == nullptr) {
            gzclose(input);
            throw std::runtime_error("Failed to open output FASTQ: " + args.output);
        }

        gzbuffer(input, GZ_BUFFER_SIZE);
        gzbuffer(output, GZ_BUFFER_SIZE);

        FastqRecord record;
        while (read_record(input, record)) {
            const std::string new_header = rewrite_header(record.header, barcode_map, args.suffix);
            gz_write_all(output, new_header);
            gz_write_all(output, record.seq);
            gz_write_all(output, record.plus);
            gz_write_all(output, record.qual);
        }

        gzclose(output);
        gzclose(input);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
