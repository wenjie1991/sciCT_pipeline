#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path


def load_well_id_matrix(matrix_csv):
    barcode_to_well = {}
    with open(matrix_csv, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            combos = [
                (row["PAGE-1-s7"], row["PAGE-1-s5"]),
                (row["PAGE-1-s7"], row["PAGE-2-s5"]),
                (row["PAGE-2-s7"], row["PAGE-2-s5"]),
                (row["PAGE-2-s7"], row["PAGE-1-s5"]),
            ]
            for combo in combos:
                barcode_to_well[combo] = row["Well-ID"]
    return barcode_to_well


def parse_barcodes_from_header(header):
    first_part = header.strip().split(" ")[0]
    parts = first_part.split("_")
    if len(parts) >= 4:
        return parts[-4], parts[-3], parts[-2], parts[-1]
    return None, None, None, None


def rewrite_fastq(input_fastq_gz, output_fastq_gz, barcode_dict, suffix=""):
    with gzip.open(input_fastq_gz, "rt") as infile, gzip.open(output_fastq_gz, "wt") as outfile:
        while True:
            header = infile.readline()
            if not header:
                break

            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()

            i7, i5, s7, s5 = parse_barcodes_from_header(header)
            new_header = header

            if s7 and s5:
                read_suffix = ""
                if "/" in s5:
                    s5, read_suffix = s5.split("/", 1)
                    read_suffix = f"/{read_suffix}"
                well_id = barcode_dict.get((s7, s5))
                if well_id:
                    prefix = header.strip().split(" ")[0].rsplit("_", 4)[0]
                    illumina_part = header.strip().split(" ", 1)[-1] if " " in header else ""
                    if illumina_part:
                        new_header = f"{prefix}:{i7}_{i5}_{well_id}{suffix}{read_suffix} {illumina_part}\n"
                    else:
                        new_header = f"{prefix}:{i7}_{i5}_{well_id}{suffix}{read_suffix}\n"

            outfile.write(new_header)
            outfile.write(seq)
            outfile.write(plus)
            outfile.write(qual)


def main():
    parser = argparse.ArgumentParser(
        description="Rewrite FASTQ headers by replacing s7/s5 barcodes with Well-ID."
    )
    parser.add_argument("--input", required=True, help="Input FASTQ(.gz) file")
    parser.add_argument("--output", required=True, help="Output FASTQ(.gz) file")
    parser.add_argument("--matrix", required=True, help="Barcode matrix CSV")
    parser.add_argument("--suffix", default="", help="Optional suffix appended to Well-ID")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    barcode_dict = load_well_id_matrix(args.matrix)
    rewrite_fastq(args.input, args.output, barcode_dict, suffix=args.suffix)


if __name__ == "__main__":
    main()
