#!/usr/bin/env bash

set -euo pipefail

input_fastq="$1"
output_fastq="$2"
swap_index_ends="${3:-false}"
is_index_read="false"

if command -v pigz >/dev/null 2>&1; then
    DECOMPRESS_CMD=(pigz -dc)
    COMPRESS_CMD=(pigz)
else
    DECOMPRESS_CMD=(gzip -dc)
    COMPRESS_CMD=(gzip -c)
fi

base_name="$(basename "$input_fastq")"
if [[ "$base_name" =~ (^|_)([RI])([12])([._]|$) ]]; then
    read_group="${BASH_REMATCH[2]}"
    read_pair="${BASH_REMATCH[3]}"
    if [[ "$read_group" == "I" ]]; then
        is_index_read="true"
    fi
else
    echo "Could not infer read/index type from filename: $base_name" >&2
    exit 1
fi

"${DECOMPRESS_CMD[@]}" "$input_fastq" | \
awk -v readn="$read_pair" -v swap_index_ends="$swap_index_ends" -v is_index_read="$is_index_read" '
BEGIN { FS = OFS = ":" }
NR % 4 == 1 {
    sub(/\r$/, "", $0)

    # Match the original sciCT behavior: drop the 8th colon-delimited field
    # and replace the resulting double colon with the normalized read label.
    if (NF >= 8) {
        $8 = ""
        sub(/::/, " " readn ":", $0)
    } else {
        sub(/[[:space:]].*$/, "", $0)
        $0 = $0 " " readn ":1:N:0:1"
    }
}
NR % 4 == 2 {
    if (swap_index_ends == "true" && is_index_read == "true" && length($0) >= 16) {
        prefix = substr($0, 1, 8)
        middle = substr($0, 9, length($0) - 16)
        suffix = substr($0, length($0) - 7, 8)
        $0 = suffix middle prefix
    }
}
NR % 4 == 0 {
    if (swap_index_ends == "true" && is_index_read == "true" && length($0) >= 16) {
        prefix = substr($0, 1, 8)
        middle = substr($0, 9, length($0) - 16)
        suffix = substr($0, length($0) - 7, 8)
        $0 = suffix middle prefix
    }
}
{ print }
' | "${COMPRESS_CMD[@]}" > "$output_fastq"
