#!/usr/bin/env bash

set -euo pipefail

input_fastq="$1"
output_fastq="$2"

base_name="$(basename "$input_fastq")"
if [[ "$base_name" =~ (^|_)([RI][12])([._]|$) ]]; then
    read_pair="${BASH_REMATCH[2]}"
else
    echo "Could not infer read/index type from filename: $base_name" >&2
    exit 1
fi

pigz -dc "$input_fastq" | \
awk -v readn="$read_pair" '
BEGIN { OFS="" }
NR % 4 == 1 {
    header = $0
    sub(/\r$/, "", header)

    split(header, parts, " ")
    first = parts[1]

    if (length(parts) > 1) {
        rest = substr(header, length(first) + 2)
        $0 = first " " readn ":" rest
    } else {
        $0 = first " " readn ":N:0:1"
    }
}
{ print }
' | pigz > "$output_fastq"
