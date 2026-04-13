#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
src="${project_dir}/src/rewrite_fastq_barcodes.cpp"
out="${project_dir}/bin/rewrite_fastq_barcodes_cpp"

mkdir -p "${project_dir}/bin"

compiler="${CXX:-}"
if [[ -z "${compiler}" ]]; then
    if command -v g++ >/dev/null 2>&1; then
        compiler="g++"
    elif command -v x86_64-conda-linux-gnu-c++ >/dev/null 2>&1; then
        compiler="x86_64-conda-linux-gnu-c++"
    else
        echo "No C++ compiler found. Install g++ or use the conda profile." >&2
        exit 1
    fi
fi

"${compiler}" -O3 -std=c++17 -o "${out}" "${src}" -lz

chmod +x "${out}"
