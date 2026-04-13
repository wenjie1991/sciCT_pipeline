#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
image="${script_dir}/cuttag-preprocess.sif"
definition="${script_dir}/cuttag-preprocess.def"
mode="${1:-sudo}"

if command -v apptainer >/dev/null 2>&1; then
    runtime="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    runtime="singularity"
else
    echo "Neither apptainer nor singularity is available in PATH." >&2
    exit 1
fi

case "${mode}" in
    sudo)
        exec sudo "${runtime}" build "${image}" "${definition}"
        ;;
    remote)
        exec "${runtime}" build --remote "${image}" "${definition}"
        ;;
    plain)
        exec "${runtime}" build "${image}" "${definition}"
        ;;
    *)
        cat >&2 <<EOF
Usage: $0 [sudo|remote|plain]

Modes:
  sudo      Build on a local Linux machine where you have sudo privileges. Default.
  remote    Build with the Singularity remote builder.
  plain     Build directly, useful if your installation allows unprivileged builds.
EOF
        exit 2
        ;;
esac
