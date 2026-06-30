#!/usr/bin/env bash
#
# setup_ncbi.sh - one-time NCBI SRA configuration for the plasmid discovery pipeline.
#
# This automates the bits that otherwise have to be remembered and typed by hand:
#   1. SRA Toolkit cache/temp location  (vdb-config /repository/user/main/public/root)
#   2. NCBI API key for sra-tools        (vdb-config /NCBI/api_key)
#   3. NCBI API key for EDirect          (export NCBI_API_KEY in your shell profile)
#
# Why both #2 and #3: prefetch/fasterq-dump read the key from vdb-config, while
# esearch/efetch (EDirect) read it from the NCBI_API_KEY environment variable.
#
# An API key raises the E-utilities rate limit from 3 to 10 requests/second and
# is free from https://www.ncbi.nlm.nih.gov/account/settings/ ("API Key Management").
#
# Usage:
#   ./setup_ncbi.sh                         # interactive prompts
#   ./setup_ncbi.sh -k <API_KEY> -c <DIR>   # non-interactive
#
# Re-running is safe; it overwrites the same settings.

set -euo pipefail

API_KEY=""
CACHE_DIR=""
PROFILE="${HOME}/.bashrc"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -k, --api-key KEY    NCBI API key (default: prompt)
  -c, --cache-dir DIR  SRA Toolkit cache/download directory (default: prompt)
  -p, --profile FILE   Shell profile to append NCBI_API_KEY to (default: ~/.bashrc)
  -h, --help           Show this help

With no -k/-c the script prompts interactively. Press Enter at a prompt to skip
that setting (e.g. leave the cache dir unchanged).
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -k|--api-key)   API_KEY="$2"; shift 2 ;;
        -c|--cache-dir) CACHE_DIR="$2"; shift 2 ;;
        -p|--profile)   PROFILE="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if ! command -v vdb-config >/dev/null 2>&1; then
    echo "ERROR: vdb-config not found. Activate the conda env first:" >&2
    echo "       mamba env create -f environment.yml && conda activate plasmidextractor" >&2
    exit 1
fi

# Prompt for anything not supplied on the command line.
if [[ -z "$API_KEY" ]]; then
    read -r -p "NCBI API key (Enter to skip): " API_KEY
fi
if [[ -z "$CACHE_DIR" ]]; then
    read -r -p "SRA cache/download directory [e.g. /data/sra-cache] (Enter to skip): " CACHE_DIR
fi

# 1 + 2. SRA Toolkit configuration via vdb-config.
if [[ -n "$CACHE_DIR" ]]; then
    mkdir -p "$CACHE_DIR"
    vdb-config --set "/repository/user/main/public/root=${CACHE_DIR}"
    echo "[ok] SRA cache root set to: $CACHE_DIR"
fi

if [[ -n "$API_KEY" ]]; then
    vdb-config --set "/NCBI/api_key=${API_KEY}"
    echo "[ok] NCBI API key written to vdb-config (used by prefetch/fasterq-dump)"

    # 3. Make the same key available to EDirect via the environment.
    if grep -q "NCBI_API_KEY" "$PROFILE" 2>/dev/null; then
        echo "[skip] NCBI_API_KEY already present in $PROFILE - not duplicating"
    else
        printf '\n# NCBI API key for EDirect (esearch/efetch)\nexport NCBI_API_KEY=%s\n' "$API_KEY" >> "$PROFILE"
        echo "[ok] export NCBI_API_KEY appended to $PROFILE"
        echo "     run 'source $PROFILE' or open a new shell to load it"
    fi
fi

echo
echo "Done. Verify with:  vdb-config --cfg | grep -E 'api_key|root'"
