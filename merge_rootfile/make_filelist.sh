#!/bin/bash

# =========================
# Usage:
#   ./make_filelist.sh <input_dir> <output_txt>
#
# Example:
#   ./make_filelist.sh /eos/cms/store/group/phys_heavyions/anstahll/CERN/pO2025 files.txt
# =========================

INPUT_DIR="$1"
OUT_TXT="$2"

if [[ -z "$INPUT_DIR" || -z "$OUT_TXT" ]]; then
  echo "Usage: $0 <input_dir> <output_txt>"
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "ERROR: directory does not exist: $INPUT_DIR"
  exit 1
fi

echo "Scanning directory:"
echo "  $INPUT_DIR"
echo "Writing file list to:"
echo "  $OUT_TXT"

# Find all ROOT files recursively, sort for reproducibility
find "$INPUT_DIR" -type f -name "*.root" | sort > "$OUT_TXT"

NFILES=$(wc -l < "$OUT_TXT")
echo "Done. Found $NFILES ROOT files."