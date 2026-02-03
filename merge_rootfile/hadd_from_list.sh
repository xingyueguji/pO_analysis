#!/bin/bash

# =========================
# Usage:
#   ./hadd_from_list.sh <filelist.txt> <output_root>
#
# Example:
#   ./hadd_from_list.sh files.txt /eos/cms/store/group/phys_heavyions/zheng/pO_2025.root
# =========================

FILELIST="$1"
OUTROOT="$2"

if [[ -z "$FILELIST" || -z "$OUTROOT" ]]; then
  echo "Usage: $0 <filelist.txt> <output_root>"
  exit 1
fi

if [[ ! -f "$FILELIST" ]]; then
  echo "ERROR: file list not found: $FILELIST"
  exit 1
fi

OUTDIR=$(dirname "$OUTROOT")

if [[ ! -d "$OUTDIR" ]]; then
  echo "ERROR: output directory does not exist: $OUTDIR"
  exit 1
fi

NFILES=$(wc -l < "$FILELIST")

if [[ "$NFILES" -eq 0 ]]; then
  echo "ERROR: file list is empty"
  exit 1
fi

echo "Merging $NFILES ROOT files into:"
echo "  $OUTROOT"

# Use xargs to avoid command-line length limits
xargs -a "$FILELIST" hadd -f -j 4 "$OUTROOT"

echo "Done."