#!/bin/bash

# =========================
# Usage:
#   ./hadd_from_list.sh <filelist.txt> <output_root> [njobs]
#
# Example:
#   ./hadd_from_list.sh files.txt /eos/cms/store/group/phys_heavyions/zheng/pO_2026_May_26/out.root
#   ./hadd_from_list.sh files.txt /eos/.../out.root 8     # use 8 parallel jobs
# =========================

FILELIST="$1"
OUTROOT="$2"
NJOBS="${3:-4}"   # number of parallel hadd jobs (-j), default 4

if [[ -z "$FILELIST" || -z "$OUTROOT" ]]; then
  echo "Usage: $0 <filelist.txt> <output_root> [njobs]"
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

# strip CR (Windows line endings) and blank lines so the count and the merge agree
CLEANLIST=$(mktemp)
sed 's/\r$//' "$FILELIST" | grep -v '^[[:space:]]*$' > "$CLEANLIST"

NFILES=$(wc -l < "$CLEANLIST")
if [[ "$NFILES" -eq 0 ]]; then
  echo "ERROR: file list is empty"
  rm -f "$CLEANLIST"
  exit 1
fi

echo "==================================================="
echo "Input list : $FILELIST  ($NFILES files)"
echo "Output     : $OUTROOT"
echo "Parallel   : -j $NJOBS"
echo "==================================================="

START=$(date +%s)

# -f : overwrite existing output
# -j : run NJOBS merge subprocesses in parallel
# @  : read input paths from the file (no command-line length limit)
hadd -f -j "$NJOBS" "$OUTROOT" "@$CLEANLIST"
STATUS=$?

END=$(date +%s)
rm -f "$CLEANLIST"

if [[ $STATUS -ne 0 ]]; then
  echo "ERROR: hadd exited with status $STATUS"
  exit $STATUS
fi

echo "Done. Merged $NFILES files -> $OUTROOT  (${SECONDS}s wall, hadd $((END-START))s)"