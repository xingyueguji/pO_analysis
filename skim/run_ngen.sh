#!/usr/bin/env bash
# skim/run_ngen.sh
#
# Compute N_gen (sum of the per-event generator weight over ALL events, no
# selection) for every MC sample, and save it to:
#   output/ngen.txt    -- human-readable table
#   rootfile/ngen.root -- h_ngen / h_nraw + per-sample TParameter<double>
#
# This is the generator-level denominator for the cross-section normalization.
# It is orthogonal to run_all.sh (which skims with selection); the MC file paths
# come from the same single source of truth (ResolveMCSample in skim_common.h).
#
# Usage:
#   ./run_ngen.sh
#
# Run it once after the MC files are in place; re-run only if the MC changes.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p rootfile output logs
log="logs/count_ngen.log"

echo "============================================================"
echo "Counting N_gen = Sum(HiTree::weight) over all MC samples ..."
echo "============================================================"

if root -l -q -b 'count_ngen.C+' >"$log" 2>&1; then
  echo "[OK] wrote output/ngen.txt and rootfile/ngen.root   (log: $log)"
  # Echo the resulting table so the run is self-documenting.
  [[ -f output/ngen.txt ]] && { echo; cat output/ngen.txt; }
  exit 0
else
  echo "[FAIL] count_ngen.C returned non-zero -- see $log"
  exit 1
fi
