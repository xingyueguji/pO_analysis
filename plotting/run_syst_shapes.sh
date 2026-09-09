#!/bin/bash
# plotting/run_syst_shapes.sh -- LHE shape-systematics diagnostics with a log.
#
#   ./run_syst_shapes.sh [met|leppt_mt40|all]      (default: leppt_mt40)
#
# Runs syst_shapes.C for the discriminant(s): per-region Up/Down-over-nominal
# plots of nPDF / qcdScale / alphaS for every MC process of the Combine inputs,
# the per-region overlays of the nominal signal with all 106 EPPS21 member
# templates (members/), the integral-shift summaries vs rapidity bin, and the INCLUSIVE CONSISTENCY
# tables (member-level inclusive from the skim twins vs the all-events
# reference in skim/output/lhe_weights.txt; sum-of-per-bin Up/Down = what one
# collapsed nuisance implies; LHAPDF-vs-pOLhe::Hessian per-bin closure).
# ROOT prints the tables to stdout only -> logs/syst_shapes_<disc>.log is the
# record. bash-3.2 safe.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p logs
arg="${1:-leppt_mt40}"
case "$arg" in
  all) DISCS="met leppt_mt40" ;;
  met|leppt_mt40) DISCS="$arg" ;;
  *) echo "[ERR] usage: $0 [met|leppt_mt40|all]"; exit 2 ;;
esac
# pre-build once (no ACLiC race if several discs run)
root -l -b -q -e '.L syst_shapes.C+' >/dev/null 2>&1 || { echo "[ERR] syst_shapes.C does not compile"; root -l -b -q -e '.L syst_shapes.C+' 2>&1 | tail -20; exit 1; }
rc=0
for d in $DISCS; do
  log="logs/syst_shapes_${d}.log"
  echo "=== syst_shapes($d)  ->  $log"
  if root -l -b -q "syst_shapes.C+(\"$d\")" >"$log" 2>&1; then
    grep -E "^\[INCL\]|^\[CLOSURE\]|^\[SUMMARY\]|^\[WARN\]|^\[ERR|^\[OK\]" "$log" || true
  else
    echo "[ERR] syst_shapes($d) failed -- tail of $log:"; tail -20 "$log"; rc=1
  fi
done
exit $rc
