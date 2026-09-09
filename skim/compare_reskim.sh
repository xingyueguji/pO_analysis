#!/bin/bash
# skim/compare_reskim.sh -- run compare_hists.C on every *_hist.root of a backup
# directory against the current rootfile/ (bit-identity gate after a re-skim).
#   ./compare_reskim.sh [backup-dir]     (default rootfile_pre_lhe)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
bk="${1:-rootfile_pre_lhe}"
[ -d "$bk" ] || { echo "[ERR] no backup dir $bk"; exit 2; }
root -l -b -q -e '.L compare_hists.C+' >/dev/null 2>&1 || true   # build once
bad=0
for f in "$bk"/*_hist.root; do
  b=$(basename "$f")
  [ -f "rootfile/$b" ] || { echo "[SUMMARY] $b: NOT in rootfile/ (skipped)"; continue; }
  out=$(root -l -b -q "compare_hists.C+(\"$f\",\"rootfile/$b\")" 2>&1 | grep -E "SUMMARY|DIFF|MISSING|ERR" || true)
  echo "$out" | sed "s#$bk/##; s# vs rootfile/[^:]*:#:#"
  # a file is bad when its SUMMARY counts DIFF/MISSING > 0 (the counts are
  # always printed, so do not grep for the words), or on an ERR line
  nd=$(echo "$out" | sed -n 's/.*DIFF \([0-9][0-9]*\)  MISSING \([0-9][0-9]*\).*/\1 \2/p')
  [ -n "$nd" ] || nd="1 1"
  set -- $nd
  { [ "$1" -ne 0 ] || [ "$2" -ne 0 ] || echo "$out" | grep -q "\[ERR\]"; } && bad=1
done
[ $bad -eq 0 ] && echo "ALL OLD HISTOGRAMS IDENTICAL" || echo "SOME FILES DIFFER -- see above"
exit $bad
