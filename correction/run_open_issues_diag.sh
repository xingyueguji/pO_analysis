#!/usr/bin/env bash
# correction/run_open_issues_diag.sh
#
# Run the diagnostics of the two open issues of the analysis note
# (AN-26-030, Section "Outstanding issues") and KEEP THE LOG:
#   (A) the charge-dependent electron/muon difference of the W yields,
#   (B) the Z->ee vs Z->mumu normalization behind the per-flavour r_Z.
#
# open_issues_diag.C reads only stored outputs (the Combine inputs, the Z skim
# histograms, the per-flavour cross-section CSV) -- no ntuple loop, a few
# seconds. ROOT prints to stdout only, hence this wrapper (same convention as
# run_qcd_abcd.sh). Re-run after run_combine_inputs.sh or run_observables.sh.
#
# Usage:  ./run_open_issues_diag.sh
# Output: logs/open_issues_diag.log   (THE record quoted by the AN)
#         plots/open_issues/            e/mu vs eta (W) and Z data/MC figures
#
# bash-3.2 safe (macOS stock /bin/bash).
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs

log="logs/open_issues_diag.log"
printf '[RUN ] open_issues_diag -> %s\n' "$log"
if root -l -b -q 'open_issues_diag.C+' > "$log" 2>&1; then
  status="OK"
else
  status="FAIL"
fi

# Headlines: the edge-bin e/mu rows, the double ratio, the Z ratios.
grep -E 'y0\+y11|double ratio|A_ch\(incl\) e - mu|e/mu of data/MC|mean\[80,100\)' "$log" | sed 's/^/       /' || true
printf '[%-4s] log: %s\n' "$status" "$log"
[ "$status" = "OK" ]
