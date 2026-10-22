#!/usr/bin/env bash
# correction/run_qcd_ysplit.sh
#
# Test of the rapidity split of the W signal-region multijet template (step 1
# of the charge-split follow-up, AN-26-030 sec:openissues:charge) and KEEP THE
# LOG: per lab/FB rapidity bin and charge, the fit's qcd_abcd template (the
# y-inclusive ABCD total split like the low-MET excess) vs the per-y ABCD
# estimate B_y*C40_y/D_y from the y-binned (iso, m_T) planes of the W skim,
# and the per-y transfer factors T_y = B_y/D_y, R_y = C40_y/D_y.
#
# qcd_ysplit.C reads only stored outputs (the W skims, the leppt_mt40 Combine
# inputs, qcd_abcd_{mu,ele}.root for a cross-check) -- a few seconds. ROOT
# prints to stdout only, hence this wrapper (same convention as run_qcd_abcd.sh).
# Needs the W skims written since 2026-10-19 (h_iso_mt_*_y<i>[_FB]); re-run
# after a W re-skim or run_combine_inputs.sh.
#
# Usage:  ./run_qcd_ysplit.sh [r_W]      (default 1.2 = the fitted W scale)
# Output: logs/qcd_ysplit.log            (THE record)
#         plots/qcd_ysplit/              figures + qcd_ysplit.csv
#
# bash-3.2 safe (macOS stock /bin/bash).
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs

rW="${1:-1.2}"
log="logs/qcd_ysplit.log"
printf '[RUN ] qcd_ysplit (r_W = %s) -> %s\n' "$rW" "$log"
if root -l -b -q "qcd_ysplit.C+(${rW})" > "$log" 2>&1; then
  status="OK"
else
  status="FAIL"
fi
# root exits 0 even when a macro bails: any [ERROR] line is a failure too
grep -q '^\[ERROR\]' "$log" && status="FAIL"

# Headlines: the closure lines, the edge-bin (y0+y11) rows, the chi2 lines and
# the edge rows of the e/mu table (template vs per-y ABCD split).
grep -E '^--- |\[closure\]|y0\+y11: abs|chi2 vs 1|^=== e/mu|^  W[pm]  bin|^    y0\+y11 \||^     y0\.\.y3 \||^      incl \||\[sum\]' "$log" \
  | sed 's/^/       /' || true
printf '[%-4s] log: %s\n' "$status" "$log"
[ "$status" = "OK" ]
