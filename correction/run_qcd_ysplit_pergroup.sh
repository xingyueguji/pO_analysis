#!/usr/bin/env bash
# correction/run_qcd_ysplit_pergroup.sh
#
# BEFORE vs AFTER of the rapidity-resolved in-fit ABCD (QCD_MODE=abcdy,
# 2026-10-19) and KEEP THE LOG: per rapidity group (fwd y10-11 | bwd y0-1 |
# cen y2-9, correction/abcd_ygroups.h), flavour x charge, lab and fb, the SR
# multijet of the leppt_mt40 fit with one ABCD per charge (qcd_abcd) vs one
# ABCD per group (qcd_abcdy) with the group ABCD stat error, the group transfer
# factors T_g, R_g, and the e/mu of the prefit r_Q = (D-B-Q)/W with both.
#
# qcd_ysplit_pergroup.C reads only stored outputs (the leppt_mt40 Combine
# inputs, rootfile/qcd_abcd_{mu,ele}.root) -- a few seconds. ROOT prints to
# stdout only, hence this wrapper (same convention as run_qcd_ysplit.sh).
# Order: ./run_qcd_abcd.sh, then ../plotting/run_combine_inputs.sh, then this.
#
# Usage:  ./run_qcd_ysplit_pergroup.sh
# Output: logs/qcd_ysplit_pergroup.log      (THE record)
#         plots/qcd_ysplit_pergroup/        figures + qcd_ysplit_pergroup.csv
#
# bash-3.2 safe (macOS stock /bin/bash).
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs

log="logs/qcd_ysplit_pergroup.log"
printf '[RUN ] qcd_ysplit_pergroup -> %s\n' "$log"
if root -l -b -q "qcd_ysplit_pergroup.C+" > "$log" 2>&1; then
  status="OK"
else
  status="FAIL"
fi
# root exits 0 even when a macro bails: any [ERROR] line or a failed
# self-check is a failure too
grep -q '^\[ERROR\]' "$log" && status="FAIL"
grep -q 'MISMATCH' "$log" && status="FAIL"

# Headlines: the self-checks, the per-group table of every set, and the
# summed rows of the e/mu tables (edge, backward, groups, inclusive).
grep -E '^--- |\[check\]|^  (fwd|bwd|cen)  |^  sum  |^=== e/mu|^  W[pm]  bin|^ +(y0\+y11|y0\.\.y3|fwd|bwd|cen|incl) \|' "$log" \
  | sed 's/^/       /' || true
printf '[%-4s] log: %s\n' "$status" "$log"
[ "$status" = "OK" ]
