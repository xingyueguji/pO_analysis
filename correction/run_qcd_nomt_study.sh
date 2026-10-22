#!/usr/bin/env bash
# correction/run_qcd_nomt_study.sh
#
# The multijet study of the NO-m_T lepton-pT fit (2026-10-21, user decision to
# make the no-m_T pT fit the nominal) and KEEP THE LOG: per flavour and charge,
# the composition of the no-m_T signal region, the pT-dependent transfer factor
# F(pT) at m_T < 30 per rapidity group, the per-bin anti-iso shapes, the
# m_T factorization, the expected precision (Asimov Fisher matrix) and the data
# stand-in fits of the candidate in-fit schemes A / C / R, and the e/mu of the
# fitted r -- every m_T item with the PF MET and with the charged-only MET.
#
# qcd_nomt_study.C reads only the W skims (the y-binned (iso, pT, m_T) cubes
# written by the 2026-10-21 skim.C) -- under a minute. ROOT prints to stdout
# only, hence this wrapper (same convention as run_qcd_abcd.sh).
#
# Usage:  ./run_qcd_nomt_study.sh
# Output: logs/qcd_nomt_study.log        (THE record)
#         plots/qcd_nomt_study/          figures + qcd_nomt_study.csv
#
# bash-3.2 safe (macOS stock /bin/bash).
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs

log="logs/qcd_nomt_study.log"
printf '[RUN ] qcd_nomt_study -> %s\n' "$log"
if root -l -b -q 'qcd_nomt_study.C+' > "$log" 2>&1; then
  status="OK"
else
  status="FAIL"
fi
# root exits 0 even when a macro bails: any [ERROR] line is a failure too
grep -q '^\[ERROR\]' "$log" && status="FAIL"

# Headlines: the slopes and group tests, the SR multijet m_T split, the
# expected <r> precision per scheme, the fitted T_g closure, the e/mu edges.
grep -E '^=============== |\[FPT\]   (incl|fwd |bwd |cen ):|group dependence|SR multijet|\[FISHER\]   |T_g fitted|^\[EMU\]' "$log" \
  | sed 's/^/       /' || true
printf '[%-4s] log: %s\n' "$status" "$log"
[ "$status" = "OK" ]
