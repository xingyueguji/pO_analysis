#!/usr/bin/env bash
# correction/run_recoil_res.sh
#
# Measure the MC MET-resolution factor on the Z recoil (u_perp width, data vs
# DY MC, per flavour) and KEEP THE LOG (2026-10-20). recoil_res.C reads only
# the Z skim histograms -- a few seconds -- and writes rootfile/recoil_res.root,
# the input of skim/met_res.h: run this BEFORE the W MC skim (the skim is
# FATAL without it) and again after any Z re-skim.
#
# Usage:  ./run_recoil_res.sh
# Output: logs/recoil_res.log     (THE record)
#         rootfile/recoil_res.root, plots/recoil_res/
#
# bash-3.2 safe (macOS stock /bin/bash).
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs

log="logs/recoil_res.log"
printf '[RUN ] recoil_res -> %s\n' "$log"
if root -l -b -q 'recoil_res.C+' > "$log" 2>&1 && grep -q '^\[DONE\]' "$log"; then
  status="OK"
else
  status="FAIL"
fi

grep -E '^\s+(mu |ele): k =|k = sigma_data|closure:' "$log" | sed 's/^/       /' || true
printf '[%-4s] log: %s\n' "$status" "$log"
[ "$status" = "OK" ]
