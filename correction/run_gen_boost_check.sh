#!/usr/bin/env bash
# correction/run_gen_boost_check.sh
#
# Is the W/Z hard process of the MC boosted to the pO lab frame? Run this on
# every new MC production BEFORE switching to it (the Aug-20 production fails:
# CLAUDE.md, TODO "MC rapidity boost").
#
# gen_boost_check.C prints the record -- the Born-event beam energies, the
# rapidity of the NN CM frame in the MC record, the oxygen side in data (ZDC)
# and MC (gen record), the data shift scan -- and one [VERDICT] line per
# check. ROOT writes all of it to stdout only, hence this wrapper (same
# convention as run_charge_flip.sh): one log per mode under correction/logs/.
#
# Usage:  ./run_gen_boost_check.sh [gen|side|scan|all]      (default: all)
#   gen   generator level, MC only (~2 min)            -- the check itself
#   side  + beam direction in data and MC (~4 min)
#   scan  + data vs shifted MC on the current leppt_mt40 Combine inputs;
#         meaningful only once those inputs are built from the MC under test
#   all   the three
#
# A NEW production, before repointing skim/skim_common.h:
#   MC_WP=/path/Wp_mu.root MC_WM=/path/Wm_mu.root MC_DY=/path/DY_mu_Z.root \
#     ./run_gen_boost_check.sh gen
# (unset = pOSkim::ResolveMCSample, the production in use; DATA_FILE overrides
#  pOSkim::kDefaultDataFile for mode side)
#
# Outputs: logs/gen_boost_check_<mode>.log   full console output (THE record)
#          rootfile/gen_boost_check.root     histograms + scan graphs
#          plots/gen_boost/                  born_beam_energy, delta, delta_vs_r,
#                                            data_shift_scan
#
# bash-3.2 safe (macOS stock /bin/bash): no associative arrays, no mapfile.
set -euo pipefail

cd "$(dirname "$0")"

MODE="${1:-all}"
case "$MODE" in
  gen|side|scan|all) ;;
  *) echo "usage: $0 [gen|side|scan|all]" >&2; exit 2 ;;
esac

mkdir -p logs rootfile plots

echo "[BUILD] compiling gen_boost_check.C ..."
if ! root -l -b -q -e '.L gen_boost_check.C+' > logs/gen_boost_check_build.log 2>&1; then
  echo "[FAIL] compilation failed -- see logs/gen_boost_check_build.log" >&2
  exit 1
fi

log="logs/gen_boost_check_${MODE}.log"
args="\"${MODE}\",\"${MC_WP:-}\",\"${MC_WM:-}\",\"${MC_DY:-}\",\"${DATA_FILE:-}\""
printf '[RUN ] %s -> %s\n' "$MODE" "$log"
SECONDS=0
if root -l -b -q "gen_boost_check.C+(${args})" > "$log" 2>&1 && grep -q '^\[DONE\]' "$log"; then
  status="OK"
else
  status="FAIL"
fi

grep -E '^\[CONFIG\]|^\[INPUT\]|^\[RESULT\]|^\[VERDICT\]|^\[WARN\]|^\[ERR|^\[FATAL' "$log" 2>/dev/null | sed 's/^/       /' || true
printf '[%-4s] %s  %ds   log: %s\n' "$status" "$MODE" "$SECONDS" "$log"
echo "============================================================"
if [ "$status" = "OK" ]; then
  echo "[DONE] ran OK -- the [VERDICT] lines above are the answer"
else
  echo "[DONE] the macro failed -- check $log" >&2
  exit 1
fi
