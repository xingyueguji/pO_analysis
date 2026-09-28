#!/bin/bash
# plotting/run_combine_inputs.sh -- build the structured Combine inputs, with logs.
#
#   ./run_combine_inputs.sh [mu|ele|both]      (default: both)
#
# Runs the two macros that write the fit inputs, for the requested flavour(s):
#   mtandmet.C+(isElec)      -> plots[/Elec]/combine_input_W.root              (met, backup)
#                               plots[/Elec]/combine_input_W_leppt_mt40.root  (PRIMARY, + 6 CR dirs)
#   dileptonpeak.C+(isElec)  -> plots[/Elec]/combine_input_Z.root
# each with its `<input>_systs.txt` sidecar (the list of shape systematics
# actually written -- the fork's card generator READS these, never assumes).
#
# Why a wrapper and not `root -l -b -q 'mtandmet.C+(false)'` directly: these
# macros print numbers that get quoted downstream -- the ABCD QCD totals, the
# in-fit A0 = B0*C40/D0 renormalization and its rescale factor, the per-y QCD
# split, the region/systematic counts per output file -- and ROOT prints them to
# stdout ONLY, so a bare run discards the record. Same rule as run_qcd_abcd.sh /
# run_isolation.sh / run_syst_shapes.sh (see CLAUDE.md "When working on this
# repo"). Logs are named by THIS script so they stay consistent with every other
# stage: logs/{mtandmet,dileptonpeak}_<chan>.log.
#
# Order: run correction/run_qcd_abcd.sh FIRST (mtandmet embeds its QCD
# templates + abcd_counts), and re-run this after ANY re-skim.
#
# bash-3.2 safe (macOS stock /bin/bash), per the repo convention.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p logs

arg="${1:-both}"
case "$arg" in
  both) CHANS="mu ele" ;;
  mu|ele) CHANS="$arg" ;;
  *) echo "[ERR] usage: $0 [mu|ele|both]"; exit 2 ;;
esac

# isElec flag per channel (case-function, NOT an associative array: bash 3.2)
is_elec() { case "$1" in mu) echo false ;; ele) echo true ;; esac; }

# Pre-build both macros once, so two flavours cannot race on the ACLiC
# artifacts (skim/run_all.sh has been bitten by exactly that).
for m in mtandmet dileptonpeak; do
  root -l -b -q -e ".L ${m}.C+" >/dev/null 2>&1 || {
    echo "[ERR] ${m}.C does not compile:"
    root -l -b -q -e ".L ${m}.C+" 2>&1 | tail -20
    exit 1
  }
done

rc=0
for c in $CHANS; do
  e=$(is_elec "$c")
  for m in mtandmet dileptonpeak; do
    log="logs/${m}_${c}.log"
    echo "=== ${m}(${c})  ->  $log"
    if root -l -b -q "${m}.C+(${e})" >"$log" 2>&1; then
      # Headlines worth seeing without opening the log. "ERROR check plot N" is
      # benign macro debug output, not a failure -- excluded on purpose.
      grep -E "^\[QCD|^\[INFO\] Saved|^\[INFO\] shape systematics" "$log" | sed 's/^/    /' || true
      nw=$(grep -c "^\[WARN\]" "$log" || true)
      if [ "$nw" -gt 0 ]; then
        # The "<proc> is empty -> flooring central bin" warnings are expected
        # (empty wtau under the Z peak); report anything else in full.
        grep "^\[WARN\]" "$log" | grep -v "is empty -> flooring" | sed 's/^/    /' || true
        echo "    ($nw [WARN] lines total, $(grep -c 'is empty -> flooring' "$log" || true) of them the expected empty-template flooring)"
      fi
    else
      echo "[ERR] ${m}(${c}) failed -- tail of $log:"; tail -20 "$log"; rc=1
    fi
  done
done

if [ "$rc" -eq 0 ]; then
  echo "=== Combine inputs written (re-run correction/run_qcd_abcd.sh first if the skim changed):"
  ls -1 plots/combine_input_*.root plots/Elec/combine_input_*.root 2>/dev/null | sed 's/^/    /'
fi
exit $rc
