#!/usr/bin/env bash
# skim/run_gen_xsec.sh -- the GEN fiducial cross sections, with a log (2026-10-13).
#
#   ./run_gen_xsec.sh
#
# Runs gen_xsec.C over all generated events of the four W and two DY MC files
# (paths from ResolveMCSample in skim_common.h) and writes
#   rootfile/gen_xsec.root          h_gen_sig_{Wp,Wm}[_mt40][_FB][_epps21], h_gen_sig_Z[_epps21], ...
#   output/gen_xsec_fid.txt         per-bin sigma_gen, NO m_T cut      (met fit)
#   output/gen_xsec_fid_mt40.txt    per-bin sigma_gen, m_T > 40        (leppt_mt40 fit)
# The two sidecars are what the fork's run_pO_fits.sh --contour bakes into the
# sigma_W-as-POI workspace: re-upload them (sync_lxplus.sh upload) after a re-run.
# Log: logs/gen_xsec.log = THE record of the quoted sigma_fid numbers.
# Discriminant-independent and reco-independent: re-run only when the MC
# production, the fiducial definition or mc_norm.h's sigma changes. ~10 min.
# bash-3.2 safe (macOS /bin/bash).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p rootfile output logs
log="logs/gen_xsec.log"

echo "[run_gen_xsec] gen_xsec.C -> $log"
# ROOT exits 0 even when a macro bails, so judge by the log and the outputs
root -l -b -q 'gen_xsec.C+' >"$log" 2>&1 || { echo "[FAIL] root returned non-zero -- see $log"; exit 1; }
if grep -q '^\[ERROR\]\|\[FATAL\]\|Error in <' "$log"; then
  echo "[FAIL] errors in $log:"; grep '^\[ERROR\]\|\[FATAL\]\|Error in <' "$log" | head -20; exit 1
fi
for f in rootfile/gen_xsec.root output/gen_xsec_fid.txt output/gen_xsec_fid_mt40.txt; do
  [ -s "$f" ] || { echo "[FAIL] $f not written -- see $log"; exit 1; }
done
# the headlines: sigma_fid per set, the Z, the member-0 closure and any warning
grep 'sigma_fid\|member0\|\[WARN\]' "$log" || true
echo "[OK] rootfile/gen_xsec.root, output/gen_xsec_fid.txt, output/gen_xsec_fid_mt40.txt (log: $log)"
