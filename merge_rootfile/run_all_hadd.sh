#!/bin/bash
set -euo pipefail

OUT=/eos/cms/store/group/phys_heavyions/zheng/pO_2026_May_26
mkdir -p "$OUT"   # use 'eos mkdir -p "$OUT"' if needed for EOS

# --- Data (list name and output name differ, so handle separately) ---
echo "=== Processing DATA ==="
./hadd_from_list.sh May_26_DATA.txt "$OUT/Data_May_26.root"

# --- MC: txt files (without .txt extension); output .root matches the name ---
SAMPLES=(
    MC_DYee_low_May26
    MC_DYee_May26
    MC_DYmu_low_May26
    MC_DYmu_May26
    MC_DYtau_low_May26
    MC_DYtau_May26
    MC_Wm_ele_May26
    MC_Wm_mu_May26
    MC_Wm_tau_May26
    MC_Wp_ele_May26
    MC_Wp_mu_May26
    MC_Wp_tau_May26
)

for s in "${SAMPLES[@]}"; do
    echo "=== Processing $s ==="
    ./hadd_from_list.sh "$s.txt" "$OUT/$s.root"
done

echo "All done."