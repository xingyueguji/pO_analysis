# Sandbox: prefit r with pO charge-dependent electron SFs (2026-10-20)

A copy of the repo (no .git, backups, docs, plot outputs). The repo itself was not touched.

1. Reproduced the repo's prefit r here from the copied skim files: `prefit_r.csv`
   for leppt_mt40 and met_mt40 identical to the repo's.
2. `skim/po_charge_sf.h` + 4 edits in `skim/skim.C` (skim_Wel only): the nominal
   electron SF of the leading electron = pO reco SF(q, cell) x [iso ? pO ID+iso
   SF(q, cell) : 1] x the MB trigger SF, cells |eta_SC| < 1.6 / fwd <= -1.6 /
   bwd >= +1.6, values = the cut-and-count SFs of
   `correction/tnp_charge_ele.C` (reco: SC-energy probe), read from
   `skim/sf/tnp_charge_ele.csv`. The Up/Down twins keep their relative size.
3. Re-skimmed the 7 electron W MC samples (`skim/run_all.sh Wel ...`, logs in
   `skim/logs/`, the `[CHARGESF]` block = the SFs and <SF_EGM> vs <SF_pO> per
   cell); data unchanged. The EGM-SF skims are kept in `skim/rootfile_egm/`.
4. `correction/run_qcd_abcd.sh ele`, `plotting/run_combine_inputs.sh ele` (QCD
   and the ABCD subtraction follow the new MC; the LHE Up/Down twins are absent
   -- not needed for the prefit r).
5. `plotting/compare_chargesf.C` -> `plotting/plots/chargesf/{rpre,emu}_{Wp,Wm,W}`,
   `compare_chargesf.csv`, log `plotting/plots/compare_chargesf.log`.
