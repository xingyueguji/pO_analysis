# pO_analysis — W/Z boson cross-section measurement in pO collisions

## What this repo is

A CMS-style analysis extracting W and Z boson cross-sections and forward/backward
ratios from proton-oxygen (pO) collision data. Four channels are measured:
W→μν, W→eν, Z→μμ, Z→ee. Tau channels (W→τν, Z→ττ) are treated as backgrounds
and the infrastructure exists for them but is not yet wired up end-to-end.

This repo produces histograms. The actual fits run in a separate checkout of
HiggsAnalysis-CombinedLimit (see "Downstream fit" below).

The code is ROOT C++ macros (not python / coffea / RDataFrame). Everything is
driven by bash, with no Snakemake / Make / config files — parameters are
hardcoded inside the macros.

## Pipeline (three stages)

```
EOS NanoAOD-like ntuple (ggHiNtuplizer/EventTree, ~2.2M events)
        │
        ▼
[skim/]      Event selection + cutflow → per-sample ROOT histograms
        │
        ▼
[plotting/]  Data/MC overlays, QCD sideband fit, ROC curves → publication plots
        │
        ▼
[analysis/]  Charge asymmetry & forward/backward ratio extraction → fit inputs
        │
        ▼
HiggsCombine datacards + workspaces → cross-section fit
```

### Stage 1 — `skim/`

Entry point: [skim/run_all.sh](skim/run_all.sh)

```bash
./run_all.sh <channel> [samples]
#   channel: Zmm | Zee | Wmu | Wel | all
#   samples: Data DY Wp Wm DYtau Wptau Wmtau   (default: all 7)
```

Driver invokes ROOT with the unified dispatcher:
```
root -l -q -b 'skim.C+(kWmu, "<DATA_FILE>", kWp)'
```

Layout (post-refactor):
- [skim/skim_common.h](skim/skim_common.h) — shared utilities, enums (`SampleType`, `ChannelType`), kinematics, cutflow helpers, pO event filters, etc. Wrapped in `namespace pOSkim`.
- [skim/skim.C](skim/skim.C) — one file, four entry functions plus a dispatcher:
  - `skim_Wmu(file, sample)` — W→μν
  - `skim_Wel(file, sample)` — W→eν
  - `skim_Zmm(file, sample)` — Z→μμ
  - `skim_Zee(file, sample)` — Z→ee
  - `skim(channel, file, sample)` — dispatcher used by `run_all.sh`
- [skim/legacy/](skim/legacy/) — the original four per-channel macros, preserved verbatim for diffing/fallback. The unified files were byte-identical in physics output to these until 2026-07-06, when the relIso definition gained the Δβ PU correction (legacy keeps the uncorrected sum — expect small selection diffs).
- [skim/count_ngen.C](skim/count_ngen.C) + [skim/run_ngen.sh](skim/run_ngen.sh) — compute N_gen (Σ gen weight over **all** events, no selection) per MC sample → `skim/rootfile/ngen.root` + `skim/output/ngen.txt`. The cross-section-normalization denominator. See "MC normalization" below.
- [skim/gen_xsec.C](skim/gen_xsec.C) — **GENERATOR-level W σ per rapidity bin (2026-08-05; FIDUCIAL since 2026-08-12)**: loops ALL generated events of the 4 W MC files (no reco, no selection; needs the `vector<vector<int>>` dictionary like the skim), picks the gen charged lepton of the sample's flavour+charge (highest-pT, W-ancestor preferred; falls back with a WARN — the ntuple mother chain often omits the W), histograms its LAB η in `kYEdges`(+FB) with the POWHEG weight → `skim/rootfile/gen_xsec.root`. **Binned in `y = −η_lab`, the skim's convention** ("p-going (−Z) = forward", `skim.C:835,1497`) — **BUG FIXED 2026-08-12**: it previously filled raw +η, so gen bin i paired with the MIRRORED reco region y_i. Caught by the per-bin (A·ε) diagnostic (mirrored pairing made A·ε charge-dependent — μ W⁻ ran 1.27→0.73 across η — instead of the charge-symmetric detector response it must be; reversing the gen bins collapsed the W⁺/W⁻ spread from 0.54 to 0.06). Effect: per-bin σ was mirrored, inclusive σ shifted ≤2% (W⁻ met 40.15→39.68, leppt_mt40 44.86→43.89); counts/charge_asym/FBratio never read this file and were unaffected. **`h_gen_sig_{Wp,Wm}[_FB]` are FIDUCIAL since 2026-08-12**: gen (bare, post-FSR) lepton pT > 25 (`kFidPtMin` = the skim's nominal cut; |η_lab| < 2.4 implicit in the binning window; NO m_T cut — one fiducial serves all discs), so **σ_meas,i = r_i × σ_gen-fid,i directly** (= yield-based extraction with MC A·ε, algebraically identical; kA_O and kSigma cancel between r and σ_gen — see the r×σ_gen plan under "Future plans"). `h_gen_tot_{Wp,Wm}[_FB]` keep the old no-pT-cut definition — but **NOT a true total**: the ntuple gen collection `HiGenParticleAna/hi` is itself filtered at **pT > 5 GeV and |η| < 2.5** (measured 2026-08-12), so fid/tot = 0.675 (W⁺) / 0.755 (W⁻) is a pT>25/pT>5 ratio, not an acceptance. **σ_fid is unaffected** (the whole fiducial sits inside the filter and Σw runs over all events) — the filter is also why gen_xsec.C WARNs that 25–37% of events have "no gen lepton" (those are below 5 GeV or beyond |η| 2.5, outside the fiducial either way). Normalization is the weighted FRACTION × `mc_norm.h` σ: σ_i = kA_O·kSigma·Σw_i/Σw_all — **unit-proof, because the July-29 `weight` is σ in pb (⟨w⟩≈6376) while `kSigma_*` is nb**; the pipeline was always consistent since σ cancels in k_s = A·σ·L/Σw. Discriminant-independent one-time producer consumed by `xsec_fiducial_comb` (gen overlay). Gen storage window extends past |η|=2.4 (nonzero over/underflow → edge bins safe); per-flavour σ: fid(pT>25, |η|<2.4) W⁺ 49.8 / W⁻ 39.6 nb; in-window no-pT-cut 73.7/52.5; total 102.0/87.4.
- [skim/mc_norm.h](skim/mc_norm.h) — single source of truth for the absolute MC→data scale `k_s = σ·L/N_gen` (`pONorm::MCScale`); consumes `ngen.root`. See "MC normalization" below.

Isolation studies (ROC curves, working-point tuning) are per-flavour and live in
`correction/` (moved out of `skim/`):
- [correction/isolation_mu_tight.C](correction/isolation_mu_tight.C) — muon (current: pure TightID Δβ scan)
- [correction/isolation_ele.C](correction/isolation_ele.C) — electron
- [correction/isolation.C](correction/isolation.C) — muon legacy multi-cone/multi-def scan (kept untouched)
- plotted by [correction/plot_iso_summary.C](correction/plot_iso_summary.C) / [correction/PlotsIsoROC.C](correction/PlotsIsoROC.C) / [correction/PlotIsoROC_ele.C](correction/PlotIsoROC_ele.C). See the `correction/` section below.

**relIso definition (since 2026-07-06):** every isolation cut in the pipeline uses
the **Δβ PU-corrected** PF relIso (MuonPOG convention),
`(ch + max(0, neu + pho − 0.5·PU)) / pT`, computed from the
`{mu,ele}PF{Ch,Neu,Pho,PU}Iso` branches — `pOSkim::RelIsoPF` in `skim_common.h`
is the single source (5-arg; null PU vector ⇒ uncorrected fallback + `[WARN]`).
Verified to reproduce the ntuplizer's `elePFRelIsoWithDBeta` exactly. Cut values
unchanged (μ 0.15 / e 0.095 — both re-confirmed as Youden-J optima under Δβ).

**Selection (8-step cutflow in the W macros)**:
1. ≥1 PF lepton with pT > 25 GeV
2. pO event-quality filters (auto-detected from branches)
3. Trigger: `HLT_OxyL1SingleMuOpen_v1` (muon) or `HLT_OxyL1SingleEG10_v1` (electron) — L1-seeded paths (skim.C:362, 1048; the older "HLT_PAL3Mu12/PAL3Ele12" in this doc was WRONG, corrected 2026-08-12)
4. DY veto: reject events with an OS dilepton pair, both legs ID'd + isolated
   (μ: pT>15, Tight, relIso<0.15; e: pT>10, eleMVAIdWP95, relIso<0.095), mll in (80,110)
5. ≥1 Tight ID lepton
6. Leading lepton pT > 25 GeV
7. Leading lepton PF relIso < 0.15
8. Trigger matching on leading lepton (ΔR < `trigMatchDR` = **0.4 since 2026-08-12**, was 0.1; skim.C:283/963, against the objects in `hltobject/HLT_Oxy*` = **L1-granularity candidates**). **Only the two W channels do this** — `skim_Zmm`/`skim_Zee` wire the trigger-object branches but never match on them.

**FIXED 2026-08-12: `trigMatchDR` 0.1 → 0.4 in both W channels + `correction/njet_WZ.C:194` (which replicates the selection), W re-skim and full downstream re-run DONE.** Post-fix cutflow step 7→8: μ⁺/μ⁻ **99.99%/99.99%** (was 93.88%/98.82%), e⁺/e⁻ **99.89%/99.88%** (was 97.79%/98.89%) — charge-symmetric in both flavours. **Data W candidates: μ 5203 → 5480 (+5.3%), e 7065 → 7186 (+1.7%)**; `njet_WZ` data counts still equal the skim's N[8] exactly (5480 / 7186), so the two independent selection implementations remain event-identical. ABCD moved with it: iso-pass QCD μ⁺/μ⁻ 413.7/453.8 → **491.6/473.2** (the artificial μ⁺ deficit is gone from the sideband too — ratio 0.91 → 1.04), e⁺/e⁻ 1887.2/1693.9 → **1937.9/1723.1**; MET-plane T μ 0.1895/0.1795, e 0.3704/0.3204. `ptmt_scan` baseline S/√(S+B) μ 53.1 → **54.40**, e 39.6 → **39.85** (same conclusions: optimum still at the pT floor with m_T > 40; e best (20, 50) = 50.27). **The defect it removes:** The saved trigger object's φ is at the muon station, `muPhi` is the vertex-direction φ, and the difference decomposes **exactly** (measured in barrel pT slices) as

&nbsp;&nbsp;&nbsp;&nbsp;`Δφ(q, pT) = c + q·k/pT`, with **k = 2.55 rad·GeV** (solenoid bending out to R ≈ 4.5 m = the muon station; the fitted k is constant to 3% over pT 25–200) and a **charge-independent offset c = +0.0060 rad** (constant to 5% over the same range).

Both terms are understood. `k`: `muPhi` is the momentum direction **at the vertex**, the object φ is the muon's azimuth **at the muon station**, and in the 3.8 T solenoid the track's azimuth changes continuously between the two — the bend *direction* is set by F = qv×B, hence the ±q. `c`: the saved objects are **L1-granularity candidates** (pT quantized in 0.5 GeV steps; φ takes exactly **576 distinct values** = the μGMT grid 2π/576 = 10.9 mrad; η likewise on a ~11 mrad grid) — reporting the cell's lower edge instead of its centre biases φ by half a cell = **5.45 mrad**, matching the measured c = 6.0 mrad.

The bending term alone is charge-*antisymmetric* and would cost both charges equally — **`c` is what breaks the symmetry**: it adds to the μ⁺ bending and cancels part of the μ⁻ one, so |Δφ⁺| − |Δφ⁻| = 2c ≈ **+0.012 rad at every pT**. That 12 mrad is decisive only because ΔR < 0.1 slices straight through the bending scale at low pT (Δφ ≈ 0.09–0.10 rad at pT 25–30 in the barrel): matching efficiency by pT slice, μ⁺/μ⁻ = **0.478/0.905** (25–30), 0.956/0.988 (30–35), ≥0.985/0.993 above 35. The lepton pT spectra are NOT the cause (⟨1/pT⟩ differs by only 2% between charges). Net: ~10% charge-dependent barrel-only inefficiency, **entirely concentrated at pT 25–30** — i.e. right on the Jacobian turn-on, the worst place for the lepton-pT discriminant. Endcap is unaffected (Δφ ≈ 0.04, both charges ≥0.998). Cutflow confirms it inclusively (step 7→8: μ⁺ −5.72%, μ⁻ −1.11%; every other step agrees between charges to <0.2%), and a direct per-cut MC scan shows reco+TightID+iso is flat in η and charge-symmetric at 0.95, so step 8 is the sole source. It is what makes the μ⁺ (A·ε)(η) dip to 0.81 in the barrel while μ⁻ stays at 0.92 (endcap 0.96 both) in `xsec_fiducial_diag`. **Electrons are much less affected** (step 8: e⁺ −1.56%, e⁻ −0.79%) because L1 EG objects are ECAL-position-based, like the offline electron φ. Applied to data and MC alike, so it largely cancels in r — **but it rides directly on the charge-asymmetry observable and discards ~5.7% of μ⁺ events in a stats-limited measurement.** **The applied fix: `trigMatchDR` = 0.4** (chosen 2026-08-12; max |Δφ| is ~0.10 at threshold, so 0.4 clears it with margin). Measured in MC on the barrel pT 25–30 slice, matching efficiency μ⁺/μ⁻ goes 0.479/0.906 → **1.000/1.000** (ΔR<0.3 already saturates; |Δη|<0.1-only gives 0.996/0.996). NB correcting only the half-cell offset would fix the *charge* bias (diff −0.427 → −0.045) but still lose ~28% of BOTH charges at threshold — the window had to widen regardless. NB the +0.006 rad offset `c` is itself worth understanding — if it is a geometry/propagation artifact it need not be identical in data and MC, which would turn this into a direct charge-asymmetry bias rather than a cancelling one.

**Output**: per-sample ROOT files in `skim/rootfile/` (W: `{channel}_pO_PFMet_{sample}_hist.root`;
Z: `ZToMuMu_pO2025_*` / `ZToEE_pO2025_*`). W files contain 12 rapidity-binned MT
histograms (`h_mt_Wp_y0..y11`, `h_mt_Wm_y0..y11`), MET distributions, isolation
histos, **plus leading-lepton kinematics `h_lepPt/h_lepEta/h_lepPhi`**
(2026-07-29: leading lepton after the FULL selection incl. iso + trigger match,
charges combined; deliberately the same names as the Z-channel histos so
`correction/dataMC_kinematics.C` reads W and Z files uniformly). **Both `skim_Wmu` and `skim_Wel`** additionally store the ABCD inputs for
the QCD/low-MET background: per-charge 2D histos `h_iso_met_{mu,ele}{Plus,Minus}`
(relIso × PF MET) and `h_iso_mt_{mu,ele}{Plus,Minus}` (relIso × m_T), filled for
every event passing the full W selection **except** the isolation cut, so the
relIso axis spans iso-pass (signal) and the anti-iso sideband. Binning: relIso
0–1.0 in **0.005** steps (`kNIsoAB = 200`); MET 0–120/2 GeV; m_T 0–200/2.5 GeV.
Regions are chosen downstream by projecting these (no re-skim) in
`correction/qcd_abcd.C`. **The 0.005 width is load-bearing** (fixed 2026-07-30):
every boundary in use must be a bin EDGE — the electron cut 0.095 (= 19×0.005),
the muon 0.15, the anti-iso edges 0.20/0.30/0.60/0.65. With the previous 0.01
bins, 0.095 fell mid-bin and `qcd_abcd.C`'s `FindBin(isoCut−eps)` silently
integrated relIso < 0.10, leaking 170 events that fail `skim_Wel`'s own iso cut
into the ABCD iso-pass region and inflating T + every electron QCD template by
~4% (muon was always exact). Any new cut value applied by projection must land
on an edge.
**2026-07-29 additions (both W channels):** a third ABCD plane
`h_iso_pt_{mu,ele}{Plus,Minus}` (relIso × leading-lepton pT, 0–100/2 GeV) —
supplies the anti-iso lepton-pT *shape* for the QCD pT template (NOT used for
the 2×2 counting; relIso is pT-correlated) — and the rapidity-binned
per-charge lepton-pT histos `h_leppt_W{p,m}_y{0..11}` (+`_FB`), 2 GeV bins,
the pT twins of `h_met_*`/`h_mt_*` feeding the display-only lepton-pT stacks
in `mtandmet.C`. **2026-07-30:** joint per-charge (pT × m_T) 2Ds for the
cut-pair scan (`correction/ptmt_scan.C`): `h_pt_mt_*` (iso-pass) and
`h_pt_mt_antiiso_*` (the anti-iso sideband, relIso in the qcd_abcd windows,
with its actual MET/m_T — the scan's QCD model); pT 1 GeV bins, m_T 2.5 GeV.
**2026-08-04: these four scan planes are filled down to leading-lepton
pT > 20** (`kPtScanFloor` in `skim_common.h`, = the scan grid's lower edge) so
the scan can probe LOWERING the nominal pT cut. They are the ONLY histograms
with the relaxed floor — every other histogram AND the printed cutflow keep
the nominal pT > 25 selection (`passPtNominal` gates the fills; `hasPF25`
gates N[1..5] so the cutflow counts are unchanged; verified bin-identical on
the nominal set after the re-skim).
**The `_mt40` set (2026-07-30) = the LEPTON-pT-DISCRIMINANT selection**
`pT > 25 && m_T > 40` (local `mtCutForPtDisc` in both W functions; the m_T cut does
the QCD suppression that the MET *shape* does in the MET fit). Twins of the
existing histos with that cut applied: `h_leppt_mt40_W{p,m}_y{0..11}`(+`_FB`),
`h_lepPt_mt40`/`h_lepEta_mt40`/`h_lepPhi_mt40`, and the ABCD plane
`h_iso_pt_mt40_{mu,ele}{Plus,Minus}`. The nominal (no-m_T-cut) histograms, the
MET/m_T histos and the cutflow are untouched — the W skim still applies no
MET/m_T cut. Z files contain dilepton mass (`hMass`, `hMass_extended`, `hMass_vipul`)
plus kinematics of the dilepton system and its leptons (`h_Zpt/h_Zeta/h_Zy/h_Zphi`,
`h_lepPt/h_lepEta/h_lepPhi` — both legs; `h_Zy` is lab-frame rapidity), filled for
iso-selected OS pairs in the Z peak [60,120] GeV (consumed by
`correction/dataMC_kinematics.C`). **`skim_Zmm` additionally** stores the hadronic
recoil for the MET correction (Sec 6 of AN2017_058): `h_uPar`/`h_uPerp` (1D
inclusive) and `h_uPar_qT`/`h_uPerp_qT` (2D, recoil component vs q_T), where
`u = −MET − q_T` and `q_T` is the dimuon pT. u-axis = AN's 2 GeV/c binning; q_T
axis fine so the fit binning can be chosen later by projecting (no re-skim). PF
tree is **soft-gated** (loud `[WARN]`, empty recoil histos if `pftree` absent).
**`skim_Zee` stores the same four recoil histos** (`u = −MET − q_T`, `q_T` = the
dielectron pT), wired identically (soft-gated PF tree); consumed by
`correction/recoil_raw_ele.C`.

Per-job logs in `skim/logs/`, cutflow text in `skim/output/`.

### Stage 2 — `plotting/`

Data/MC overlay, background estimation, plots.

- [plotting/mtandmet.C](plotting/mtandmet.C) — MT and MET data/MC plots. **MET stacks (muon AND electron) include the data-driven ABCD QCD** from `correction/rootfile/qcd_abcd_{mu,ele}.root` (run `correction/qcd_abcd.C[+(true)]` first): added to the inclusive MET stack (rigorous) and split across the per-y bins by the per-bin low-MET excess (data−EWK at MET<30) proxy, keeping the inclusive template shape (per-y sums back to inclusive). **Also writes the structured Combine input `plots[/Elec]/combine_input_W.root`** — one TDirectory per fit region (`Wp_lab_y3/`, `Wm_fb_y7/`, `Wp_incl/`, `W_incl/`, …) each with the 6 absolute templates `data_obs/signal/z/ztau/wtau/qcd` (consumed by the fork's `run_pO_fits.sh` — see "Downstream fit"). Channel auto-selected by `isElec`. **Lepton-pT stacks (2026-07-29, DISPLAY ONLY):** the same absolute data-vs-stack comparisons in leading-lepton pT, in **two variants** driven by one 2-entry `varStem[]`/`varTag[]` table (`kVarNom`/`kVarMt40`): `plots[/Elec]/leppt/` = plain W selection (`h_leppt_*` + `qcd_pt_*`), and **`plots[/Elec]/leppt_mt40/` = the pT-DISCRIMINANT selection pT>25 && m_T>40** (`h_leppt_mt40_*` + `qcd_pt_mt40_*`, 2026-07-30). Both: per-y `leppt[_mt40]_W{p,m}_y*` + `_FB` + inclusive, per-y QCD split by the same low-MET-excess weights (the m_T cut doesn't move QCD in y). **`combine_input_W.root` itself stays PF-MET-only**, but since 2026-07-30 the same writer also emits the two lepton-pT variants as **separate files** — `combine_input_W_leppt.root` (plain selection) and `combine_input_W_leppt_mt40.root` (pT>25 && m_T>40) — same 51 regions × 6 templates each, pT binning 50/0–100, consumed by the fork's `run_pO_fits.sh --disc leppt|leppt_mt40` (out-trees `pO_fit_out_leppt[_mt40]/`; `sync_lxplus.sh` uploads the variants when present and downloads all out-trees). qcd_norm stays FREE for all three (2026-07-30 decision). **Composition after m_T>40** (data / W sig / DY+τ / QCD): μ 4382 / 3704 (84.5%) / 272 (6.2%) / 226 (5.2%), closure 96%; e 4465 / 3080 (69.0%) / 151 (3.4%) / 1202 (26.9%), closure 99.3%. (Before the cut: W purity 73.5% μ / 45.0% e.) (For electrons QCD is the *dominant* low-MET component — see `qcd_abcd.C`.)
- [plotting/mtandmet_overlay.C](plotting/mtandmet_overlay.C) — MC stack overlays
- [plotting/dileptonpeak.C](plotting/dileptonpeak.C) — Z peak plots
- (QCD sideband fit and isolation ROC curves moved to `correction/` — see below)
- [plotting/observables.C](plotting/observables.C) — final charge-asym + F/B plots from the FITTED yields (Combine), with **all four** nPDF theory bands (EPPS21/nCTEQ15HQ/nNNPDF3.0/TUJU21nlo). **Simfit-aware since 2026-08-04:** `observables_comb(disc)` = the PRIMARY μ+e-combined plots of the grand fit (reads the disc-tagged `comb` graph files + `simfit/summary/comb_fitted_yields.root`, label "W → l ν / μ + e combined (simfit)") → `plots/comb/{charge_asym,FBratio}/<disc>/`; internals shared with the per-flavour path via `observables_run(chan, lepSym, outBase, …)`. **Discriminant-aware since 2026-08-03:** `observables(isElec, disc)` = per-channel (legacy fits), `observables_overlay(disc)` = merged μ+e overlay (legacy), with `disc` = `met`(default)|`leppt`|`leppt_mt40` driving BOTH the default inputs (the matching `pO_fit_out<suffix>/` fitted yields + the disc-tagged `{charge_asym,FBratio}_fit_<chan>_<disc>.root`) and the per-disc output folders `plots[/Elec]/{charge_asym,FBratio}/<disc>/`, `plots/merged/<disc>/`, so variants never overwrite each other; the tag is stamped on every plot as a 4th header line ("PF MET fit" etc.). Mapping single-sourced in [plotting/disc_variants.h](plotting/disc_variants.h) (`pODisc::Spec`, unknown tag ⇒ hard error; `pODisc::GraphFile` = tagged name with met-only legacy-untagged fallback). Reads `RpO_rootfile/RpO_FB_graphs.root` for theory (optional). Fitted yields read `h_yield_*` first (`h_mt_*` legacy fallback). Whole Module-5 chain in one command: [analysis/run_observables.sh](analysis/run_observables.sh). The old "Projection with Electrons" pseudo-band is gone (real e overlaid instead).
- [plotting/xsec_fiducial.C](plotting/xsec_fiducial.C) — fiducial W cross sections from the fit-summary CSVs (σ_fid = N_fit/L, stat-only, NO eff/acc correction — μ vs e gap = channel efficiency). `xsec_fiducial(disc)` = W⁺/W⁻/W incl, μ+e overlaid (reads `<chan>_summary.csv`); `xsec_fiducial_diff(isElec, disc)` = dσ_fid/dη_CM per channel (reads `<chan>_W_yields.csv`; η bins from `g_chargeAsym` in the disc-tagged charge-asym file). Same `disc` scheme as `observables.C` → `plots/xsec/<disc>/`; run by `run_observables.sh`. **`xsec_fiducial_comb(disc)` — THE r×σ_gen MEASUREMENT since 2026-08-12** (2026-08-05..12 it was the N/L + prefit/gen-overlay view): **σ_meas,i = r_i × σ_gen-fid,i, PER LEPTON FLAVOUR** (μ/e-shared r's ⇒ σ(W→ℓν) under lepton universality; the old "×2 = μ+e summed" display convention is gone). Reads per-bin `r/rErr` from `comb_W_yields.csv` and σ_gen-fid from `gen_xsec.root`; inclusive W⁺/W⁻/W = Σᵢr_iσ_i with the full r-covariance recovered as cov_Y_ij/(S_i·S_j) from `h_cov_yield` × the `signal_prefit` column (diagonal-rErr fallback + WARN; gen MC-stat ~0.3%/bin neglected vs ~6% fit errors). Plots show measured points vs the gen-fiducial r=1 expectation (green diamonds / dashed lines — measured÷dashed per bin IS r_i; the prefit-open-marker overlay was dropped, that role is absorbed). Also writes the σ table `plots/comb/xsec/<disc>/xsec_comb.csv` (per-bin + inclusive rows; incl `r` col = gen-weighted mean r). **The COUNT record is untouched** — `comb_W_yields.csv`/`comb_summary.csv` keep the fitted yields (echoed as a console cross-check), and charge_asym/FBratio stay count-based by design. First numbers, per flavour (post η-convention fix): met σ(W⁺)=50.36±0.91 / σ(W⁻)=39.68±0.80 / W=90.03±1.21 nb (eff. r≈1.01); leppt_mt40 54.20±1.37 / 43.89±1.08 / 98.09±1.75 (eff. r 1.09/1.11) — **the ~9% met↔leppt_mt40 spread = the e pT-tail mismodeling absorbed differently, now THE visible inter-discriminant systematic**. → `plots/comb/xsec/<disc>/{W_fiducial,W_dsigma_deta_comb}`. **`xsec_fiducial_diag(disc)` — PER-FLAVOUR EXTRACTION DIAGNOSTIC (2026-08-12)**: two plots (μ, e), each with 4 series × {W⁺,W⁻,W} — σ_gen-fid, σ_reco=S_f/L, r×σ_gen (the measurement, identical in both since r is shared), r×σ_reco=N_fit,f/L (count-based) — so each gen/reco pair differs by exactly (A·ε)_MC,f, printed on the plot. Per-flavour S_f comes from THIS repo's `plots[/Elec]/combine_input_W<suffix>.root` (`W{p,m}_lab_y*/signal` integrals; μ+e reproduce the CSV `signal_prefit` exactly), r-errors via the same `sumWithCov`; Σ_flavour of r×σ_reco reproduces the fork's `*_sum_yield/L` exactly (printed cross-check) → `W_xsec_diag_{mu,ele}` + `xsec_diag.csv`. **(A·ε)_MC: μ 0.90/0.94/0.92 vs e 0.76/0.77/0.76 (W⁺/W⁻/W, met)** — the μ/e gap is the electron ID+ECAL-crack loss, i.e. the visual case for the ECAL-gap TODO; leppt_mt40 shifts them to μ 0.89 / e 0.74 (the m_T cut's efficiency). **The η-binned diagnostic (`W_dsigma_diag_<flav>_<charge>` + the `AxEps_vs_eta` summary) resolves the crack directly**: both electron curves dip to ≈0.62–0.64 in exactly the two |η_lab| ∈ [1.2,1.6] bins (containing 1.4442–1.566), with no such feature in muon — and, being detector response, the W⁺ and W⁻ curves now lie on top of each other (that charge-symmetry is what exposed the η-convention bug above; **use it as the standing sanity check on this plot**). Muon A·ε is ~0.96 at both η edges and dips to 0.81(W⁺)/0.92(W⁻) centrally — the edge rise is pT-threshold in-migration (the lepton pT spectrum is softer and falling at 25 GeV at large |η|, so smearing adds more than it removes), which is why it is charge-dependent while the crack is not. The bin-by-bin correction is not an unfolding, so per-bin σ carries this migration; the inclusive σ and the count-based observables do not. Both `xsec_fiducial_comb` and `_diag` share `loadCombIngredients` + `sumWithCov` (single-source of the G/S/R/cov ingredients). Run by `run_observables.sh` (comb chain). The μ-vs-e comparison plots remain legacy-only (meaningless under shared r's).
- [plotting/postfit_incl.C](plotting/postfit_incl.C) — **rapidity-INCLUSIVE simfit postfit stacks (2026-08-05)**: `postfit_incl(disc)` sums the 12 LAB bins per flavour for W⁺/W⁻/W (6 plots) → `plots/comb/postfit_incl/<disc>/postfit_{mu,ele}_{Wp,Wm,W}`. Built WITHOUT fitDiagnostics: since every simfit parameter is a pure normalization (no shape nuisances), postfit shape ≡ prefit template × fitted scale, so the stacks are reconstructed exactly from `combine_input_W*.root` × the fitted params in `comb_W_yields.csv` (r per bin, r_Z, per-channel qcd_norm) — runs locally on downloaded summaries. **This shortcut breaks if shape systematics are ever added** (then use shapes_fit_s). Same cosmetics as the per-bin postfit plots (pull pad); info box: data/total, Baker-Cousins χ², r_Z. First results (leppt_mt40): μ W-incl χ²/ndf 1.46 (p 0.03); **e W-incl 3.43 (p≈0) — systematic data excess over the model across the pT tail 45–95 GeV**, i.e. the electron pT-shape mismodeling aggregated (links: e-scale/calibration shift, unembedded MC, anti-iso shape tilt; the per-bin fits partially absorb it via r/qcd_norm). Run by `run_observables.sh` (comb chain).
- [plotting/plotZcurve.C](plotting/plotZcurve.C), [plotting/plotRpOtheory.C](plotting/plotRpOtheory.C) — `plotRpOtheory.C` PRODUCES the theory graphs (`filelist_theory.txt` → `RpO_rootfile/RpO_FB_graphs.root`, channel/yield-independent); plotZcurve = Z curve plot
- [plotting/CMS_lumi.C](plotting/CMS_lumi.C), [plotting/plotting_helper.C](plotting/plotting_helper.C) — style / helpers. **Pull sub-pad (2026-07-19):** `SaveNicePlot1D_WithBkg` grew an opt-in bottom pull pad (`ps.pullPad`, default off): per-bin `(data − MC_total)/σ` with the Poisson-correct `σ² = MC_total + σ_MCstat²` (NOT the observed data error — that blows up empty-data bins against a small-error MC tail), drawn as zero-anchored bars (`"B"`), red dashed 0-line + dotted ±2σ guides, auto-symmetric y-range `max(3, 1.15·max|pull|)` (fix via `ps.pullYRange`). The canvas stays NEAR-SQUARE (height × `ps.pullCanvasScale`=1.125, so 800×900); the main pad compresses slightly — everything is pad-relative so the layout scales consistently. Enabled in `mtandmet.C` (all MT/MET stacks, both channels) and the fork's `draw_postfit_pO.C`; `dileptonpeak.C` left single-pad (flip `ps.pullPad` there to add it).

### Stage 3 — `analysis/`

Final observable extraction from the rapidity-binned histograms:

- [analysis/analysis_helpers.h](analysis/analysis_helpers.h) — shared header. `pOAnalysis::YieldInRange`, `AsymErr`, `RatioErr`, and the `kPORapidityShift = 0.3466` constant (single source of truth for the pO→CM boost). **2026-08-04:** `AsymErr`/`RatioErr` take an optional trailing `cov` argument (default 0 = the old independent-yield formulas exactly) for the simfit covariance cross terms.
- [analysis/charge_asym.C](analysis/charge_asym.C) — `A = (N+ − N−) / (N+ + N−)` vs rapidity (12 bins, y ∈ [−2.4, 2.4]). Reads `h_cov_yield` if present (simfit `comb_fitted_yields.root`) and includes cov(N+, N−) per bin; absent → legacy behavior.
- [analysis/FBratio.C](analysis/FBratio.C) — forward/backward ratio in |y_CM|. Lab-frame `yEdges` chosen symmetric around Δy so the CM-frame bins are symmetric around 0 (required for F/B pairing). Reads `h_cov_yield_FB` if present: Var(F)/Var(B) get the within-sum cov terms (which `Yield::operator+` can't know) and the ratio error gets cov(F, B); absent → legacy behavior.
- [analysis/run_observables.sh](analysis/run_observables.sh) — **Module-5 driver (2026-08-03; simfit-aware 2026-08-04)**: `./run_observables.sh [met|leppt|leppt_mt40|all]` runs the whole fitted-yields→observables chain for ONE W-discriminant variant, carrying the disc tag through every filename/folder. Two conditional blocks per variant: **PRIMARY simfit chain** (when `pO_fit_out<suffix>/simfit/summary/comb_fitted_yields.root` exists): charge_asym + FBratio on it (covariance-aware) → `../skim/rootfile/{charge_asym,FBratio}_fit_comb_<disc>.root` → `observables_comb(disc)` → `plots/comb/{charge_asym,FBratio}/<disc>/` → `xsec_fiducial_comb(disc)` + `xsec_fiducial_diag(disc)` → `plots/comb/xsec/<disc>/`; **legacy chain** (when both `<chan>_fitted_yields.root` exist): charge_asym + FBratio (both channels → `..._fit_<chan>_<disc>.root`), then `observables(±, disc)`, `observables_overlay(disc)`, `xsec_fiducial[_diff](…, disc)` → `plots[/Elec]/{charge_asym,FBratio}/<disc>/`, `plots/merged/<disc>/`, `plots/xsec/<disc>/` (the μ-vs-e xsec comparison stays legacy-only — it IS the efficiency diagnostic, meaningless under shared r's; the comb chain has its own `xsec_fiducial_comb`). Strict single-variant errors only when NEITHER block has inputs; `all` skips such variants. Post-checks outputs (root exits 0 even when a macro bails). bash-3.2-safe; `FORK_TEST` env overrides the fork location.

### `correction/` — corrections & studies

Orthogonal to the main skim→plotting→analysis flow: code that derives/validates
corrections and background estimates. **All macros are run from `correction/`**
and write outputs there (`correction/plots/`, `correction/rootfile/`). The two
that need `plotting_helper.C` `#include "../plotting/plotting_helper.C"`; the
isolation study is self-contained (writes/reads `correction/rootfile/`); the QCD
fit reads the *main* W skim output at `../skim/rootfile/`.

- [correction/qcd_abcd.C](correction/qcd_abcd.C) — **current** data-driven QCD/low-MET background via the ABCD method, **both muon and electron** (`qcd_abcd.C+` = muon, `qcd_abcd.C+(true)` = electron). Plane = relIso × PF MET (also stores relIso × m_T). Regions: iso-pass (`relIso < isoCut`) / iso-fail (anti-iso window `[isoFailLo, isoFailHi)`) × MET high/low (`metCut`); all four counts are QCD-only (`data − Σ EWK·k_s`, absolute MC from `mc_norm.h`). Reports the transfer factor `T = B/D`, the signal-region QCD `A = B·C/D`, and writes the iso-pass QCD MET/m_T **template** (anti-iso shape × T) to `correction/rootfile/qcd_abcd_{mu,ele}.root` for Combine, plus closure overlays (data vs EWK ± ABCD QCD). Channel diffs: electron uses `isoCut=0.095`, anti-iso `[0.20,1.0)` (muon `0.15`, `[0.30,1.0)`) and the `Wp_ele/Wm_ele/DYee` MCScale labels. Boundaries are projections → retune `ABCDConfig` with no re-skim.
**Lepton-pT template (2026-07-29):** also writes `qcd_pt_{mu,ele}{Plus,Minus}` =
(anti-iso pT shape from `h_iso_pt_*`) × (the **MET-plane** T) — relIso×pT is
correlated for QCD, so the pT plane is never counted with the 2×2; template
total = the MET-plane iso-pass QCD by construction. **`qcd_pt_mt40_*`
(2026-07-30)** is the same built from `h_iso_pt_mt40_*`, i.e. the QCD template
of the pT-discriminant selection (pT>25 && m_T>40), same MET-plane T (QCD
isolation efficiency assumed independent of the recoil variable) ⇒ its total
IS the ABCD prediction for QCD surviving the m_T cut. **m_T>40 keeps only
≈26%/25% of the μ⁺/μ⁻ QCD and ≈33%/34% of the e⁺/e⁻ QCD** (μ 413+454→111+116,
e 1881+1687→623+579). Ships an anti-iso
slice-stability shape check (`antiiso_shape_pt_*`): μ stable; e shows a mild
tilt (sideband pT slightly soft vs iso-pass) — the known systematic knob. **Result:** muon QCD is small (T≈0.17; ≈3% of the high-MET signal region). Electron QCD was originally huge with the loose `eleCutIdWP95` (T≈0.30; ≈24–26%), which **motivated the ID switch to `eleMVAIdWP95`** (see the isolation/ID study + FIXMEs); after the switch the electron QCD dropped to **T≈0.32–0.38, ≈10–11%** of the high-MET signal region (iso-pass QCD ≈10037→3634, W signal only −6.5%, S/B 0.34→0.87). Closure is good in both channels. (Electron T/QCD came down a further ~4% on 2026-07-30 with the relIso bin-edge fix — see the `kNIsoAB = 200` note above. **MET-plane** values, μ⁺/μ⁻ and e⁺/e⁻: μ T = 0.2001/0.1832 **unchanged** by construction, e T = 0.3967/0.3368 → **0.3802/0.3243**; iso-pass QCD μ 413.7+453.8, e 1969.5+1758.9 → **1887.2+1693.9**. NB the m_T-plane has its own T (μ 0.178/0.163, e 0.343/0.323) — don't quote those as the MET-plane numbers.)
- [correction/ptmt_scan.C](correction/ptmt_scan.C) — **(pT × m_T) cut-pair scan** (2026-07-30, iso NOT scanned — stays at nominal): which (lepton-pT, m_T) cut combination optimizes the W selection. Signal = absolute W MC; background = the **anti-iso sideband** (full selection except iso, relIso in the qcd_abcd windows [0.30,1.0) μ / [0.20,1.0) e, ACTUAL MET/m_T) EWK-subtracted and shape-normalized to the measured iso-pass QCD total, + non-W EWK counted directly. (An earlier MET<5 "low-MET proxy" convention was **removed same day**: conditioning on MET caps the proxy's m_T at ~2√(pT·MET)≈30 GeV — a scan of m_T must not use a MET-conditioned background. Residual anti-iso caveat: mild iso–pT shape correlation, bounded by the `antiiso_shape_pt_*` slice checks.) Reads the joint `h_pt_mt_{mu,ele}{Plus,Minus}` (iso-pass) / `h_pt_mt_antiiso_*` (sideband) 2Ds from the skim. Outputs S/√(S+B) & S/B maps (optimum ★, tentative (25,40) ✚), sweep ROCs + AUC, m_T profiles → `plots/ptmt_scan_{mu,ele}/` + `rootfile/ptmt_scan_{mu,ele}.root`. **Result (anti-iso, final):** pT>25 always best (raising it only cuts the Jacobian peak; pT-sweep AUC ~0.71 vs m_T-sweep ~0.94); μ plateau flat over m_T 40–45 (S/√(S+B) 56.90–56.96, ε_S 95–97%; optimum (25,42.5)); e peaks at (25,52.5): 37.8→48.5, S/B 0.82→4.7, ε_S 89.7% (at (25,45): 47.7, ε_S 95.2%). **Channel-uniform compromise (25, 45)** sits within 0.1% (μ) / 1.7% (e) of the per-channel optima. Sideband stats: 4579 μ / 10239 e events. `root -l -q 'ptmt_scan.C+'` / `+(true)`. NB an m_T precut suits the lepton-pT-fit path; it is NOT compatible with the current MET-shape fit (it would remove the QCD-constraining low-MET region). **2026-08-04: pT grid extended down to 20** (skim scan planes now filled for pT > `kPtScanFloor` = 20 — see Stage 1; `kPtRef = 25` keeps ε_S/baseline quoted vs the CURRENT selection, so ε_S > 100% = signal gained). Result: WITH the m_T cut the FOM keeps rising as pT drops — μ optimum moves to **(20, 42.5)**: S/√(S+B) 56.99→58.14 (+2.0%), S +5.4% (3676→3874), S/B 7.6→6.9, still rising at the 20 floor; e is ~flat: (20, 50) 49.86 vs (25, 50) 49.21 (+1.3%), i.e. same S as the old (25, 40) point at S/B 4.2 vs 2.7 — but B(mT>40) grows +55% and the e sideband tilt caveat grows with the 20–25 slice (sideband: 7700 μ / 25902 e events, of which ~15k e sit in pT 20–25). WITHOUT an m_T cut (= the current MET-shape-fit selection) lowering pT only hurts: μ 53.1→52.1, e 39.6→33.6 (e QCD B ×2.4). So pT→20 is a real (small) win for the μ lepton-pT-discriminant path, a wash for e, and NOT advisable for the MET-fit selection. NB pT≥25 FOM values shifted slightly vs the 2026-07-30 numbers (μ ≤0.2%, e ~1.5% up, e best-mT 52.5→50) because the sideband-shape normalization `qcdTot/sideTot` is now anchored over the extended pT>20 plane — the iso–pT-correlation systematic in action, not a code change (per-slice transfer r = QCD_isopass/sideband: μ 0.208 in 20–25 vs 0.219 at ≥25, −5%; e 0.245 vs 0.345, −29%; the mixed anchor rescales the ≥25 QCD by ×0.98 μ / ×0.82 e). ROC/AUC now normalized to the loosest point (20, 0) (NOT comparable to the old AUCs); per-pT-cut optimum table (pT 20–30) printed in the report.
- [correction/qcd_sideband_fit_and_extrapolate.C](correction/qcd_sideband_fit_and_extrapolate.C) — **superseded** by `qcd_abcd.C`. QCD shape from anti-iso sideband, Rayleigh-like fit + linear shape-parameter extrapolation to signal iso; did not behave at pO statistics. Kept for reference.
- [correction/isolation_mu_tight.C](correction/isolation_mu_tight.C) (muon, **current**), [correction/isolation_ele.C](correction/isolation_ele.C) (electron) — isolation ROC study, both on the **Δβ-corrected** relIso since 2026-07-06. Signal = OS Z-window pairs; backgrounds = SS pairs (≈empty for μ) and single-lepton MET<5 (QCD-like). `isolation_mu_tight.C` is the fresh, pure muon study (2026-07-06): ONE ID (`muIDTight`, pT>15, |η|<2.4), ONE variable (branch-based Δβ relIso = the skim's `RelIsoPF`), continuous 200-point scan, dbeta + nodbeta(reference) tags → `rootfile/IsoStudyOutputs_muon_tight.root`; run `root -l -q 'isolation_mu_tight.C+("/Users/zhenghuang/pO_2026_May_26/Data_May_26.root")'`. **Electron `isolation_ele.C`** scans the 8 ID WPs × 4 MVA-Iso WPs AND a continuous-relIso scan per ID (pass the current data file — default path is stale EOS). The older multi-cone/multi-def muon scan [correction/isolation.C](correction/isolation.C) is kept untouched for reference (uncorrected iso). **Conclusions (re-derived under Δβ, 2026-07-06 on May-26 data):** μ TightID AUC_MET=0.943, J(QCD) optimum 0.156 ⇒ relIso<0.15 stands (ε_sig=0.934, ε_QCD=0.130; Δβ vs uncorrected nearly identical — pO pileup tiny, PU-iso nonzero for only ~20% of leptons); e `MVAIdWP95` AUC_MET=0.909, J optimum exactly 0.095 ⇒ cut stands (`CutWP95` remains the worst ID, AUC 0.868).
- [correction/PlotsIsoROC.C](correction/PlotsIsoROC.C) (`PlotsIsoROC.C+(false)` for the tight set), [correction/PlotIsoROC_ele.C](correction/PlotIsoROC_ele.C) — isolation ROC curves → `correction/plotsROC[_ele]/`; plus `plotsROC_ele/roc_MET_continuous_allID.png` (one-off overlay of the 8 IDs' continuous-relIso QCD-ROC).
- [correction/plot_iso_summary.C](correction/plot_iso_summary.C) — **per-fixed-ID summary in one consistent cosmetic**: a 2-pad canvas (LEFT = continuous efficiency vs Δβ-relIso cut for signal/QCD/SS with the operating cut marked + ε_sig/ε_QCD; RIGHT = the corresponding ROC + AUC with the operating-point star). Electron → one per ID (`plotsROC_ele/summary_<ID>.png`, reads the continuous scan from `isolation_ele.C`); muon → single tight-ID summary (`plotsROC/summary_muon.png`, reads `IsoStudyOutputs_muon_tight.root` from `isolation_mu_tight.C` — since 2026-07-06; before that the old `ggbranch` scan). `root -l -q 'plot_iso_summary.C+'`.
- [correction/dataMC_kinematics.C](correction/dataMC_kinematics.C) — Data vs signal-MC overlay + ratio pad, shape-normalized. **Z channels** (`"Zmm"`/`"Zee"`): DY MC vs data for the Z kinematics (`h_Zpt/h_Zeta/h_Zy/h_Zphi`, `h_lepPt/h_lepEta/h_lepPhi`, `hMass`) — the Data/MC boson-pT check. **W channels** (`"Wmu"`/`"Wel"`, added 2026-07-29): leading-lepton `h_lepPt/h_lepEta/h_lepPhi` after the full W selection, data vs W⁺+W⁻ MC combined with the relative `k_s` from `mc_norm.h` (shape-norm cancels the absolute scale; backgrounds NOT subtracted — and since the W selection has NO MET cut, QCD here is ≈19% μ / ≈50% e of the plotted sample (data−EWK, 2026-07-29; the familiar 3%/10% are high-MET-region figures), so expect a low-pT data excess — it's a shape check). Probes whether lepton pT is well-enough modeled to serve as an alternative fit discriminant (Jacobian peak, MET-free). **`"Wmu_mt40"`/`"Wel_mt40"` (2026-07-30)** run the same check on the pT-discriminant selection (pT>25 && m_T>40) via the `_mt40` histos → `plots/dataMC_W{mu,el}_mt40/`; with QCD down to ~5% (μ) / ~28% (e) there, the shape comparison is far less background-dominated. `root -l -q 'dataMC_kinematics.C+("Wmu")'` etc.; W data file has no sample suffix (`<base>_hist.root`).
- [correction/recoil_raw.C](correction/recoil_raw.C) — raw look at the Z→μμ hadronic recoil (`u_par`/`u_perp`) for the MET recoil correction: data vs DY MC, inclusive and in q_T slices (projected from the 2D histos), **plus a printed entry-count table per q_T slice** to judge statistics/binning before fitting. Slice edges = editable `kQtEdges` array (projections → no re-skim to change). **Checked 2026-07-02: recoil looks stable for now** — no correction applied, the planned `recoil_fit.C` (double-Gaussian per q_T bin → μ(q_T)/σ(q_T), AN Sec 6) is deferred; the MET data/MC discrepancy is attributed to the unembedded MC samples (no pO underlying event), not recoil. `root -l -q 'recoil_raw.C+'`.
- [correction/recoil_raw_ele.C](correction/recoil_raw_ele.C) — electron-channel twin of `recoil_raw.C`: same raw-recoil look for Z→ee (reads `ZToEE_pO2025_*`, writes `plots/recoil_Zee/`). `root -l -q 'recoil_raw_ele.C+'`.
- [correction/njet_WZ.C](correction/njet_WZ.C) — **jet multiplicity in W-/Z-tagged events** (2026-08-05). Loops the NTUPLES directly (skim outputs carry no jet info), replicating the skim selections exactly: W = the full 8-step selection (verified event-identical — data counts equal the skim cutflow N[8]: 5203 μ / 7065 e), Z = first iso-selected OS pair in the peak [60,120] (356 μμ / 248 ee; event counted once). Jets = `ak4PFJetAnalyzer/t` (entry-aligned with EventTree, equal-entry-counts enforced; verified run/evt-identical on July-29), calibrated `jtpt` > 30 GeV, |`jteta`| < 2.1, cleaned by ΔR > 0.4 against the selected lepton(s) (both Z legs — essential for electrons, which always double as a PF jet); no jet ID (PF-composition branches exist if wanted). Data vs signal MC (W: Wp+Wm; Z: DY), gen-weighted, shape-normalized `SaveDataMCRatio` overlays; W channels also fill the `_mt40` njet twin (pT>25 && m_T>40 — data 4382 μ / 4465 e, matching the known composition counts). Outputs `rootfile/njet_<chan>.root` (per-sample raw + `*_mc` combined at the ABSOLUTE `k_s` scale — the combined totals close on the known expectations: Zμμ 372.1 vs data 356, Zee 252.5 vs 248, W μ 3821.7 / e 3178.7) and `plots/njet_<chan>/{njet[,njet_mt40],jetpt,jeteta}` (header lower-right, `ps.yHeadroom` opt-in added to `plotting_helper.C` for the log-pad top room — legacy 1.4 default untouched elsewhere). Second arg = plots-only mode: `njet_WZ.C+("all", true)` redraws every plot from the stored rootfiles in seconds (cosmetics iterations, no ntuple loops). **Result (2026-08-05):** Z njet is well modeled (⟨njet⟩ data/MC: μμ 0.183/0.160, ee 0.177/0.167, frac(≥1 jet) ≈ 14%); W data sits far above signal MC (μ 0.291/0.138, e 0.372/0.146) — that's the jet-rich QCD in the no-MET-cut W sample, shrinking under m_T>40 (μ 0.222, e 0.300); the remaining excess = residual QCD (not subtracted) + real W+jets vs the unembedded MC. `root -l -q 'njet_WZ.C+'` (all four channels) or `njet_WZ.C+("Wmu")` etc.

The shared ratio-pad helper `SaveDataMCRatio(...)` lives in [plotting/plotting_helper.C](plotting/plotting_helper.C).

### `merge_rootfile/` — input prep

Utilities to scan EOS and `hadd` raw ntuple files into the **per-sample** ROOT
files that `skim/run_all.sh` reads (now the May-26 production — see "Input data";
earlier it was a single consolidated `pO_2025.root`):
- [merge_rootfile/make_filelist.sh](merge_rootfile/make_filelist.sh) — discover files on EOS
- [merge_rootfile/hadd_from_list.sh](merge_rootfile/hadd_from_list.sh) — merge them
- `MC_*.txt`, `DATA_pass_*.txt`, `version_*.txt` — per-sample filelists

## Input data

**Per-sample** ROOT files (the **July-29 production** since 2026-07-30; May-26
before that), not a single consolidated file. Single source of truth:
`kDefaultDataFile` (data) and `ResolveMCSample` (MC) in `skim/skim_common.h`
(`run_all.sh` greps the data path; override per-run via the `DATA_FILE` env var).

Currently pointing at the **local** copy under `~/pO_2026_July_29/` — written as
absolute paths (`/Users/zhenghuang/pO_2026_July_29/…`) because ROOT's
`TFile::Open` doesn't reliably expand `~`. To read off EOS/lxplus again, repoint
`kDefaultDataFile` + `inputBase` (in `ResolveMCSample`) to
`root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2026_July_29/`.
Files: `July_29_DATA.root` + 9 MC (`July_29_MC_DY_{mu,ele,tau}_Z.root`,
`July_29_MC_W{p,m}_{mu,ele,tau}.root`). No low-mass DY files in this production.
The new ntuples carry extra branches vs May-26 — harmless: the fast skim
enables only the branches it uses. **N_gen labels stay canonical** (`Wp_mu`,
`DYee`, `Wp_tau`, …): `count_ngen.C::LabelFromFname` maps the new filenames to
the legacy labels (`July_29_MC_DY_ele_Z.root` → `DYee`), so every
`pONorm::MCScale(label)` call site is production-independent. Merged from EOS
filelists by `merge_rootfile/run_all_hadd_July_29.sh` (glob-autodetects
`July_29_*.txt`, skips already-merged outputs; `FORCE=1` to redo). The May-26
local copy remains at `~/pO_2026_May_26/`.

Trees per file: `ggHiNtuplizer/EventTree` (main physics tree),
`hiEvtAnalyzer/HiTree` (gen `weight`), `hltanalysis/HltTree`, and PF candidate /
muon / electron collections.

## Downstream fit (separate repo)

Fits run from a fork of HiggsAnalysis-CombinedLimit at
`/Users/zhenghuang/HiggsAnalysis-CombinedLimit`.

**The working branch is `zheng/po-analysis`** — that's where all fit code lives.
`main` in that repo is stock unmodified Combine v10.6.0 (a tracking branch, not
a working branch). If the checkout is on `main`, switch with
`git checkout zheng/po-analysis` before doing anything, or read files via
`git show origin/zheng/po-analysis:<path>`.

### Current pipeline — `test/run_pO_fits.sh` (grand simultaneous fit since 2026-08-04)

Run under `cmsenv`:
```bash
./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|simfit|all] [--dry-run] [--no-postfit] [--draw-only] [--asimov]
```

**`simfit` is the DEFAULT mode (2026-08-04) — the GRAND SIMULTANEOUS FIT.**
One likelihood per binning variant (lab, fb — same events rebinned, so fitted
separately): all 48 W channels ({mu,ele} × {Wp,Wm} × y0..11, channel names
`<F>_<C>_<B>_y<i>`) + BOTH Z peaks (`mu_Z_incl`, `ele_Z_incl`). 2N+1 = 25 POIs:
**`r_<C>_y<i>`** (24) scales that bin's W-related MC (`signal`+`wtau`) in BOTH
flavours' channels (μ/e SHARED — lepton universality; relative μ/e acc×eff from
MC, lepton SFs still not applied), and **`r_Z`** (the "+1") is ONE global scale
on all DY-related MC (`z`/`ztau` in every W channel + `zsig`/`ztau` under both
peaks — DY rapidity dependence trusted from MC, normalization pinned by the
peaks). `qcd_norm_<channel>` free per W channel (48). `w`/`wtau` under the Z
peaks FROZEN at absolute MC (measured 0.03–0.06 events vs 372/252 signal —
decision 2026-08-04). Implemented WITHOUT combineCards.py:
`my_script/make_pO_simfit_cards.sh` writes one 50-channel card + the
multiSignalModel map file per variant (`t2w_maps_simfit_<B>.txt` = THE model
definition; regex maps like `(mu|ele)_Wp_lab_y3/(signal|wtau)$:r_Wp_y3[1,0,10]`
— the trailing `/...` keeps y1 from matching y10, the `$` keeps `z` from
matching `ztau`); `text2workspace.py -P ...:multiSignalModel` +
`combine -M FitDiagnostics --skipBOnlyFit` (no `--rMin/--rMax` — there is no
POI named `r`). `--asimov` adds a prefit-Asimov (`-t -1`) closure fit per
variant: every POI must return 1, checked/PASS-FAIL'd by the extraction.
`my_script/extract_pO_simfit.C` reads the one `fit_s` per variant → all POIs +
the r-correlation matrix → `simfit/summary/comb_W_yields.csv` (first 9 columns
= legacy layout), `comb_summary.csv` (POIs, status/covQual, covariance-
propagated Wp/Wm/W inclusive sums, Asimov rows), and `comb_fitted_yields.root`:
`h_yield_W{p,m}_y{0..11}(_FB)` with yields = r×(S_mu+S_ele) (+`h_mt_*` aliases)
PLUS the 24×24 yield-covariance TH2Ds **`h_cov_yield` (lab) / `h_cov_yield_FB`
(fb)**, fixed order [Wp_y0..11, Wm_y0..11] (axis labels set). Postfit plots per
grand-fit channel (draw_postfit_pO.C grew optional trailing args poi/dy/qcd
name + ndfParams; defaults = legacy behavior). Statistical motivation: the
legacy per-bin scheme re-fits the SAME Z data 48× per flavour ignoring the
induced correlations, and never combines μ/e; the grand fit fixes both and
supplies the covariance the downstream error propagation now uses. μ/e-shared
r's ⇒ per-flavour observables are 100% correlated — the comb result is THE
result; per-flavour results are independent only from the legacy fits.

**The legacy per-flavour pipeline below (`perbin|incl|combined`) is UNCHANGED
and kept runnable for comparison** (user: "keep the original version"). Mode
`all` = legacy + simfit (channel `both`). The `simfit` mode ignores the channel
argument (always μ+e); its out-tree is `pO_fit_out<suffix>/simfit/`
{datacards, fits/simfit_{lab,fb}, postfit, summary}.
(`--draw-only` redraws postfit plots from an existing run — needs only `root` +
the `fits/` tree, e.g. after cosmetic `draw_postfit_pO.C` changes. **`--disc
met|leppt|leppt_mt40`** (2026-07-30, default met) selects the W discriminant:
reads `combine_input_W[_leppt[_mt40]].root`, writes to
`pO_fit_out[_leppt[_mt40]]/`, sets the postfit x-title; datacards/fit model
identical, qcd_norm free for all three. **Downstream continuation of the same
tag (2026-08-03):** `analysis/run_observables.sh <disc>` runs the whole
Module-5 observables chain per variant into disc-tagged files/folders — see
Stage 3 and the `observables.C`/`xsec_fiducial.C` bullets.)
It (1) locates this repo's **structured** Combine inputs, (2) generates all
datacards, (3) runs `text2workspace`+`combine -M FitDiagnostics` per fit region
into a clean output tree, (4) extracts fitted signal yields, (5) draws postfit
plots in the same cosmetics as `plotting/mtandmet.C`.

**Full step-by-step runbook (both repos):** this repo's `README.md` is the
single end-to-end procedure (skim → ngen → ABCD QCD → structured inputs → fit →
charge_asym/FBratio → observables). The fork's `test/README_pO_fits.md` is now
just the fit-stage quick reference pointing back to it.

Scripts on `zheng/po-analysis`:
- `test/run_pO_fits.sh` — master driver (bash-3.2 safe; `--dry-run` builds
  datacards without `cmsenv`). Legacy output tree: `test/pO_fit_out/<chan>/
  {datacards,fits/<region>/,postfit,summary}`; simfit tree:
  `test/pO_fit_out<suffix>/simfit/`.
- `test/my_script/make_pO_simfit_cards.sh` — simfit 50-channel datacards +
  multiSignalModel maps (the 25-POI model definition), lab + fb.
- `test/my_script/extract_pO_simfit.C` — simfit extraction: POIs + covariance →
  `comb_*` CSVs + `comb_fitted_yields.root` (+`h_cov_yield[_FB]`), Asimov
  closure check.
- `test/my_script/make_pO_datacards.sh` — generates **all** datacards: per-region
  W (`Wp/Wm` × `lab/fb` × `y0..11` = 48, each a **two-channel card fitted
  simultaneously with `Z_incl`**; falls back to standalone W-only if
  `combine_input_Z.root` is missing), per-charge `Wp_incl`/`Wm_incl`, `W_incl`,
  `Z_incl` (standalone), and the simultaneous `WZ` card.
- `test/my_script/extract_pO_yields.C` — reads each region's `fit_s`, computes the
  fitted signal yield (`r`×prefit-signal-integral, error `rErr`×same) → `summary/
  <chan>_W_yields.csv`, `<chan>_summary.csv`, and `<chan>_fitted_yields.root`
  containing single-bin `h_mt_W{p,m}_y{0..11}` (lab) and `..._FB` (FB) histos —
  the **exact names + Sumw2 errors** that `analysis/charge_asym.C` and
  `analysis/FBratio.C` read (just pass the file as their first arg → fitted yields
  replace raw yields, no macro edits).
- `test/my_script/draw_postfit_pO.C` — generic postfit data/MC plotter routed
  through the synced `plotting_helper.C::SaveNicePlot1D_WithBkg` (α=0.65 fills,
  width-3 outlines, `CMS_lumi`) → identical look to `mtandmet.C`.
- `test/my_script/plotting_helper.C` — **synced** from this repo's
  `plotting/plotting_helper.C` (was an older solid-fill copy; `CMS_lumi.C` was
  already identical). Re-synced 2026-07-19 with the pull-pad helper;
  `draw_postfit_pO.C` sets `ps.pullPad=true`, so postfit plots gain the pull
  sub-pad on the next fit run (the local `pO_fit_out/` tree has no `fits/`
  dirs left, so `--draw-only` can't regenerate them without re-fitting).

**Structured inputs** (produced HERE, consumed by the fork):
- `plotting/mtandmet.C` → `plots[/Elec]/combine_input_W.root`: one **TDirectory
  per fit region** (`Wp_lab_y3/`, `Wm_fb_y7/`, `Wp_incl/`, `W_incl/`, …), each
  holding the 6 **absolute** templates `data_obs/signal/z/ztau/wtau/qcd` (MET
  discriminant; per-y ABCD `qcd` split from the inclusive template). `signal` =
  both W MC samples in that reco-charge region; `wtau` = Wptau+Wmtau.
- `plotting/dileptonpeak.C` → `plots[/Elec]/combine_input_Z.root`: a `Z_incl/`
  dir with `data_obs/signal/w/wtau/ztau` (mass peak, absolute).

**Legacy fit model** (TWO-PARAMETER scheme, 2026-07-01; superseded by the
simfit 25-POI model above but still what `perbin|incl|combined` run;
systematics deferred):
- **Two MC scales per fit**: the POI **`r` = all W-related MC** (W `signal` +
  `wtau`; in simultaneous cards also the `w`/`wtau` backgrounds under the Z
  peak) and **`dy_norm` = all DY-related MC** (`z` + `ztau`; in simultaneous
  cards also the Z signal `zsig`). The relative composition WITHIN each group
  is LOCKED by the absolute `k_s` templates. Implemented by reusing the name
  `r` as a `rateParam` on the W backgrounds (signal is index 0 → `r` scales it
  automatically). Fitted W yield = `r`×signal-prefit (NO extra factor — the old
  `r`×`eff_lumi` is gone). Discriminant = **PF MET shape**.
  (Was: ONE shared `r` on all MC; before that `r` + separate `ewk_norm`.)
- Standalone `Z_incl` card: roles flip — the POI `r` IS the DY scale
  (`signal`+`ztau`), the W backgrounds `w`/`wtau` get a free `w_norm`.
- data-driven ABCD `qcd` → its own free `qcd_norm` rateParam.
- **Simultaneous cards — the `WZ` card AND every per-bin W card** (each per-bin
  card is a two-channel fit: the W region + `Z_incl`): `r` scales the W-related
  MC in both channels, the shared `dy_norm` scales the DY-related in both — the
  high-purity Z peak pins `dy_norm` (this REPLACES the old shared `eff_lumi`,
  which was degenerate with a separate DY scale). The fitted per-bin yields for
  `charge_asym`/`FBratio` come from these simultaneous per-bin fits. Sanity:
  `dy_norm ≈ 1` (`w_norm ≈ 1` for `Z_incl`). Postfit plots print `r`, `DY
  norm`/`W norm`/`QCD norm` in the info box (`draw_postfit_pO.C`).

All templates are **absolutely** normalized (`k_s = A·σ·L/N_gen`) — NO area
normalization (the old `make_combine_input*.C` `data_int/mc_sum` rescale is gone).

**Superseded** (kept for reference, not in the new pipeline): `test/run_fit.sh`,
`test/my_script/make_combine_input{,_Z}.C` (area-normalized + Rayleigh `pdfbkg`),
`testdatacard_{inclusive,Zmumu,Zee}.txt`, `draw_postfit_{inclusive,Zmumu,Zee}.C`.

Fit method: `combine -M FitDiagnostics --saveShapes --saveWithUncertainties`.

## Current state of the work

Active (in flight, recent commits):
- **GRAND SIMULTANEOUS FIT — simfit (2026-08-04, NEW DEFAULT)** — the fit
  procedure changed completely: one likelihood per binning variant (lab/fb)
  with all 48 (flavour, charge, y) W channels + both Z peaks; 25 POIs
  (`r_<C>_y<i>` μ/e-shared + global `r_Z`); implemented end-to-end (fork:
  `make_pO_simfit_cards.sh`/`extract_pO_simfit.C`/driver `simfit` mode +
  `--asimov` closure; repo: covariance-aware `charge_asym`/`FBratio`,
  `observables_comb`, run_observables simfit chain). **Validation pending:**
  needs the first lxplus run (`./run_pO_fits.sh --asimov`) → check Asimov
  closure PASS, fit status 0/covQual 3, r's vs the legacy per-bin values, and
  `r_Z ≈ 1`. Legacy per-bin pipeline intentionally kept runnable ("wrong way"
  reference: it re-fits the same Z data per bin). Details in "Downstream fit".
  DY vetoes included) and both isolation studies now use the Δβ PU-corrected
  relIso; cuts 0.15/0.095 re-confirmed as optima. **Requires a full re-skim
  (all channels, all samples) + downstream re-run** (ABCD QCD, plots, Combine
  inputs) — existing `skim/rootfile/` outputs still carry the old definition.
- **Electron isolation + background** — working point confirmed
  (`eleMVAIdWP95` + Δβ relIso < 0.095); PU correction for electron isolation
  DONE via the Δβ switch above.
- **Tau channels as background** — sample enums (`kWptau`, `kWmtau`, `kDYtau`)
  and `Wptau/Wmtau/DYtau` filelists exist; selection code largely there but
  validation pending. Commit `447a2b5` ("before adding tau") marks the boundary.
- **Corrections still WIP** — lepton scale factors and momentum
  scale/smearing are not applied yet. Don't assume MC is fully corrected
  when reading skim output. (The per-event *generator* weight IS now
  applied — see "Generator event weight" below — but these reco-level
  corrections are separate and still missing.) **Recoil: checked
  2026-07-02 via the raw-recoil study (`correction/recoil_raw[_ele].C`)
  and looks stable for now — no recoil correction applied; `recoil_fit.C`
  deferred. The MET data/MC shape discrepancy is understood to come from
  the current MC samples being UNEMBEDDED (no pO underlying event embedded
  in the simulation), not from recoil miscalibration.**
- **ABCD QCD background (muon + electron, done 2026-06-23)** — replaces the
  Rayleigh shape-extrapolation. Both `skim_Wmu`/`skim_Wel` store the relIso × MET
  (and × m_T) 2D per charge (`h_iso_{met,mt}_{mu,ele}{Plus,Minus}`);
  `correction/qcd_abcd.C[+(true)]` does `N_A = N_B·N_C/N_D` on QCD-only
  (EWK-subtracted) counts and emits the iso-pass QCD template, wired into the
  `mtandmet.C` MET stacks for both channels. Electron QCD is ~10× the muon and
  is the dominant low-MET background there. Remaining: use these templates in
  the Combine fork datacards.

Future plans / on the TODO list (not yet started):
- **σ extraction as r × σ_gen-fid — DONE (2026-08-12)** — the measured cross
  section is now quoted as fitted-signal-strength × gen fiducial σ instead of
  N_fit/(L·A·ε); implemented in `xsec_fiducial_comb` (see the Stage-2 bullet).
  Algebraically identical when A·ε comes from the same MC
  (r×σ_gen = N_fit/L/(A·ε)_MC), but kA_O and kSigma cancel between r and σ_gen
  (result independent of the assumed theory σ and the A=16 scaling), and future
  SFs on the reco templates propagate into r automatically. Counts CSVs and
  charge_asym/FBratio deliberately stay count-based. Residual model dependence
  concentrates at the pT=25 threshold (lepton scale — the known e-scale shift
  matters here) and in-fiducial efficiency.
- **ECAL transition gap veto (TODO added 2026-08-12)** — exclude
  1.4442 < |η_SC| < 1.566 in the electron RECO selection, both data and MC
  (skim_Wel + skim_Zee + the veto legs), so the measurement stops relying on MC
  modeling of the crack: no SF can fix an *acceptance* hole, only efficiency
  inside the fiducial can be corrected from data (T&P). Decision needed when
  implementing: whether the electron GEN fiducial also excludes the gap (μ/e
  fiducials then differ — awkward for the shared-r combined σ) or the common
  |η|<2.4 fiducial is kept and the small gap extrapolation is quoted as a
  systematic (leaning to the latter). Requires an electron re-skim + downstream
  re-run.
- **Tau cross-section measurement** — currently τ is only a background to W/Z.
  Adding W→τν / Z→ττ as their own signal channels is on the wishlist (tagged
  as "probably just an attempt").

Known FIXMEs (mostly in the electron path):
- `DrawDielectronPeak.C` — ECAL coverage comment says 2.4 but should be 2.5
  with the 1.4442–1.566 transition gap excluded
- **W→eν electron ID: DONE (2026-06-24)** — the isolation/ID study
  (`correction/isolation_ele.C`, now with a continuous-relIso scan) showed the old
  `eleCutIdWP95` was the *worst* QCD rejector; `skim_Wel` now uses **`eleMVAIdWP95`**
  (leading lepton + DY veto), continuous relIso < 0.095 kept (Youden-J optimum). This
  cut the electron ABCD QCD ~64% (≈10037→3634 iso-pass; signal −6.5%).
  **Electron gates HARMONIZED (2026-07-02):** `skim_Zee` switched to `eleMVAIdWP95`
  too (both legs; Z-window data 276→247, −10.5%, consistent with the tighter ID),
  and the `skim_Wel` DY veto's iso gate switched from the integer `eleMVAIsoWP95`
  WP to the **continuous `RelIsoPF < 0.095`** (the study characterized the
  continuous cut; the MVA-Iso WP was shown not to track it). Every electron ID
  gate in the skim is now `eleMVAIdWP95` and every iso gate is continuous relIso.
  Requires an electron re-skim (Wel + Zee, all samples) + downstream re-run.
- **Electron ΔβPU correction: DONE (2026-07-06)** — all relIso (both flavours,
  skim + studies) now Δβ-corrected; the 0.095 optimum was re-confirmed under the
  new definition (J_MET optimum exactly 0.095, AUC_MET 0.909). As anticipated,
  the correction changes very little at pO pileup.

Not in repo: no tests, no CI, no config files. Many recent commit messages are
placeholder ("xx") so `git log` is not a reliable narrative — read the diffs.

## When working on this repo

- Don't expect type-checking or a test suite. Sanity check by re-running the
  relevant `skim/run_all.sh <channel>` invocation and reading `skim/logs/`.
- After the skim refactor, common scaffolding lives in `skim_common.h` and
  the four channel selections live as functions in `skim/skim.C`. Most
  cross-channel changes (new cut, new histogram) can be made once in the
  header; only physics differences need to go in the per-function bodies.
  When in doubt, diff against `skim/legacy/` — that's the pre-refactor
  reference.
- Histogram naming is load-bearing — `analysis/charge_asym.C` and
  `analysis/FBratio.C` read by name, and so do the Combine input scripts.
  Rename in skim → must rename in analysis AND in the Combine fork.
  **Naming audit (2026-07-30):** names must state the quantity they hold.
  Fixed then: the ABCD region printout said "MET-high/low" even on the m_T
  plane (`ABCDConfig::metCut` → `yCut`, labels now derived from the plane);
  the fitted-yield containers were called `h_mt_*` although the fit
  discriminant is PF MET → the fork now writes **`h_yield_W{p,m}_y*(_FB)`**
  as primary with `h_mt_*` kept as a deprecated alias, and
  `charge_asym.C`/`FBratio.C` prefer `h_yield_*` and fall back (so old and new
  files both work in both directions). Same dual-name scheme for the output
  graphs: **`g_chargeAsym`, `g_RFB_{sum,Wp,Wm}`** primary + `_mt` alias;
  `observables.C` resolves via a `GetGraph(file, name, legacy)` helper and its
  plots are now `chargeAsym.png` / `RFB_{sum,Wp,Wm}.png`. `useMT` still
  genuinely selects m_T vs MET **for raw skim files only** — it is ignored when
  `h_yield_*` is present. Also fixed: a real crossed assignment in
  `correction/PlotsIsoROC.C` (the SS and MET efficiency graphs were drawn under
  each other's legend entry).
- Output paths in the skim functions are referenced again in plotting and in
  the Combine input scripts (`make_combine_input*.C`). Moving files mid-pipeline
  silently breaks downstream stages.
- `.so` / `_ACLiC_dict_rdict.pcm` files in `skim/` are ROOT's compiled-macro
  cache (now gitignored). Delete them if the macro signature changes and
  ROOT picks up the stale build.
- **Always create output directories before writing.** Output dirs (`rootfile/`,
  `plots/`, `output/`, `logs/`) are gitignored and absent on a fresh checkout
  (e.g. a pull on lxplus), so any code that writes must `gSystem->mkdir(dir,
  kTRUE)` (C macros) or `mkdir -p` (shell) first. This is already done in
  `skim/skim.C` (creates `rootfile/`), `skim/run_all.sh` (`rootfile output
  logs`), and every `correction/` macro (their `plots/` / `rootfile/`). When
  adding a new writer, do the same — and zombie-check input files you open.
- **Never run two `run_all.sh` invocations in parallel from `skim/`.** They share
  the ACLiC build artifacts (`skim_C.so`, `skim_C_ACLiC_dict.*`) in that one
  directory, so two simultaneous first-jobs race to rebuild them and one dies
  with `Error in <ACLiC>: Executing 'rootcling ...' failed!`. The driver reports
  it only as `failures=1` in its summary line, and because the compile (not the
  physics) failed, the previous output file is silently left in place — easy to
  mistake for success. Either run the channels sequentially, or pre-build once
  (`root -l -b -q -e '.L skim.C+'`) before launching parallel jobs. Hit
  2026-07-30 with `run_all.sh Wmu` + `run_all.sh Wel` started together.
- **Shell scripts must be bash-3.2 compatible.** The user runs locally on macOS,
  whose stock `/bin/bash` is **3.2** — no `declare -A` associative arrays (they
  crash with a cryptic `unbound variable` under `set -u`). `skim/run_all.sh` uses
  `case`-functions (`channel_enum`/`sample_enum`) instead. Keep new scripts free of
  bash-4-only features so they run on both the Mac and lxplus. (Homebrew bash 5 is
  at `/opt/homebrew/bin/bash`, but don't depend on it.)
- **Clone histograms read from a TFile before mutating them** (`Rebin`,
  `Scale`, normalization, `Add`, …). `file->Get("h")` returns the *file-owned*
  object; rebinning/scaling it in place mutates the shared histogram, so a
  repeated `Get` of the same name returns an already-mutated object (e.g.
  double-rebinned) and normalization comes out wrong. Always
  `h = (TH1D*)src->Clone(Form("%s_copyN", name))` with a unique name, then
  `h->SetDirectory(nullptr)`, and operate on the clone. This was the
  **`plotting/dileptonpeak.C` normalization bug** the user fixed in commit
  `0f15b45` (the `getRebinned` lambda rebinned the file-owned histogram in
  place). `plotting_helper.C::SaveDataMCRatio` and the Combine-input scripts
  already clone; audit any other macro that reads-then-mutates.
- **Fast skim: every tree disables all branches, then enables only the ones
  used.** All TTrees in `skim/skim.C` (all four channels) call
  `t->SetBranchStatus("*", 0)` and then `SetBranchStatus(name, 1)` for each
  branch read — `EventTree` is huge, so reading only ~5–15 branches instead of
  everything is a big speedup. **Pairing is load-bearing:** after `"*",0`, a
  branch that you `SetBranchAddress` but forget to `SetBranchStatus(name,1)` is
  **silently not read** — the variable keeps a stale value and the physics is
  wrong with no error. So always set status AND address together for every
  branch. For *optional* trees (e.g. `HiTree` in the Z channels, gated on
  `haveHiTree`) guard the `"*",0` with the have-flag, since the pointer may be
  null. (Brought the Z channels in line with the W channels 2026-06; W was
  already done.)

## Sumw2 / error tracking (audited May 2026)

ROOT histograms need `Sumw2()` to track per-bin sum-of-weights-squared,
otherwise per-bin errors are `sqrt(N)` (raw Poisson) and become wrong as
soon as the histogram is filled with weights or `Scale()`d.

State of the pipeline:

- **`skim/skim.C`** — all 84 TH1Ds call `Sumw2()` explicitly right after
  construction. Additionally, the `skim(...)` dispatcher calls
  `TH1::SetDefaultSumw2(kTRUE)` defensively, so any *future* histogram
  added here will be Sumw2'd by default even if the author forgets.
- **`plotting/`** — every macro that touches histograms either reads them
  from skim outputs (where Sumw2 is set) or uses `Clone()` / `Add()` /
  `Integral()`, all of which propagate the Sumw2 array correctly.
- **`analysis/charge_asym.C`, `analysis/FBratio.C`** — error handling is
  now Sumw2-aware (fixed 2026-05-29, together with the gen-weight change
  below). `analysis/analysis_helpers.h` defines a `Yield {value, error}`
  struct; `YieldInRange` returns it via `TH1::IntegralAndError` (reads the
  stored σ², not √N), and `AsymErr` / `RatioErr` take `Yield`s and do full
  linear error propagation. These reduce *exactly* to the old Poisson forms
  when the input is unweighted, so they are correct whether or not the MC
  weights turn out to be unity. `Yield::operator+` combines independent
  yields (errors in quadrature) — used for F = W⁺+W⁻ in `FBratio.C`.

## Generator event weight (added 2026-05-29)

`skim/skim.C` now applies the per-event generator weight to MC. The weight is
`hiEvtAnalyzer/HiTree::weight` (a `Float_t`; the same tree also carries
`pthat`, `ProcessID`, `ttbar_w`). All four channel functions:

- wire the branch, gated on `has_genWeight = isMC && [haveHiTree &&]
  HasBranch(tHi, "weight")` (the `haveHiTree` clause only in the Z channels,
  which tolerate a missing HiTree); warn once if an MC file lacks it;
- define a per-event `const double w = has_genWeight ? (double)genWeight : 1.0;`
  right after the `GetEntry` block — **data is always filled with weight 1**;
- pass `, w` to every one of the 19 `Fill()` calls (MT, MET, FB variants,
  QCD iso-sidebands, Z mass histos).

The cutflow `N[]` counters are intentionally left as **raw integer event
counts** (diagnostics, not yields). This is *only* the per-event gen weight —
absolute cross-section normalization (Σw, lumi·xsec) is still handled
downstream and the datacards remain shape-only. **TODO for the user:** verify
the `weight` distribution in the MC (all == 1 → weighting is a no-op; spread or
negative weights → it matters). If lepton SFs / recoil / pileup weights are
added later, fold them into the same `w`.

## MC normalization — N_gen and the absolute scale (added 2026-06-23)

To make MC comparable to data (and the samples to each other), each MC sample
`s` is scaled by `k_s = σ_s · L_int / N_gen,s`. The per-event gen weight (above)
is folded into the skim histograms; `k_s` is the *absolute* normalization layered
on top, applied **downstream** (plotting / Combine-input time, not in the skim —
keeps the skim raw and re-skim-free).

- **N_gen** — `skim/count_ngen.C` (run via `skim/run_ngen.sh`) sums
  `hiEvtAnalyzer/HiTree::weight` over **all** events (no selection) of each of the
  9 MC files → `skim/rootfile/ngen.root` (`h_ngen`/`h_nraw` + per-sample
  `TParameter<double> ngen_<label>`) and a readable `skim/output/ngen.txt`. N_gen
  is per physical file (deduped across the mu/ele/tau routes of `ResolveMCSample`).
  Its `⟨w⟩`/`N_neg` columns also serve the gen-weight verification TODO above.
- **k_s** — `skim/mc_norm.h` is the single source of truth: `kLumi_invnb = 46.5`
  (pO data lumi, nb⁻¹), `kA_O = 16` (Oxygen A-scaling, σ_pO = A·σ_NN), per-process
  `kSigma_Wp/kSigma_Wm/kSigma_DY` (nb), and `pONorm::MCScale("Wp_mu")` →
  `A·σ·L/N_gen` reading `ngen.root`. Returns 1.0 + warns if a σ is unset (safe
  no-op). One σ per process serves all three lepton-flavour files (universality);
  the σ here are the POWHEG **per-nucleon-nucleon** cross sections read straight
  from the weights (this production's `⟨w⟩ = σ`), so ×A=16 lifts them to pO.

**Decision (2026-06-23):** go absolute — ship **fixed, absolutely-normalized**
templates; the downstream fit floats only the overall rate (signal strength). This
makes absolute normalization of the *backgrounds* mandatory — a frozen background
at the wrong scale biases the signal.

**Current state (2026-06-23): σ resolved, A-scaled, both W+Z stacks drawn
ABSOLUTE.** The POWHEG per-event weight IS a cross-section weight (`⟨w⟩ = σ`), so σ
was read straight off `skim/output/ngen.txt` — σ(W⁺→ℓν)=6.376, σ(W⁻→ℓν)=5.464,
σ(DY→ℓℓ,m>50)=1.175 nb (identical across e/μ/τ → universality; ~0.9% NLO negatives
folded into Σw) — and filled into `mc_norm.h`. Because ⟨w⟩=σ, `k_s = A·σ·L/N_gen`
reduces to `A·L/N_raw`. `plotting/mtandmet.C` (W) and `plotting/dileptonpeak.C` (Z)
scale each per-sample MC histo by `MCScale(label)` **before** the W⁺/W⁻ `Add`, and
both set `ps.normBkgToData=false` so the stacks are drawn **absolute** — no area
norm (`dileptonpeak`'s own `dataInt/mcTotal` block was removed). First W→μν MET
stack validated it: composition physical (W±≫DY≫τ), absolute level ≈ data
(~4400 W + bkg vs 5103), low-MET data excess = QCD (not in MC), residual high-MET
overshoot: initially suspected uncalibrated MET/recoil, but the recoil check
(2026-07-02) came out stable — now attributed to the **unembedded MC** (no pO
underlying event in the simulation; NOT efficiency). The Z mass peak (no MET)
is the clean lepton-efficiency/scale + absolute-norm cross-check. **Style:** stacked
plots restyled to translucent fills (α=0.5) + thick color-matched outlines
(width 2); muon `h_mt_inclusive` reverted from signal-only back to the full stack.
**Escape hatch:** if absolute MC is ~16× off, set `kA_O=1.0` (lumi was per-NN L_NN).
**Untouched (still raw / shape-only):** the skim, the Combine fork's
`make_combine_input{,_Z}.C`, and the `correction/` `SaveDataMCRatio` checks.
Low-mass DY skipped. NB: per-sample skim outputs must be (re)built on the local
files (`run_all.sh`) before the plots produce anything. See README "Module 2b" and
memory `project_mc_normalization.md`.

## Pre-existing asymmetries between channels — user-confirmed status

User reviewed these on 2026-05-25:

- **DY-veto pT threshold**: 15 GeV in `skim_Wmu`, 10 GeV in `skim_Wel` —
  **intentional**, do not harmonize.
- **Isolation cut**: 0.15 in `skim_Wmu`, 0.095 in `skim_Wel` — **intentional**;
  both re-confirmed as Youden-J optima under the Δβ-corrected definition
  (2026-07-06), so the values are now settled.
- **DY-veto gate type**: ~~continuous PF relIso (muon) vs integer
  `eleCutIdWP95` / `eleMVAIsoWP95` (electron)~~ — was pre-existing; **harmonized
  2026-07-02**: both flavours now gate veto legs on continuous PF relIso
  (μ < 0.15, e < 0.095) with the flavour's tight ID (`eleMVAIdWP95` for e).
- **DY-veto mass-window comment** ("mll > 30 GeV") — **bug, FIXED in the
  refactor**. Active `skim.C:107, 192` correctly says `mll in (80, 110)`.
  The misleading comment only remains in `skim/legacy/*.C`.
- **`tHi->GetEntry(ie)` unguarded in `skim_Zmm`** — **bug, FIXED**. Now
  gated on `haveHiTree` matching `skim_Zee`'s pattern.
- **`requiredAncestorPdg`** in the gen-matching helper —
  **intentional, reserved for later use**. Keep the parameter.

## Guard audit (May 2026) — applied fixes

Added 9 FATAL blocks across the four channel functions to prevent segfaults
on files missing mandatory trees (`skimanalysis/HltTree`,
`HiGenParticleAna/hi` when `isMC`, `hltobject/HLT_*`,
`hltanalysis/HltTree`). Also: removed unused `tHi` and trigger-object-ID
parameters from helpers in `skim_common.h`; reordered `HasBranch` before
`SetBranchStatus` for PFIso branches in `skim_Zmm`; added top-of-loop
null-checks for `{mu,ele}{Pt,Eta,Phi,Charge}` in all four channel event loops.

**Deferred** (will re-audit when user is ready):
- `HasBranch` guards around electron ID/iso branch wiring (`eleCutIdWP95`,
  `eleMVAIdWP80`, etc.) in `skim_Wel`. Will FATAL on missing when applied.
- Startup-time warnings when `RelIsoPF` or `ComputePFMET` degrade to the
  null path (currently silent).

See memory `project_guard_audit.md` for the full audit + line references.
