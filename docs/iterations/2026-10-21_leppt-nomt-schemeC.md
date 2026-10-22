# Iteration 2026-10-21b/c — the no-m_T lepton-pT fit (scheme C) as THE nominal

**Status:** built, inputs written, cards validated offline; **NOT fitted**
(lxplus). Builds on `2026-10-21_noreco-nomt-study.md` (electron Reco SF →
one ±10% nuisance; the multijet design study). The last fitted iteration is
`2026-10-19_abcdy-mode.md`; its leppt_mt40 results stay valid for that model.

Two rounds the same day:
- **b** — scheme C built with ONE multijet level per rapidity group.
- **c** (this state) — after the user checked the "uniform T within a group"
  assumption: (1) the **e/γ-calibrated PF MET** in every skim (§3), (2) the
  multijet level **per rapidity bin**, (3) the residual lnN split into a
  **window** term (coherent per charge) and a **transfer** term (per group),
  (4) a set of **diagnostic plots** (§5). The b version was never fitted; its
  numbers are kept in §8 for reference.
- **c2 (2026-10-22)** — the four transfer-κ choices settled with the user
  after a dedicated study (§4b): r_W stays 1.2, κ = dev ⊕ stat, the transfer
  is **measured pooled over all 12 bins** (window m_T 30–40 μ / **30–50 e**)
  and that one κ goes to every group, each with its own θ; the group shape's
  **tail above 50 GeV is the pooled all-bins anti-iso shape** scaled to the
  group's integral.

## 1. What was asked

- *"Yes: build the no-m_T nominal with scheme C (SR split at m_T = 30 sharing
  r, multijet = T_g per group × slope × group anti-iso shape × per-bin
  anti-iso count), PF MET."*
- On the uniform-T-within-group assumption: *"Can we check this? confirm with
  me before changing."* → checked: real and significant for electrons (§4).
- *"For the MET correction, we know muons are fine … Correct matched PF
  candidates if it is an electron, using electron pT / electron raw pT …
  photon pT and photon raw pT … Keep m_T = 30 GeV as the separation."*
- *"Since we know the T is not uniform within group, why don't we do a
  per-rapidity-bin shape, or a per-rapidity-bin nuisance?"* → agreed design:
  per-group shape, per-bin free level, window lnN coherent + transfer lnN per
  group (*"I agreed, let's do it"*), plus *"diagnostic plots, just like what
  we have for abcd before … so we can verify these potential issues."*

## 2. The model (c)

- **Discriminant** `leppt` (fork `--disc leppt`, default; out-tree
  `pO_fit_out_leppt/`): iso-pass leading-lepton pT 24–100 (2 GeV bins),
  pT > 25, **no m_T cut**. Each (flavour, charge, rapidity bin) is **two
  channels**, `<F>_<C>_<B>_y<i>_lo` (m_T < 30 = region B) and `_hi`
  (m_T ≥ 30), sharing the bin's r. MET (now e/γ-calibrated, §3) enters only
  through the lo/hi split.
- **Multijet of bin y, channel c** = T_y × exp(s (pT − 40)/10) × the GROUP's
  anti-iso (QCD-only) pT shape of channel c × the bin's anti-iso count
  N_{y,c}/N_{g,c} (`correction/qcd_abcd_nomt.h runNoMt`). Groups fwd y10–11 |
  bwd y0–1 | cen y2–9 (lab and fb) for the SHAPE and the transfer term only.
- **Fit parameters** (fork `QCD_MODE=mtsplit`): per (flavour, charge, BIN) a
  free level `qcd_sT_<F>_<C>_y<i>` on both channels (prefit 1 = T_y), and the
  anti-iso anchors `qcd_sD_<F>_<C>_y<i>` (lo + the 1-bin CR
  `<F>_<C>_<B>_y<i>_CRD`, anti-iso m_T < 30) / `qcd_sH_…` (hi + `_CRH`,
  anti-iso m_T ≥ 30); CR `ewk` frozen (lumi only). The lo channel (77–84%
  multijet) measures sT in the fit; W is only ~2% of the W MC there, so sT is
  not degenerate with r.
- **Electron slope**: one shared s per charge (all bins), from the m_T < 30
  profile likelihood with every bin its own level; shape nuisance
  `qcdFF_ele_<C>` at ŝ ± 2 × its 68% interval. Muon: flat (s = 0).
- **Residual lnN on the m_T ≥ 30 multijet only**, two terms:
  - **window** `qcd_rate_<F>_<C>` = the anti-iso window [0.3,0.6) / [0.6,1.0)
    moving the predicted hi multijet with B fixed — one method choice, so
    COHERENT over the 12 bins of a charge (sidecar `#! kappaWin <C> <k>`);
  - **transfer** `qcd_xfer_<F>_<C>_<g>` = F(m_T 30–X)/F(m_T < 30) with
    X = 40 (muon) / 50 (electron), **measured on all 12 bins pooled**, W MC at
    r_W = 1.2 (`kNoMtXferR`; r_W = 1 printed too), κ = 1 + √(dev² + stat²);
    the SAME κ for the three groups, but one UNCORRELATED θ per group
    (sidecar `#! kappaXfer <B> <g> <C> <k>`, identical k per group). Since
    c2; in c it was each group's own 30–40 value. The groups' own values are
    printed as diagnostics with their χ² against the pooled one.
- **Tail (since c2):** above pT = 50 GeV each group's anti-iso shape is
  replaced by the pooled all-bins shape, scaled to the group's own integral
  above 50 (normalization and everything below 50 unchanged).
- **Parameter count (grand, per binning):** 25 POIs; 144 free multijet scales
  (12 bins × 2 charges × 2 flavours × sT/sD/sH; b had 36); 28 constrained:
  lumi, ttbar_xsec, 4 window, 12 transfer, 10 shape (nPDF, qcdScale, alphaS,
  muSF, eleSF, eleReco, metResMu, metResEle, qcdFF_ele_Wp/Wm).
- **Channels:** 96 W SR + 2 Z + 96 anchor CRs = 194 (grand), 97 per flavour.
- **Fiducial:** the gen set WITHOUT m_T (`gen_xsec_fid.txt`),
  `pODisc::GenFid("leppt")` = tag "".

## 3. The e/γ-calibrated PF MET (c)

In DATA the forest electron (and photon) energy carries a data-only scale
correction (median elePt/eleRawPt 1.021 barrel / 1.094 endcap) that the PF
candidates do not, so the PF MET of electron events is biased along the
electron (`correction/edge_study/`). Fix (user recipe):

- the PF electron (type 2) closest to each reco electron within ΔR < 0.05 is
  scaled by elePt/eleRawPt; the PF photon (type 4) closest to each reco photon
  by phoEt/phoRawEt; ratio window [0.5, 2] (outside: not applied, counted);
  each PF candidate corrected at most once (electrons first); data and MC
  alike; all four skim channels + the two replicating macros that use MET
  (`trig_eff_mb.C` mt40, `njet_WZ.C` mt40).
- Single source `skim_common.h`: `ComputePFMETEGCalib`,
  `BindEGammaMetBranches`, `MetCalibStats` (`[MET-EG]` block in every job log);
  `PO_MET_EGCALIB=off` restores the raw PF MET (verified bit-identical to the
  previous skims).
- Muons: no momentum correction exists, so nothing changes for them except
  through photons in their events.

**Effect (skim logs, `[MET-EG]`):**

| sample | electrons with a PF match | ⟨elePt/eleRawPt⟩ | photons matched | ⟨phoEt/phoRawEt⟩ | ⟨\|MET\|_calib − \|MET\|_PF⟩ |
|---|---|---|---|---|---|
| W→eν data (step 8) | 81.2% | 1.049 | 22.5% | 1.055 | +0.53 GeV |
| W⁺→eν MC | 96.5% | 0.997 | 5.7% | 0.997 | −0.11 GeV |
| Z→ee data (peak) | 97.2% | 1.040 | 4.9% | 1.035 | +0.43 GeV |
| W→μν data | — | 1.007 | 79.5% | 1.035 | +0.02 GeV |

(The lower data match rate is the multijet: a fake electron is a PF hadron.)
Data W→eν candidates: m_T > 40 e⁻ 1826 → 1911 (+4.7%), e⁺ 2332 → 2435
(+4.4%); |η| > 2.0: e⁻ 301 → 345 (+14.6%), e⁺ 419 → 459 (+9.5%); muons
unchanged (2751 → 2748). Prefit leppt_mt40 e/μ (`prefit_r`, lab): W⁻
inclusive 0.859 → 0.906, **y0 0.571 → 0.752 ± 0.130, y11 0.597 → 0.729 ±
0.131**, W⁺ inclusive 0.964 → 1.016 — the electron W⁻ edge deficit shrinks
but does not close. MET resolution from the Z recoil (`recoil_res`): k_e
1.2645 → 1.2571, k_μ 1.1030 → 1.1037 — the Z→ee recoil-width excess is NOT
this bias (the two electrons cancel). Trigger SF (mt40): e 1.0010 → 1.0011,
μ unchanged.

## 4. Numbers (prefit, lab unless noted)

**Why per-bin levels — within-group LRT (one level per group vs per-bin T_y,
at the shared slope):** e⁺ fwd p = 0.006, bwd 0.047, cen 0.088; e⁻ fwd 0.030,
bwd 0.026, cen 0.072; μ⁺ 0.61 / 0.65 / 0.25; μ⁻ 0.71 / 0.72 / **0.0055**
(cen). E.g. e⁺ fwd y10 0.727 [0.634, 0.833] vs y11 1.272 [1.090, 1.483];
e⁻ bwd y0 1.212 vs y1 0.749.

**T_y (m_T < 30, 68%):**

| | y0 | y1 | y2 | y3 | y4 | y5 | y6 | y7 | y8 | y9 | y10 | y11 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| μ⁺ | 0.236 | 0.282 | 0.162 | 0.168 | 0.113 | 0.185 | 0.086 | 0.142 | 0.109 | 0.204 | 0.237 | 0.275 |
| μ⁻ | 0.219 | 0.258 | 0.279 | 0.251 | 0.137 | 0.094 | 0.092 | 0.140 | 0.201 | 0.162 | 0.162 | 0.185 |
| e⁺ | 1.031 | 0.649 | 0.400 | 0.354 | 0.413 | 0.382 | 0.145 | 0.361 | 0.383 | 0.571 | 0.727 | 1.272 |
| e⁻ | 1.212 | 0.749 | 0.500 | 0.550 | 0.287 | 0.179 | 0.257 | 0.317 | 0.370 | 0.457 | 0.683 | 1.086 |

(e: at the shared slope; intervals in `correction/logs/qcd_abcd_<lep>.log`,
lines `[NOMT-T]`, and in the `nomt_T_y_*` plots.)

**Slope:** e⁺ +0.257 [+0.181, +0.333], e⁻ +0.281 [+0.204, +0.357] per 10 GeV;
flat rejected (p = 0.0014 / 0.0005); one slope for all groups fine (p = 0.57 /
0.60). Muon shared slope (diagnostic only): μ⁺ +0.174 [+0.078, +0.260], flat
p = 0.096; μ⁻ +0.065, p = 0.62 → flat kept. (μ⁺ bwd's own slope runs to the
scan edge −1.5: two bins, low statistics.)

**Residual κ (round c: transfer per group from m_T 30–40 — superseded by the pooled values of §4b):**

| | window | xfer lab fwd / bwd / cen | xfer fb fwd / bwd / cen |
|---|---|---|---|
| μ⁺ | 1.069 | 1.581 / 1.517 / 1.217 | 1.736 / 1.564 / 1.218 |
| μ⁻ | 1.037 | 1.978 / 1.631 / 1.270 | 1.887 / 1.349 / 1.226 |
| e⁺ | 1.032 | 1.154 / 1.227 / 1.278 | 1.263 / 2.085 / 1.164 |
| e⁻ | 1.041 | **1.666** / 1.188 / 1.216 | 1.603 / 1.546 / 1.354 |

The muon transfer κ's are mostly the statistical error of F(m_T 30–40)
(±0.2–0.55: little multijet there). The electron ones are real deviations in
places: **e⁻ fwd F(30–40)/F(<30) = 1.62 ± 0.24 (2.6σ)**, e⁺ fb bwd 1.98 ±
0.48. For the muon the hi multijet is 7–9% of the data, so a large κ costs
little.

**SR composition (`[NOMT-IN]`, lab sums):**

| | data lo | multijet lo | data hi | W hi | multijet hi | data/prefit hi |
|---|---|---|---|---|---|---|
| μ⁺ | 368 | 285.1 (77%) | 2877 | 2089.0 | 199.5 (7%) | 1.166 |
| μ⁻ | 344 | 276.3 (80%) | 2410 | 1665.9 | 212.7 (9%) | 1.181 |
| e⁺ | 610 | 507.1 (83%) | 2713 | 1548.1 | 696.1 (26%) | 1.150 |
| e⁻ | 592 | 500.1 (84%) | 2194 | 1237.2 | 695.2 (32%) | 1.076 |

The calibration moved ~90 electron data events per charge from lo to hi
(b: e⁺ 693 / 2630, e⁻ 686 / 2100). If the electron r equalled the muon's
(hi: 1.196 μ⁺, 1.222 μ⁻), the hi multijet would need ×1.07 (e⁺) and
**×0.83 ± 0.07 (e⁻)** — was ×0.74 in b.

## 4b. The transfer-κ choices (c2, 2026-10-22)

Scratch study (`xfer_options.C` in the session scratchpad, read-only),
then implemented in `qcd_abcd_nomt.h`. Lab binning.

| choice | evidence | adopted |
|---|---|---|
| r_W of the EWK subtraction | e: W 4–18% of the iso-pass 30–40 data, ratio moves ≤ 0.01 between r_W 1.0 and 1.2; μ: W 25–34%, moves ~0.1 (μ⁺ 1.08 → 0.98); the μ m_T ≥ 30 data prefer r ≈ 1.20 / 1.22 | **1.2** |
| κ recipe | pooled: dev ⊕ stat vs max(dev, stat) differ by ≤ 0.03 | **dev ⊕ stat** |
| window | e: W 14–37% still at 40–50, twice the events, nearer ⟨m_T⟩ ≈ 48 GeV of the m_T ≥ 30 multijet; e⁻ fwd 1.62 ± 0.24 (30–40) → 1.28 ± 0.14 (30–50); μ: W 50–57% at 40–50 and r_W-dependent | **30–40 μ, 30–50 e** |
| per group vs pooled | groups agree with the pooled value (30–40 μ⁺ p 0.62, μ⁻ 0.09; 30–50 e⁺ 0.92, e⁻ 0.21); per group the μ value is statistics (±0.2–0.55) | **pooled κ, one θ per group** |
| forward spikes | every group's tail agrees with the pooled tail (m_T ≥ 30: p 0.21–0.96); e⁺ fwd pT ≥ 60 multijet 16.0 ± 4.6 → 16.0 ± 2.6; e⁺ bwd 1.4 → 9.3 (data 31 vs W + other 10.0), e⁻ bwd 6.9 → 11.8 (33 vs 13.9) | **pooled tail above 50 GeV** |

**Transfer κ now** (`[NOMT-XFER] … POOLED`): lab μ⁺ 1.168 (0.983 ± 0.167),
μ⁻ 1.148 (0.995 ± 0.148), e⁺ 1.183 (1.167 ± 0.073), e⁻ 1.104 (1.078 ±
0.069); fb μ⁺ 1.177, μ⁻ 1.158, e⁺ 1.204, e⁻ 1.162. Window κ unchanged
(e⁻ 1.041 → 1.042 through the templates). The m_T < 30 tail tests in two
groups (e⁻ bwd p ≈ 0, μ⁺ fwd p 0.01) rest on 2–5 anti-iso events, where a
Pearson χ² means nothing; those tails hold ≤ 2% of the group's multijet.
Multijet totals: μ unchanged (199.5 / 212.7 hi), e hi 696.1 → 697.3,
695.2 → 697.3. Validation: every other Combine input and sidecar
bit-identical (`compare_hists.C` vs `backups/pre_xferpool/`); in
`combine_input_W_leppt.root` only `qcd` (145 μ / 147 e) and the electron
`qcd_qcdFF_ele_*` twins changed; the cards differ from c only in the κ values
(maps identical), parse with `DatacardParser`, 0 missing shapes. New plots
`nomt_tail_<tag>_<g>_<lo|hi>`; `nomt_xfer_*` shows the pooled band.

## 5. Diagnostic plots

`correction/plots/qcd_abcd_<lep>/nomt_fit/` (made by `run_qcd_abcd.sh`):

| plot | question |
|---|---|
| `nomt_regions_<lep><chg>` | the (iso, m_T) plane with lo (= B), hi, D, H, the gap |
| `nomt_T_y_<tag>` | T_y per bin vs the one-level-per-group bands + the LRT per group |
| `nomt_lofit_<tag>` | m_T < 30 data/prefit with the group level vs the per-bin level; what the m_T ≥ 30 data need at W × 1.2 |
| `nomt_ff_<tag>` | F(pT) per group vs the shared slope (+ each group's own slope, LRTs) |
| `nomt_xfer_<tag>` | F in m_T slices / F(m_T < 30) per group at r_W = 1 and 1.2 (the transfer κ) |
| `nomt_shape_<tag>_<g>_<lo\|hi>`, `nomt_shape_p_<tag>` | each bin's anti-iso shape vs its group's (χ², p) |
| `nomt_mt_<tag>_<g>` | iso-pass m_T per group: data vs W (r = 1) + other + multijet (T_y f(pT) × anti-iso) |
| `nomt_sr_<tag>_<g>_<lo\|hi>` | per-group pT stacks of the two SR channels (prefit, W at r = 1) |

(`<tag>` = `<lep><Plus|Minus>_<lab|fb>`.) Findings so far: per-bin anti-iso
shapes agree with their group's everywhere (all p > 0.05 in e⁺ lab), so a
per-group shape is adequate; the electron m_T ≥ 30 multijet template has
**spikes at pT > 60 in the forward group** (single anti-iso events × f(pT) ≈
3 × T_y), next to W-dominated bins — template noise the fit has no nuisance
for (no autoMCStats).

## 6. What was built (b + c)

| file | change |
|---|---|
| `skim/skim_common.h` | (c) `ComputePFMETEGCalib`, `BindEGammaMetBranches`, `MetCalibStats`, `PO_MET_EGCALIB` |
| `skim/skim.C` | (b) `h_leppt_{mtlo,mthi}_*` + twins; (c) e/γ MET in all four channels + `[MET-EG]` |
| `correction/trig_eff_mb.C`, `correction/njet_WZ.C` | (c) e/γ MET for the mt40 variants |
| `correction/qcd_abcd_nomt.h` (new, c) | `runNoMt` (per-bin level, shared slope, window + transfer κ, the plots of §5); `qcd_abcd.C` includes it and passes `plots/qcd_abcd_<lep>/nomt_fit` |
| `plotting/combine_input_nomt.C` | (b) the input; (c) per-bin CR dirs `<C>_<B>_y<i>_CR{D,H}`, sidecar `#! kappaWin` / `#! kappaXfer` |
| fork `make_pO_simfit_cards.sh` | (b) `QCD_MODE=mtsplit`; (c) per-bin sT/sD/sH + CRs, window κ on the 4 `qcd_rate` rows, 12 `qcd_xfer_*` rows, env `QCD_MTSPLIT_WIN_MU/ELE`, `QCD_MTSPLIT_XFER_MU/ELE`, sidecar `mtsplitKappas` |
| fork `run_pO_fits.sh` | (c) passes `mtsplitKappas` to the extractor (13th arg); postfit per-bin level name |
| fork `extract_pO_simfit.C` | (c) `MtSplitScale` per bin (lo = sT·sD, hi = sT·sH·κ_win^θ·κ_xfer^θ), rows `qcd_s{T,D,H}_<F>_<C>_y<i>`, `qcd_mtsplit_{lo,hi}_mult_*_y<i>`, `qcd_xfer_*_theta`, Asimov closure on all |
| `plotting/r_per_bin.C`, `compare_fits.C`, `postfit_incl.C` | mtsplit per-bin multipliers, lo + hi sums, ndf |

## 7. Validation (c)

- `PO_MET_EGCALIB=off`: every skim histogram bit-identical to the previous
  skims (before the full re-skim).
- Re-skim (`skim/logs/compare_metcalib.log` vs `backups/pre_metcalib/`): data
  files differ ONLY in the MET-dependent histograms (Z data: the 4 recoil
  histograms; W data: 391 of 508 e, 262 of 508 μ); every MC histogram moved,
  through the tiny trigger-SF change (Z→μμ ⟨TRIG⟩ identical to 5 digits,
  Σw·SF −3e-8) and the MET.
- LHE Up/Down: `max|member0 − nominal| = 0` in all 12138 templates.
- `combine_input_nomt`: every bin's SR multijet template sums to the
  `nomt_counts_*_y<i>` value (no WARN); lo data/prefit 0.99–1.00.
- Fork dry-run (`both all --disc leppt`): 194 / 97 / 97 channels; every card
  parses with `DatacardParser`; 778 (grand) nominal + 5466 Up/Down shapes
  resolve, 0 missing; 144 distinct rateParams; `qcd_rate_*` on the 12 hi qcd
  columns of the charge, `qcd_xfer_*` on the 2 / 2 / 8 hi qcd columns of the
  group; maps: every W SR signal/wtau column its own r (0 mismatches).
- Extraction on a synthetic fit_s (r 1.15, sT 0.90, sD 1.00, sH 1.05, every
  θ = −0.5): lo 0.9000, hi = 0.945 × κ_win^−0.5 × κ_xfer^−0.5 (e.g. μ⁺ y0
  0.7421, e⁻ y11 0.7176 — by hand 0.7421 / 0.7177); 169 stat parameters kept,
  28 conditioned; `postfit_incl` and `r_per_bin` run on it (fake outputs moved
  to the session scratchpad).

## 8. The b version (superseded, never fitted)

One level per group F0_g (e.g. e⁺ fwd/bwd/cen 0.975 / 0.788 / 0.351), slopes
e⁺ +0.212, e⁻ +0.321, one κ_hi per flavour (μ 1.080, e 1.168, window ⊕
F(30–40)/F(<30) at r = 1), 24 group CRs, 122 channels, raw PF MET.

## 9. Not done

- The fit (lxplus): `./sync_lxplus.sh upload` + `upload-scripts`, then
  `./run_pO_fits.sh both all --asimov` (default `--disc leppt`, QCD mtsplit),
  `./run_pO_impacts.sh --disc leppt --fit all`, `./sync_lxplus.sh download`,
  `analysis/run_observables.sh leppt`.
- (The open κ choices of c were settled in c2, §4b.)
- The other fits' inputs (leppt_mt40, met_mt40, met, Z) now also carry the
  calibrated MET; their last fits used the raw PF MET.
- `run_syst_shapes.sh` has no `leppt` mode; the AN describes the leppt_mt40
  abcdy fit.

## 10. Backups

`backups/pre_xferpool/` = the c state before c2: `correction/rootfile/qcd_abcd_{mu,ele}.root`,
`correction/logs/qcd_abcd_{mu,ele}.log`, `correction/plots/qcd_abcd_{mu,ele}/nomt_fit/`,
`plotting/plots[/Elec]/combine_input_*` + sidecars, `plotting/logs/*.log` (50 MB).
`backups/pre_metcalib/` = the state before the e/γ MET re-skim (round b
inputs): `skim/{rootfile,logs,output}`, `correction/{rootfile,logs}`,
`correction/{electron,muon}_sf_check/trigger/`, `plotting/{plots (the Combine
inputs + sidecars),logs}` — 390 MB. The b-version `runNoMt` was never
committed (it lived in the uncommitted `qcd_abcd.C`; replaced by
`qcd_abcd_nomt.h`). The step-1 backups of the day are `backups/pre_noreco/`.
