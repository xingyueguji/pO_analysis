# Iteration 2026-10-19 — rapidity-resolved in-fit ABCD (`QCD_MODE=abcdy`)

> **Paths (2026-10-20):** every backup below now lives in `backups/pre_<tag>/<area>/`
> (`<area>/pre_<tag>/` → `backups/pre_<tag>/<area>/`, `docs/AN-26-030_pre_<tag>/` →
> `backups/pre_<tag>/docs/AN-26-030/`, fork `test/pO_fit_out<sfx>_pre_<tag>/` →
> `test/backups/pre_<tag>/pO_fit_out<sfx>/`). Full list: `backups/MOVES_2026-10-20.txt`.

**Status:** built locally on 2026-10-19 as an **opt-in** mode. **Fitted on
lxplus the same day**; its observables and the comparison with the nominal
abcd fit are in §6. The nominal mode (`QCD_MODE=abcd`) and every
input object of the pending nominal fit (`2026-10-19_fbqcd.md`) are unchanged,
so that record and its fit stay valid. Backups: tag `pre_abcdy_mode`, §7.

## 1. What the mode is

User decision, after the step-1 test of the multijet rapidity split
(`correction/qcd_ysplit.C`, CLAUDE.md "Active"): run the in-fit ABCD
**separately for three rapidity groups**, using the same bin indices for the
lab and the FB binning:

| group | bins | |
|---|---|---|
| `fwd` | y10, y11 | group 1, most forward (proton-going, y = −η > 0) |
| `bwd` | y0, y1 | group 2, most backward |
| `cen` | y2 … y9 | group 3 |

The groups are defined in one place, `correction/abcd_ygroups.h`.

- Each group g gets its own control regions (B: iso-pass, m_T < 30; C40:
  anti-iso, m_T > 40; D: anti-iso, m_T < 30). They come from the y-binned
  skim planes `h_iso_mt_<lep><chg>_y<i>[_FB]`. Each group also gets its own
  free scales `qcd_s{B,C,D}_<F>_<C>_<g>` and its own transfer, so its SR
  multijet is A_g = B_g·C40_g/D_g × sB_g·sC_g/sD_g.
- Three things stay inclusive: the split *inside* a group (the low-MET excess
  proxy), the lepton-pT shape, and the electron fake-factor weighting.
- The residual lnN `qcd_rate_<F>_<C>` stays one per (flavour, charge), shared
  by the groups (κ μ 1.09 / e 1.05, as in abcd).
- SR template, per bin y of group g:
  `qcd_abcdy_y = qcd_abcd base × (A_g/A0) × w_in,y`, where A0 = B0·C40/D0 is
  inclusive and w_in,y = the bin's low-MET excess divided by the excess of its
  group, in the same binning. Each group therefore carries A_g; for the
  electron it carries A_g × A0ff/A0.
- The FB groups use the FB-binned planes, i.e. the FB window's own control
  regions. Their totals are therefore the window's share by construction.
- This mode was not expected to resolve the electron W⁻ charge split.

## 2. Code

**Analysis repo**
- `correction/abcd_ygroups.h` (new): the groups.
- `correction/qcd_abcd.C`: new function `runGroupABCD`, printed last. It
  writes the per-group counts as `abcd_counts_<lep><chg>_<lab|fb>_<g>`, with
  the inclusive label set. `T_met` is set to 0, and for the electron the
  `A0ff*` labels are the inclusive values × A_g/A0. The log adds a closure
  line: the lab groups' data must equal the inclusive plane exactly.
- `plotting/mtandmet.C`, additions only, `leppt_mt40` file:
  - `qcd_abcdy` in every SR region and in the incl dirs;
  - the electron slope twins `qcd_abcdy_qcdFF_ele_{Wp,Wm}Up/Down`;
  - 36 group CR dirs `{Wp,Wm}_{lab,fb}_{fwd,bwd,cen}_CR{B,C,D}`;
  - the electron sidecar rows now read `qcdFF_ele_Wp qcd_abcd@Wp qcd_abcdy@Wp`;
  - the group sidecar `combine_input_W_leppt_mt40_abcdy.txt`, written only
    when everything above is complete (otherwise it is deleted).
- `plotting/postfit_incl.C`: `qcd_model abcdy` reads the `qcd_abcdy` template;
  the ndf counts one QCD parameter per group.
- `correction/qcd_ysplit_pergroup.C` + `run_qcd_ysplit_pergroup.sh` (new): the
  before/after diagnostics.

**Fork (`zheng/po-analysis`, `test/`)**
- `my_script/make_pO_simfit_cards.sh`, `QCD_MODE=abcdy` (leppt_mt40 only):
  - reads the group sidecar and refuses to run without it;
  - per group: CR channels `<F>_<C>_<g>_CR*`, scales, and the formula
    `qcd_abcd_<F>_<C>_<g>` repeated once per SR channel of the group;
  - QPATH `qcd_abcdy`;
  - the group-channel `w_y` maps, both nominal and σ-POI;
  - the κ sidecar gets `qcdMode abcdy` and a `qcdGroups` line;
  - side change: an invalid `QCD_MODE` value is now an error instead of
    silently becoming lnN.
- `my_script/extract_pO_simfit.C`, new 12th argument `qcdGroups`:
  - per-group `AbcdScale`;
  - the CSV qcd columns carry each bin's group multiplier;
  - `qcd_model abcdy`;
  - summary rows per group;
  - the Asimov closure covers every group scale.
- `run_pO_fits.sh`: copies the group sidecar, passes `qcdGroups` to the
  extractor, and defaults the out-tree to `pO_fit_out_leppt_mt40_abcdy/`.
- `run_pO_impacts.sh`: `QCD_MODE=abcdy` points it at the `_abcdy` tree.
- `sync_lxplus.sh`: uploads the group sidecars and downloads the `_abcdy` tree.

## 3. Validation

| check | result |
|---|---|
| `run_qcd_abcd.sh both` vs backup (`compare_hists.C`) | IDENTICAL 12 / 18, DIFF 0, NEW 12 per file; logs = pure insertion of the abcdy block |
| lab groups' data vs inclusive B / C40 / D | EXACT, all four flavour × charge |
| six Combine inputs vs `plotting/plots/pre_abcdy_mode/` | IDENTICAL 2397 / 2520 / 46 (μ met / leppt_mt40 / Z), 2397 / 2622 / 46 (e); DIFF 0, MISSING 0; NEW = 51 `qcd_abcdy` + 192 group-CR objects (μ), + 102 slope twins (e) |
| sidecars | all identical except the electron leppt_mt40 `_systs.txt`, where the two qcdFF lines gained `qcd_abcdy@<C>` |
| writer logs | only the "Saved … + 36 abcdy group CR dirs" line and 15 `[QCD-ABCDY]` lines new; dileptonpeak, prefit_r logs identical |
| fork dry-runs lnN / free / abcd / default (leppt_mt40 + met, `both all`) | 105 card / map / sidecar files identical to the pre-change ones up to the path; abcd also = the current local cards |
| abcdy cards, offline harness (DatacardParser + ModelTools + MultiSignalModel, real RooWorkspace) | 6/6 cards pass. 86 / 43 channels; 442 / 221 aligned columns; 2490 / 1221 / 1269 shape histograms resolve; the qcdFF rows sit only on `qcd_abcdy` of their own charge (12 entries each); ONE RooFormulaVar per group, shared by its 2/8/2 SR channels, evaluating sB·sC/sD; every column has the intended scale with the nominal and the σ-POI maps; r_i = 1 at init |
| `QCD_MODE=abcdy --disc met` | refused by the card generator, as intended |
| extractor, default arguments, on the 10-05 fit tree (`pO_fit_out_leppt_mt40_pre_ttbar`) | CSVs byte-identical, 104 histograms identical (grand, μ-only, e-only) |
| extractor, abcdy path, synthetic fit_s | per-bin group multipliers = κ^θ·sB·sC/sD of the group; 24 group multiplier rows; Asimov closure over 72 group scales |
| `postfit_incl.C` on the synthetic abcdy CSV | runs on `qcd_abcdy` (mechanics only) |
| diagnostics self-checks | Σ_{y∈g} qcd_abcdy = A_g × template/A0 to ≤ 7e-16; group CR dirs = the exported counts; lab data closure EXACT |

## 4. Prefit (leppt_mt40, `correction/logs/qcd_ysplit_pergroup.log`)

### Group multijet totals

OLD = `qcd_abcd` summed over the group's bins; NEW = `qcd_abcdy` ± the group
ABCD stat error.

**Lab binning**

| | fwd old → new | bwd old → new | cen old → new | total |
|---|---|---|---|---|
| μ⁺ | 31.8 → 23.4 ± 4.3 | 16.1 → 17.3 ± 4.0 | 80.3 → 84.1 ± 8.5 | 128.2 → 124.8 |
| μ⁻ | 22.9 → 16.9 ± 3.4 | 13.7 → 16.2 ± 4.3 | 87.3 → 91.3 ± 8.7 | 123.9 → 124.4 |
| e⁺ | 159.9 → 154.6 ± 18.6 | 128.3 → 104.6 ± 14.1 | 185.5 → 189.0 ± 18.3 | 473.7 → 448.2 |
| e⁻ | 164.0 → 152.4 ± 18.6 | 150.4 → 140.1 ± 17.5 | 187.1 → 180.4 ± 17.9 | 501.6 → 472.8 |

**FB binning**

| | fwd | bwd | cen | total |
|---|---|---|---|---|
| μ⁺ | 27.1 → 17.5 | 16.4 → 16.6 | 73.1 → 74.7 | 116.6 → 108.8 |
| μ⁻ | 21.8 → 15.9 | 21.2 → 26.1 | 69.1 → 70.2 | 112.2 → 112.2 |
| e⁺ | 141.5 → 122.5 | 53.1 → 59.1 | 173.8 → 175.3 | 368.5 → 356.8 |
| e⁻ | 136.2 → 125.2 | 54.1 → 67.2 | 183.4 → 165.2 | 373.6 → 357.7 |

### Group transfer factors

Relative to the inclusive values (lab). Electron T_g/T: 1.34–1.65 at the edges
and 0.58–0.63 central. R_g/R: 0.81–0.92 at the edges and 1.10–1.11 central
(the step-1 picture, now per group). Muon T_g/T: 1.03–1.53 at the edges.

### e/μ of the prefit r_Q = (D − B − Q)/W (lab, data stat. only)

| | W⁺ old → new | W⁻ old → new |
|---|---|---|
| y0+y11 | 0.943 → 0.979 | 0.568 → 0.618 ± 0.088 |
| y0..y3 | 1.044 → 1.094 | 0.732 → 0.763 |
| fwd group | 0.894 → 0.893 | 0.709 → 0.746 |
| bwd group | 1.008 → 1.099 | 0.665 → 0.715 |
| cen group | 1.087 → 1.086 | 0.978 → 0.988 |
| inclusive | 1.036 → 1.050 | 0.883 → 0.904 |

Per bin and fb in the log; figures and CSV in
`correction/plots/qcd_ysplit_pergroup/`.

## 5. How to fit (lxplus)

On the Mac:

```bash
test/sync_lxplus.sh upload
```

This ships the regenerated inputs, the two group sidecars, and the changed
fork scripts. Then on lxplus, in the fork's `test/` after `cmsenv`:

```bash
QCD_MODE=abcdy PO_PLOTS=<ana>/plotting/plots ./run_pO_fits.sh both all --asimov
QCD_MODE=abcdy ./run_pO_impacts.sh --fit all
```

The first command writes `pO_fit_out_leppt_mt40_abcdy/`. Then run
`sync_lxplus.sh download` on the Mac; it pulls the `_abcdy` tree.

What to check:
- Asimov PASS, including all 36 group scales.
- The 12 group multipliers `qcd_abcd_mult_<F>_<C>_<g>` and the
  `qcd_rate_<F>_<C>` pulls, against the abcd fit of the same inputs
  (`2026-10-19_fbqcd.md`).
- σ, A_ch and R_FB against the abcd fit.

The observables chain (`run_observables.sh`, `postfit_incl.C` defaults) reads
the nominal tree. For the abcdy tree, pass its CSV / fitDiagnostics paths
explicitly, or extend `plotting/fit_variants.h`. That is a follow-up.

## 6. Fit

**Run:** on lxplus on 2026-10-19, with the §5 commands
(`QCD_MODE=abcdy … run_pO_fits.sh both all --asimov` and
`QCD_MODE=abcdy ./run_pO_impacts.sh --fit all`). It was downloaded the same
evening to the fork's `test/pO_fit_out_leppt_mt40_abcdy/{simfit,simfit_mu,simfit_ele}/`.
The κ sidecar reads `qcdMode abcdy`, `kQcdMu 1.09`, `kQcdEle 1.05`, and
`lheSysts nPDF,qcdScale,alphaS,muSF,eleSF,qcdFF_ele_Wp,qcdFF_ele_Wm`.

**Downstream:**
- Observables: `analysis/run_observables.sh leppt_mt40_abcdy`, a new opt-in tag
  (see the end of this section). Its log is
  `analysis/logs/run_observables_leppt_mt40_abcdy.log`.
- Comparison with the nominal abcd fit: `plotting/run_compare_fits.sh`. Its log
  is `plotting/logs/compare_fits_leppt_mt40_vs_leppt_mt40_abcdy.log`, and its
  outputs (`compare_fits.csv`, the per-group multijet plots, the A_ch / R_FB
  overlays) are in `plotting/plots/compare_fits/leppt_mt40_vs_leppt_mt40_abcdy/`.

Every number below was read from those logs, from the fork's summary CSVs and
extraction logs, or from the `xsec_*` CSVs.

**Reference caveat — lab is clean, fb is not.** The reference is the nominal
abcd fit in `pO_fit_out_leppt_mt40/`, i.e. the **2026-10-16 fit**. The nominal
abcd model was not refit on the FB-window-fixed inputs; only the abcdy fit was
made on them. **User decision, 2026-10-19: the comparison is taken as it is,
abcd (10-16) vs abcdy (10-19).**

- Lab: the FB-window fix left the lab inputs bit-identical, so the lab
  comparison is a pure model comparison on the same inputs.
- fb: the reference fb fit used the pre-fix templates. `compare_fits.C` checks
  this from each fit's own `shapes_prefit`: the fb/lab prefit multijet total is
  **1.0000** in every (flavour, charge) for the reference. For the abcdy fit it
  is μ 0.872 / 0.902 and e 0.796 / 0.756 (W⁺ / W⁻). **Every fb difference below
  is therefore the combined effect of the FB-window fix and the per-group
  ABCD.** It is quoted that way. Separating the two would need an abcd fit on the
  fixed inputs; `plotting/run_compare_fits.sh` would pick that up with no code
  change.

### 6.1 Fit quality

| fit | status / covQual (lab, fb) | Asimov | statonly vs conditioned (max \|ΔrErr\|/rErr) | r_Z lab: abcd → abcdy |
|---|---|---|---|---|
| μ + e | 0/3, 0/3 | PASS, PASS | 0.0001, 0.0001 | 1.1589 → 1.1529 ± 0.0592 |
| μ only | 0/3, 0/3 | PASS, PASS | 0.0002, 0.0002 | 1.0702 → 1.0701 ± 0.0647 |
| e only | 0/3, 0/3 | PASS, PASS | 0.0005, 0.0007 | 1.3551 → 1.3475 ± 0.1036 |

- Asimov closure: max |POI − 1| = 0.0000 in all six fits. The abcdy
  extractor's closure covers every group scale.
- Profiled (σ_W, σ_Z) contours: all six are drawn. The staleness gate passes
  (scan minimum vs covariance best fit: |Δσ_W| ≤ 0.0053 nb, |Δσ_Z| ≤ 0.0027 nb).
- Ellipse vs profiled 1σ radius: −3.8 … +3.3 % (μ + e, lab), up to
  −4.9 … +4.7 % (e only, fb). The abcd fit gives the same.
- Impacts were run for all three fits (`*/impacts/`). They are not summarized here.

### 6.2 SR multijet per rapidity group (lab, events)

The yields are summed from each fit's own `shapes_prefit` / `shapes_fit_s`
(process `qcd`). Multipliers: abcd has one per charge, abcdy one per group.

| μ + e fit | bwd: abcd → abcdy | cen | fwd | total | abcd mult. | abcdy mult. bwd / cen / fwd |
|---|---|---|---|---|---|---|
| μ⁺ | 14.5 → 15.4 | 72.1 → 73.7 | 28.5 → 22.1 | 115.0 → 111.1 | 0.897 | 0.890 / 0.875 / 0.945 |
| μ⁻ | 13.0 → 17.7 | 83.0 → 85.5 | 21.8 → 16.0 | 117.8 → 119.2 | 0.951 | 1.095 / 0.936 / 0.947 |
| e⁺ | 117.3 → 97.4 | 169.6 → 182.7 | 146.2 → 140.4 | 433.0 → 420.5 | 0.926 | 0.941 / 0.977 / 0.918 |
| e⁻ | 117.7 → 111.3 | 146.4 → 151.2 | 128.4 → 116.2 | 392.6 → 378.8 | 0.782 | 0.792 / 0.835 / 0.759 |

Group multiplier errors: ±0.12–0.31 (μ), ±0.08–0.11 (e).

- The postfit electron multijet falls ×0.971 (e⁺) / ×0.965 (e⁻). The prefit
  per-group ABCD cut the totals ×0.946 / ×0.943, and the fit restored part of
  that.
- **The e⁻ multipliers stay at 0.76–0.84 in every group, not only at the
  edges.** Per group, the data still want less e⁻ multijet than the ABCD
  predicts, and the e⁺ ones 0.92–0.98.
- e-only fit: e⁻ multipliers 0.871 / 0.807 / 0.797 (abcd 0.819); totals
  e⁺ 413.7 → 398.3, e⁻ 414.6 → 393.2.
- Figures: `qcd_groups_<fit>_<bn>.png` (the multijet over the abcd prefit, per group).

### 6.3 Nuisance pulls (lab)

- **μ + e:** `qcd_rate_ele_Wm` −1.67 → **−1.37**, `qcd_rate_ele_Wp` −0.35 →
  −0.13, `qcdFF_ele_Wp` −0.23 → −0.19, `qcdFF_ele_Wm` +0.03 → +0.08,
  `eleSF` +1.17 → +1.19. Every other θ moves by ≤ 0.02.
- **e only:** `qcd_rate_ele_Wm` −1.15 → −1.00, `qcd_rate_ele_Wp` −0.47 → −0.30.
- **μ only:** every |Δθ| ≤ 0.02.

The fb pulls move more, but those differences mix in the FB-window fix:

- μ + e: `qcd_rate_ele_Wm` −2.10 → −1.05, `qcd_rate_ele_Wp` −0.63 → +0.01,
  `eleSF` 0.82 → 1.32.
- e only: `qcd_rate_ele_Wm` −1.58 → −0.86.

### 6.4 r and σ (lab, the clean comparison)

**Per-bin r (μ + e):** the mean Δr/σ_stat is +0.04 and the largest is +0.45
(r_Wm_y11, 1.320 → 1.364).
- The edge bins move up: W⁺ y0 / y1 / y10 / y11 by +0.32 / +0.40 / +0.25 / +0.26
  stat, W⁻ y10 / y11 by +0.36 / +0.45.
- The central bins move by −0.03 … −0.10 stat.
- μ only: mean +0.004, max 0.25. e only: mean +0.08, max 0.36.

| fiducial σ (nb) | abcd (10-16) | abcdy | Δ / σ_stat |
|---|---|---|---|
| μ + e, σ_W⁺ | 57.99 ± 1.01 ± 1.78 | 58.16 ± 1.01 ± 1.78 | +0.16 |
| μ + e, σ_W⁻ | 45.91 ± 0.90 ± 1.42 | 45.99 ± 0.90 ± 1.42 | +0.09 |
| μ + e, σ_W | 103.90 ± 1.35 ± 3.18 | **104.15 ± 1.36 ± 3.18** | +0.18 |
| μ + e, σ_Z | 10.163 ± 0.417 ± 0.314 | 10.110 ± 0.415 ± 0.311 | −0.13 |
| μ only, σ_W / σ_Z | 103.61 / 9.385 | 103.66 / 9.384 | +0.03 / −0.00 |
| e only, σ_W⁺ / σ_W⁻ | 60.71 / 44.36 | 61.03 / 44.92 | +0.18 / +0.36 |
| e only, σ_W / σ_Z | 105.07 / 11.883 | 105.96 / 11.817 | +0.37 / −0.08 |

**e/μ from the per-flavour fits (stat errors):**
- W⁺: 1.065 → 1.069 ± 0.039
- W⁻: 0.952 → 0.964 ± 0.041
- W: 1.014 → 1.022 ± 0.028
- Z: 1.266 → 1.259 ± 0.106
- **The (e/μ)_W⁺ / (e/μ)_W⁻ double ratio: 1.120 ± 0.063 → 1.109 ± 0.062.**

Per bin (`xsec_flavfit_ratio.csv`):
- W⁻ edges: e/μ at y0 0.706 → 0.719, at y11 0.760 → 0.794.
- W⁺ edges: y0 1.005 → 1.077, y11 0.995 → 0.958.
- Central bins: within ±0.015.

μ-vs-e χ² (stat), abcd → abcdy:
- A_ch: 13.0 → 13.4 / 12
- dσ/dη W⁻: 17.4 → 15.2 / 12
- dσ/dη W⁺: 9.4 → 11.9 / 12
- inclusive W⁺, W⁻: 4.30 → 4.00 / 2

### 6.5 A_ch (lab fit) and R_FB (fb fit)

All fiducial (r × σ_gen). Overlays: `chargeAsym_<fit>.png`, `RFB_{sum,Wp,Wm}_<fit>.png`.

**A_ch:**
- μ + e: max |Δ|/σ_stat 0.28. The two most backward bins (η_CM −2.55, −2.15)
  rise by +0.012; the two most forward fall by −0.005 / −0.008.
- e only: max 0.48 (η_CM +1.85: 0.311 → 0.272; +1.45: 0.207 → 0.178;
  −2.15: 0.235 → 0.265).
- μ only: max 0.08.
- Σ(Δ/σ_stat)² over 12 bins: 0.20 / 0.70 / 0.01 (μ + e / e / μ).

**R_FB (mixes in the FB-window fix):**
- μ + e: max 1.04 σ_stat. R_FB_sum 1.049 → 1.122 at |η_CM| 1.54 and 0.859 →
  0.920 at 1.88; W⁺ max 0.66, W⁻ 0.82.
- e only: 1.34 (R_FB_sum 0.873 → 1.012 at 1.54, 0.798 → 0.931 at 1.88).
- μ only: 0.38.

**fb-window σ (mixes in the fix):**
- σ_W (μ + e) 89.49 → 91.21 nb (+1.37 stat)
- e only 88.03 → 93.23 (+2.27 stat)
- μ only 89.92 → 90.52
- fb e/μ (W) 0.979 → 1.030

The direction matches the fix alone: the fb abcd templates carried ×1.10 (μ) /
×1.29–1.34 (e) too much multijet, and the prefit fb e/μ (W⁻) moved
0.828 → 0.928 from the fix. Most of this shift is therefore expected to be the fix.

### 6.6 Postfit χ² of the rapidity-inclusive stacks (lab, `postfit_incl`)

Baker-Cousins χ² / ndf, abcd → abcdy:

| | μ + e fit, W | single-flavour fit, W |
|---|---|---|
| e | 70.2 / 14 → 68.8 / 10 | 69.4 / 15 → 67.7 / 11 |
| μ | 54.7 / 14 → 54.7 / 10 | 54.8 / 17 → 54.7 / 13 |

The extra group freedom barely lowers the raw χ²; ndf falls by the added
multijet parameters. The remaining misfit is in the lepton-pT shape, not in the
rapidity split.

### 6.7 Conclusion

- **Lab:** the per-group ABCD moves the result by ≤ 0.45 σ_stat in any r bin.
  The inclusive σ moves by ≤ 0.18 σ_stat in the μ + e fit (σ_W 103.90 → 104.15
  nb) and by 0.37 σ_stat in the e-only fit. A_ch moves by ≤ 0.48 σ_stat.
- **It does not resolve the electron W⁻ charge split, as expected from §4:**
  - the double ratio stays at 1.11 ± 0.06;
  - the W⁻ edge e/μ goes 0.71 / 0.76 → 0.72 / 0.79;
  - the e⁻ multijet multipliers are 0.76–0.84 in every group;
  - `qcd_rate_ele_Wm` relaxes only from −1.67 to −1.37.
- **fb:** the shifts are the FB-window fix and the per-group ABCD together
  (user decision: no separate abcd refit for this comparison). The combined
  effect:
  - σ_W (fb window, μ + e) 89.49 → 91.21 nb;
  - fb e/μ (W) 0.979 → 1.030;
  - R_FB moves by up to 1.04 σ_stat (μ + e) and 1.34 σ_stat (e only).

  The direction is the one the fix alone predicts (§6.5).

**How the abcdy tree reaches the observables (added with this record):**
- `plotting/disc_variants.h` gains the downstream tag `leppt_mt40_abcdy`. Its
  out-tree suffix is `_leppt_mt40_abcdy`; `Base()`, `InputSuffix()` and
  `GenFid()` map it back to leppt_mt40's Combine inputs and gen fiducial set.
- `plotting/postfit_incl.C` and `xsec_fiducial_diag` read the inputs through
  `InputSuffix()`.
- `analysis/run_observables.sh leppt_mt40_abcdy` is opt-in and not part of
  `all`. Every output carries the full tag, so the nominal outputs are never
  touched.
- Check: the nominal chain re-run with these changes reproduced every nominal
  output.
  - PNG and CSV files byte-identical.
  - PDFs differ only in their timestamps.
  - Fiducial-yield histograms 52 / 52 IDENTICAL per fit (`compare_hists.C`).
  - All A_ch / R_FB graphs identical point by point.
  - The log is identical except for one ACLiC rebuild line.

## 7. Backups (tag `pre_abcdy_mode`, the state before this change = the 2026-10-19 fbqcd inputs)

- `correction/rootfile/pre_abcdy_mode/` (`qcd_abcd_{mu,ele}.root`)
- `correction/logs/pre_abcdy_mode/` (`qcd_abcd_{mu,ele,build}.log`)
- `correction/plots/pre_abcdy_mode/` (`qcd_abcd_{mu,ele}/`)
- `plotting/plots/pre_abcdy_mode/` (the six Combine inputs + sidecars, and the
  `mt`, `met`, `leppt_mt40`, `dilepton`, `prefit_r` plots; `Elec/` the same)
- `plotting/logs/pre_abcdy_mode/` (`mtandmet`, `dileptonpeak` × 2,
  `prefit_r_leppt_mt40`)
- The fork's work dirs were not touched: every dry-run went to a scratch dir.
