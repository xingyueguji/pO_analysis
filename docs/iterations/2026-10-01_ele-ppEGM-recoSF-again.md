# Iteration 2026-10-01 — electron Reco SF back to the pp EGM tables (the EGamma official SF)

> **Paths (2026-10-20):** every backup below now lives in `backups/pre_<tag>/<area>/`
> (`skim/rootfile_pre_<tag>/` → `backups/pre_<tag>/skim/rootfile/`, `<area>/pre_<tag>/` →
> `backups/pre_<tag>/<area>/`, fork `test/pO_fit_out<sfx>_pre_<tag>/` →
> `test/backups/pre_<tag>/pO_fit_out<sfx>/`). Full list: `backups/MOVES_2026-10-20.txt`.

**Status:** inputs regenerated and validated locally on 2026-10-01; **fitted on
lxplus on 2026-10-05** (grand + μ-only + e-only, lab + fb, Asimov, impacts; fits
09:52–11:38, impacts 17:22–17:41), downloaded and the observables chain run the
same day. **The fit reproduces the 09-29 fit byte for byte** (section 7).
Predecessor: `2026-09-30_ele-pO-TnP-recoSF.md`.

**The inputs of this iteration are bit-identical to those of
`2026-09-29_ele-ppEGM-recoSF.md`** — every electron skim histogram, every
Combine input and sidecar, every datacard (section 5). Its sections 1–2
(definition, scale factors, QCD, prefit) therefore apply unchanged, and they
were re-read from this run's files below. Its fit trees are preserved in the
fork's `test/pO_fit_out_leppt_mt40_pre_recotnp/`.

## 1. What changed (user decision, 2026-10-01)

"Use the EGamma official SF for reco, for electron again": `skim/electron_sf.h`
`kEleRecoSource = kRecoTnP → kRecoEgm`.

| | 2026-09-30 (predecessor) | 2026-10-01 (this iteration) |
|---|---|---|
| applied electron Reco SF | pO Z→ee tag-and-probe, ONE inclusive number for every electron: 0.988175 +0.01813 −0.02258 (stat) ± 0.00164 (syst) | pp EGM 2025Prompt `Reco20to75` (pT < 75) / `RecoAbove75`, per electron in (SC η, pT): 0.9138–0.9532 at 20–75 GeV, 0.9420–0.9681 at 75–500 GeV; η-symmetric |
| Z legs below 20 GeV | the same inclusive number | the 20–45 GeV cell (no `RecoBelow20`): 108618 of 3970704 Reco lookups in the DY Z→ee job |
| ±1σ of the Reco factor | stat ⊕ syst per side: 1.00638 / 0.96554 | the EGM `sfup/sfdown` per cell (errors 0.4–1.1%), coherent over all electrons |

Everything else unchanged: selection (WP90 ID+iso, crack veto), ID+iso SF (EGM
`wp90iso`), trigger SF (1.0010 +0.0035 −0.0040), the muon side, the QCD method
(in-fit ABCD, fake-factor electron template), the fit model. The pO tag-and-probe
(`correction/reco_sf_ele.C`) stays as the pO cross-check; its result file is no
longer read by the skim (only under `kRecoTnP`). The comments that described the
TnP number as the applied factor were updated (`skim/electron_sf.h`, `skim/skim.C`,
`correction/reco_sf_ele.C`, `correction/run_reco_sf.sh`, README §2.6, CLAUDE.md).

## 2. Electron scale factors (skim `[SF]` blocks)

| | ⟨RECO⟩ | ⟨ID+iso⟩ | ⟨trig⟩ | **⟨SF⟩ this** | ⟨SF⟩ 09-30 |
|---|---|---|---|---|---|
| W⁺ (nominal W selection) | 0.92861 | 0.95385 | 1.00096 | **0.88653** | 0.94348 |
| W⁻ | 0.92899 | 0.95463 | 1.00096 | **0.88763** | 0.94424 |
| DY in the W selection | 0.93123 | 0.95416 | 1.00096 | 0.88934 | 0.94378 |
| Z→ee pairs (DY) | 0.86486 | 0.91089 | 1.00013 | **0.78785** | 0.88959 |

`eleSF` nuisance on the summed weights: W⁺ reco ±0.684%, ID+iso ±0.961%, trigger
+0.352/−0.400% → **+1.231/−1.246%** (09-30: +2.109/−2.517%); Z→ee pairs reco
+1.216/−1.209%, ID+iso +2.380/−2.346%, trigger +0.014/−0.018% → **+2.673/−2.640%**
(09-30: +4.417/−5.102%). The `[SF]` blocks are identical to the 09-29 logs except
for the label "APPLIED = the pp EGM tables:" (added on 09-30).

## 3. Electron QCD (`correction/logs/qcd_abcd_ele.log`, identical to the 09-29 log)

| | e⁺ | e⁻ |
|---|---|---|
| FF slope s (68% profile likelihood) | 0.2299 [0.1429, 0.3152] (F ×1.258 / 10 GeV) | 0.2506 [0.1688, 0.3312] (×1.285) |
| SR template total A0ff | **476.6** (09-30: 468.0) | **504.7** (09-30: 496.9) |
| flat A0 = B0·C40/D0 | 423.3 | 444.7 |
| slope ±1σ twins | 506.5 / 451.7 (+6.3 / −5.2%) | 534.9 / 479.9 (+6.0 / −4.9%) |
| T_mT | 0.4917 | 0.4961 |

Residual κ_e (window) 1.045 → **1.05** in the cards; the muon output is unchanged
(κ_μ 1.09 in the cards, the log's 1.10 still open since the embedding).

## 4. Prefit (Combine inputs, `leppt_mt40`)

| | data | DY / signal MC | other MC | ratio |
|---|---|---|---|---|
| Z→μμ peak | 364 | 339.17 | 0.54 | data/MC 1.0715 ± 0.0562 |
| Z→ee peak | 225 | 171.67 | 0.20 | data/MC **1.3091 ± 0.0873** (09-30: 1.1594); Z e/μ 1.222 (09-30: 1.082) |

W, inclusive over the fit window 24–100 GeV, r-like = (data − z − ztau − qcd_abcd)/(signal + wtau):

| | data | signal | z | ztau | wtau | qcd_abcd | r-like | data/total |
|---|---|---|---|---|---|---|---|---|
| μ | 4998 | 3678.9 | 177.1 | 8.7 | 93.0 | 254.7 | 1.2083 | 1.1865 |
| e⁺ | 2299 | 1409.3 | 34.9 | 3.0 | 32.3 | 476.6 | 1.2379 | 1.1753 |
| e⁻ | 1791 | 1129.6 | 33.6 | 2.9 | 27.0 | 504.7 | 1.0806 | 1.0549 |
| e | 4090 | 2538.8 | 68.5 | 5.9 | 59.3 | 981.3 | **1.1679** (09-30: 1.1015) | 1.1194 |

W e/μ (r-like) 0.967 (09-30: 0.912). `eleSF` on the electron W signal templates
+1.29/−1.31%, on the Z→ee signal +2.67/−2.64%. Electron CR prefit QCD (data):
CRB e⁺ 595.15 (694), e⁻ 600.51 (688); CRC 860.89 (872) / 896.28 (905); CRD
1210.34 (1211) / 1210.36 (1211).

## 5. Validation done locally (2026-10-01)

- Electron re-skim `run_all.sh Wel`, `Zee`: 14/14 OK (data 106 s / 91 s, MC 28–37 s).
- `skim/compare_reskim.sh rootfile_pre_recotnp` (log `skim/logs/compare_recotnp.log`):
  **all 14 files IDENTICAL to the 09-29 state** — 1848 histograms per W file
  (incl. the LHE Up/Down twins), 31 per Z MC file, 216 / 14 in the data files;
  MISSING 0, NEW 0.
- `skim/compare_reskim.sh rootfile_pre_recoegm` (log `skim/logs/compare_recoegm.log`):
  vs the 09-30 state both data files identical, every non-empty MC histogram
  changed (the 77 / 264 / 378 identical ones in W⁺ / W⁻τ / W⁺τ and the 29 / 1 in
  the Z W⁻τ / W⁺τ files are all empty).
- `run_lhe_updown.sh Wel`, `Zee`: max|member0 − nominal| = 0 in all 12 files;
  both logs identical to the 09-29 ones.
- `run_qcd_abcd.sh both`: `qcd_abcd_ele.root` 18/18 histograms and the log
  identical to 09-29; `qcd_abcd_mu.root` 12/12 identical to both 09-29 and 09-30.
- `run_combine_inputs.sh both`: all six inputs bit-identical to the 09-29 ones
  (electron 1938 / 2161 / 37 histograms, muon 1938 / 2059 / 37); vs 09-30 only the
  electron `data_obs` templates (51 / 57 / 1) are unchanged; all six `_systs.txt`
  sidecars unchanged; the four writer logs identical to 09-29.
- `run_syst_shapes.sh all`: LHAPDF closure 2.95e-7 (met) / 5.92e-8 (leppt_mt40),
  member 0 ≡ `signal` exactly; both logs identical to 09-29.
- Fork `run_pO_fits.sh both all --dry-run`: the 21 files per set (cards, t2w maps,
  σ-POI maps, κ sidecars of the grand, μ-only and e-only fits) identical to the
  as-fitted 09-30 AND 09-29 ones once the absolute path prefix (lxplus `/afs/…`
  vs local) is substituted — 42/42 comparisons. Sidecar: κ μ 1.09, e 1.05, lumi
  1.03, `qcdMode abcd`, `lheSysts nPDF,qcdScale,alphaS,muSF,eleSF,qcdFF_ele_Wp,qcdFF_ele_Wm`.
  **No fork change.**

## 6. Expected fit outcome

Same inputs, same cards, and the fork's fit code unchanged since before the 09-29
fit (`make_pO_simfit_cards.sh` and `extract_pO_simfit.C` last modified 09-29 16:25,
`run_pO_fits.sh` 09-22; the fit ran at 17:34; only the plotting helper changed
since). The fit should therefore reproduce the 09-29 results: σ_W grand
108.60 ± 1.41 ± 3.31 nb, σ_Z 10.199 ± 0.418 ± 0.314 nb, σ_Z e/μ 1.2667 ± 0.1068
(+2.5σ), σ_W e/μ 1.0142 ± 0.0281, `eleSF` pull +1.170 (lab) / +0.820 (fb), e-only
r_Z 1.3597 ± 0.1045. A deviation would point to the lxplus environment, not to the
inputs. The 09-29 profiled-contour scans had no failed grid nodes (the 09-30 ones
lost 16 / 22).

## 7. Fit results

Fitted with `./run_pO_fits.sh both all --asimov` + `./run_pO_impacts.sh --fit all`;
downloaded trees checked first: every file in `simfit/`, `simfit_mu/`, `simfit_ele/`
is dated 2026-10-05 (nothing left over from the 09-30 download). Observables:
`analysis/run_observables.sh leppt_mt40` (log `analysis/logs/run_observables_leppt_mt40.log`).

### Identity with the 09-29 fit (the expectation of section 6, confirmed)

Compared with the preserved 09-29 trees (`test/pO_fit_out_leppt_mt40_pre_recotnp/`)
and outputs (`*pre_recotnp*`):

| | result |
|---|---|
| summary CSVs of all three fits (`comb_*`, `simfit_{mu,ele}_*`) | byte-identical |
| extraction logs `extract_simfit{,_mu,_ele}.log` (every pull, closure, cross-check) | byte-identical |
| cards, t2w maps, σ-POI maps, κ sidecars as fitted (21 files) | byte-identical |
| `*_fitted_yields.root` (3 × 104 histograms) | identical (`compare_hists.C`) |
| impacts JSONs (6) | byte-identical |
| observables CSVs (15: `xsec_comb`, `xsec_diag[_bins]`, both contour sets, `xsec_flavfit[_ratio]`, `mu_vs_e_chi2` ×2) | byte-identical |
| `fidyields_{comb,simfit_mu,simfit_ele}_leppt_mt40.root` (3 × 52 histograms) | identical |
| `run_observables_leppt_mt40.log` | identical but for one ACLiC build line and the line added on 10-01 (below) |

So every number of the 09-29 record's section 3 holds unchanged for this
iteration. Plot images differ only where the code changed since 09-29: the μ-vs-e
χ² now drawn on `W_fiducial` / `W_dsigma_deta_*` (10-01).

### Fit quality
All six fits status 0 / covQual 3; Asimov closure PASS (max |POI − 1| = 0.0000);
`--statonly` companion vs conditioned covariance: max |rErr ratio − 1| = 0.0002
(grand), 0.0002/0.0003 (μ), 0.0005/0.0006 (e). Profiled contours: all 12 scans
complete (2501 points each, one contour piece — no failed grid nodes, unlike the
09-30 fit), scan minimum vs covariance best fit ≤ 0.0085 nb; ellipse vs profiled
68% contour d² 2.124–2.448 (grand lab), 2.072–2.517 (e-only lab).

### Headline numbers (re-read from this run's files)

| fit | σ(W⁺) | σ(W⁻) | σ(W) | σ(Z) |
|---|---|---|---|---|
| grand | 60.72 ± 1.05 ± 1.85 | 47.89 ± 0.93 ± 1.47 | **108.60 ± 1.41 ± 3.31** | **10.199 ± 0.418 ± 0.314** |
| μ-only | 59.66 ± 1.28 ± 1.82 | 48.62 ± 1.16 ± 1.49 | 108.28 ± 1.73 ± 3.30 | 9.413 |
| e-only | 63.54 ± 1.87 ± 2.13 | 46.27 ± 1.64 ± 1.59 | 109.81 ± 2.49 ± 3.64 | 11.923 |

(nb, stat ± syst; grand fb σ_W 93.70 ± 1.32 ± 2.87, σ_Z 10.257.) Grand effective
r: W⁺ 1.219, W⁻ 1.209, W 1.214.

e/μ (stat): σ(W⁺) 1.0652 ± 0.0388, σ(W⁻) 0.9517 ± 0.0407, **σ(W) 1.0142 ± 0.0281**
(fb 0.9789 ± 0.0304), **σ(Z) 1.2667 ± 0.1068** (+2.5σ; fb 1.2610 ± 0.1064). New
console line since 10-01: the inclusive μ-vs-e χ² over W⁺ and W⁻ = 4.32/2 (stat,
p = 0.116), 1.79/2 (total, p = 0.408). dσ/dη μ vs e: W⁺ 9.5/12, W⁻ 17.5/12, W 13.9/12
(stat).

Grand-fit pulls (lab | fb): r_Z 1.1630 ± 0.0596 | 1.1696 ± 0.0600; **eleSF +1.170 ±
0.835 | +0.820 ± 0.850**; qcd_rate_ele_Wm −1.678 | −2.114, e⁺ −0.359 | −0.641; μ⁺
−0.304 | −0.246, μ⁻ −0.016 | −0.062; nPDF −0.180 | −0.228; qcdScale −0.267 | −0.022;
alphaS −0.064 | −0.125; muSF −0.306 | −0.270; lumi +0.016 | +0.021; qcdFF_ele W⁺
−0.213 ± 0.485 | −0.833 ± 0.602, W⁻ +0.043 ± 0.477 | −0.127 ± 0.498. e-only r_Z
1.3597 ± 0.1045, μ-only 1.0734 ± 0.0648.

Inclusive postfit BC χ²/ndf (lab): grand μ W⁺/W⁻/W 62.7/28, 43.5/28, 55.2/15; e
55.0/28, 65.8/28, 71.1/15; μ-only 62.3/31, 43.2/31, 55.2/18; e-only 52.5/29,
66.4/29, 70.2/16.

### Reading
With the EGM Reco SF the Z→ee tension is back at its 09-29 size: σ_Z e/μ = 1.27
(+2.5σ stat), the grand fit pulls `eleSF` by +1.17σ, the e-only r_Z is 1.36 against
μ-only 1.07. The 09-30 iteration (pO TnP Reco 0.988) had brought these to 1.14
(+1.4σ), −0.18σ and 1.22, at the price of σ_W e/μ 1.01 → 0.96. Neither Reco choice
moves (Z e/μ)/(W e/μ)² ≈ 1.23, and the electron charge split (double ratio 1.12 ±
0.06) is untouched by either — still the next item.

## 8. Backups

No new backup is needed for this iteration's inputs: they are identical to the
09-29 set, preserved in the `*_pre_recotnp` copies (listed in the 09-29 record).
The predecessor (09-30, pO TnP Reco SF) was backed up as `*_pre_recoegm` before
this iteration started — listed in its record.
