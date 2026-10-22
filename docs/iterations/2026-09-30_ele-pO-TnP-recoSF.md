# Iteration 2026-09-30 — electron Reco SF from our pO tag-and-probe (one inclusive number)

> **Paths (2026-10-20):** every backup below now lives in `backups/pre_<tag>/<area>/`
> (`skim/rootfile_pre_<tag>/` → `backups/pre_<tag>/skim/rootfile/`, `<area>/pre_<tag>/` →
> `backups/pre_<tag>/<area>/`, fork `test/pO_fit_out<sfx>_pre_<tag>/` →
> `test/backups/pre_<tag>/pO_fit_out<sfx>/`); the Reco tag-and-probe outputs
> (`correction/{rootfile,plots}/reco_sf_ele/`) moved to `correction/electron_sf_check/reco_sf/{rootfile,plots}/`.
> Full list: `backups/MOVES_2026-10-20.txt`.

**Status:** inputs regenerated and validated locally on 2026-09-30; **fitted on
lxplus on 2026-09-30** (grand + μ-only + e-only, lab + fb, Asimov, impacts; fits
14:31–15:02, impacts 15:48–15:55), downloaded and the observables chain run on
2026-10-01. Predecessor: `2026-09-29_ele-ppEGM-recoSF.md` — everything not
listed under "What changed" is identical to it. **Superseded on 2026-10-01** by
`2026-10-01_ele-ppEGM-recoSF-again.md` (user decision: back to the pp EGM Reco SF);
the full outputs of this iteration are preserved in the `*pre_recoegm*` backups
listed at the end. Headline numbers re-checked against those files on 2026-10-01.

## 1. What changed (user decision, 2026-09-30)

The electron **Reco** scale factor is no longer the pp EGM 2025Prompt table
(`Reco20to75` / `RecoAbove75`) but ONE inclusive number from our own pO
Z→ee tag-and-probe, `correction/reco_sf_ele.C` (probe = a photon supercluster
with no ID, pass = any forest electron on it; nominal probe variant, inclusive
fit; `correction/rootfile/reco_sf_ele/reco_sf_ele_result.root`, `recosf_noid_incl0`):

| | value |
|---|---|
| SF_reco (applied to EVERY electron: W lepton, both Z legs, anti-iso sideband) | **0.988175 +0.01813 −0.02258 (stat, MINOS) ± 0.00164 (syst)** |
| applied ±1σ (stat ⊕ syst per side) | 1.00638 / 0.96554 |
| ε_data / ε_MC | 0.9628 (+0.0177 −0.0220) / 0.9743 |
| data pairs (60–120 GeV) | 547 pass / 53 fail (≈32 fitted background in fail) |
| fit status / covQual | 0 / 3 |
| cut-and-count 80–100 GeV | 0.9898 |
| pp EGM Reco SF on the same probes (the old applied factor) | 0.9293 |

Switch: `skim/electron_sf.h` `kEleRecoSource = kRecoTnP` (was the bool
`kApplyEleRecoSF = true` ≡ `kRecoEgm`). The TnP study probes E_T^SC > 20 GeV;
Z legs below 20 GeV get the same number (108618 of 3970704 Reco lookups in the
DY Z→ee job).

Everything else is unchanged: selection, ID+iso SF (EGM `wp90iso`), trigger SF
(1.0010 +0.0035 −0.0040), muon side, QCD method, fit model. **The datacards,
t2w maps and κ sidecars rebuilt by `run_pO_fits.sh both all --dry-run` are
identical to the as-fitted 09-29 ones up to the path prefix**; all six input
sidecars are unchanged.

## 2. Electron scale factors (skim `[SF]` blocks)

| | ⟨RECO⟩ | ⟨ID+iso⟩ | ⟨trig⟩ | **⟨SF⟩ new** | ⟨SF⟩ 09-29 | ratio |
|---|---|---|---|---|---|---|
| W⁺ (nominal W selection) | 0.98817 | 0.95385 | 1.00096 | **0.94348** | 0.88653 | 1.0642 |
| W⁻ | 0.98817 | 0.95463 | 1.00096 | **0.94424** | 0.88763 | 1.0638 |
| DY in the W selection | 0.98817 | 0.95416 | 1.00096 | 0.94378 | 0.88934 | 1.0612 |
| Z→ee pairs (DY) | 0.97649 | 0.91089 | 1.00013 | **0.88959** | 0.78785 | 1.1291 |

`eleSF` nuisance (per-bin quadrature of the three sources), on the summed
weights: W⁺ reco +1.843/−2.291%, ID+iso ±0.963%, trigger +0.352/−0.400% →
**+2.109/−2.517%** (09-29: +1.231/−1.246%); Z→ee pairs reco +3.719/−4.529%,
ID+iso +2.382/−2.349% → **+4.417/−5.102%** (09-29: +2.673/−2.640%). The Reco
term is now the largest electron-SF term: statistics-dominated (53 failing
probes), coherent over all electrons.

## 3. Electron QCD (`correction/logs/qcd_abcd_ele.log`; muon output bit-identical)

| | e⁺ new | e⁺ 09-29 | e⁻ new | e⁻ 09-29 |
|---|---|---|---|---|
| FF slope s [68%] | 0.2195 [0.1313, 0.3059] | 0.2299 | 0.2416 [0.1588, 0.3230] | 0.2506 |
| SR template A0ff | **468.0** | 476.6 | **496.9** | 504.7 |
| flat A0 = B0·C40/D0 | 418.5 | 423.3 | 440.3 | 444.7 |
| slope ±1σ | +6.3 / −5.2% | +6.3 / −5.2% | +5.9 / −4.9% | +6.0 / −4.9% |
| T_mT | 0.4866 | 0.4917 | 0.4916 | 0.4961 |
| residual κ_e (window) | 4.7% | 4.5% | 4.3% | 4.1% |

→ residual κ_e 1.047, still **1.05** in the cards (`QCD_ABCD_LNN_ELE`). κ_μ
stays 1.09 (the muon log says 1.10 — open since the embedding).

## 4. Prefit (Combine inputs, `leppt_mt40`)

| | 09-29 (pp Reco) | 2026-09-30 (pO TnP Reco) |
|---|---|---|
| Z→ee data / MC | 225 / 171.87 = 1.3091 ± 0.0873 | 225 / 194.07 = **1.1594 ± 0.0773** |
| Z→μμ data / MC | 1.0715 ± 0.0562 | 1.0715 (unchanged) |
| Z e/μ | 1.222 | **1.082** |
| W→eν r-like e⁺ / e⁻ / e | 1.2379 / 1.0806 / 1.1679 | **1.1670 / 1.0200 / 1.1015** |
| W→μν r-like μ⁺ / μ⁻ / μ | 1.1955 / 1.2242 / 1.2083 | unchanged |
| W e/μ (r-like) | 0.967 | 0.912 |
| electron W signal in the fit window (W⁺ / W⁻) | 1409.3 / 1129.6 | 1500.2 / 1202.1 |
| electron SR QCD (qcd_abcd, W⁺ / W⁻) | 476.6 / 504.7 | 468.0 / 496.9 |

r-like = (data − z − ztau − qcd_abcd)/(signal + wtau) over the fit window
24–100 GeV. The CR prefit QCD counts moved only through the EWK subtraction
(e.g. CRB e⁺ 595.2 → 588.9, e⁻ 600.5 → 595.0).

**Naive expectation for the fit** (the electron r's scale by the inverse of the
MC normalization change, ⟨RECO⟩ ratio 0.940 per W electron, 0.886 per Z→ee
pair): σ_W e/μ 1.014 → ~0.95, σ_Z e/μ 1.267 → ~1.12. What a per-electron
factor f cannot move is (Z e/μ)/(W e/μ)² — f enters the Z twice and the W
once: 1.23 in the 09-29 fit (1.267/1.014²), and prefit 1.307 → 1.301 here
(1.222/0.967² vs 1.082/0.912²). So no single per-electron correction can bring
both e/μ ratios to 1; the Reco switch only redistributes the tension from the
Z to the W.

## 5. Validation done locally (2026-09-30)

- `correction/run_reco_sf.sh all` with the modified macro: `reco_sf_ele.csv`
  byte-identical to the 09-29 one, all 70 skim histograms identical.
- Unit test of `electron_sf.h`: the Reco factor is the one number at every
  (η, pT) incl. < 20 and > 75 GeV; wp90iso still looked up at the RAW SC η
  (asymmetric); both Z legs varied coherently; `Load(log, false)` + `Reco()`
  throws.
- Electron re-skim 14/14 OK; `skim/compare_reskim.sh rootfile_pre_recotnp`
  (log `skim/logs/compare_recotnp.log`): both data files bit-identical (216 /
  14 histograms), every non-empty MC histogram changed, MISSING 0 / NEW 0 after
  `run_lhe_updown.sh Wel Zee` (max|member0 − nominal| = 0 in all 12 files).
- `run_qcd_abcd.sh both`: muon rootfile (12/12) and log identical.
- `run_combine_inputs.sh both`: muon inputs bit-identical (1938 + 2059 + 37
  histograms); electron: only `data_obs` unchanged; sidecars unchanged.
- `run_syst_shapes.sh all`: LHAPDF closure 3e-7, member 0 ≡ `signal` exactly.
- Fork `run_pO_fits.sh both all --dry-run`: cards/maps/sidecars identical to
  the as-fitted 09-29 ones (mod paths).
- Presentation only (user request): `correction/plots/reco_sf_ele/reco_sf_vs_{eta,et}`
  now draw the asymmetric MINOS intervals (e.g. endcap +0.037/−0.050, previously
  shown as ±0.043) via the new opt-in asymmetric errors of
  `plotting_helper.C` `OverlaySeries`; all other helper-based plots regenerated
  byte-identical; fork copy re-synced (uploaded with the next `sync_lxplus.sh upload`).

## 6. Fit results

Fitted with `./run_pO_fits.sh both all --asimov` + `./run_pO_impacts.sh --fit all`
(extractor call: κ μ 1.09, e 1.05, lumi 1.03, abcd, the 7 shape nuisances).
Confirmed to be the new inputs: the fit's prefit Z→ee signal is 193.84 (09-29:
171.67), the μ side unchanged at 339.17. Every number below is read from the
extraction logs, `plotting/plots/{comb,flavfit}/xsec/leppt_mt40/*.csv`, the
impacts JSONs and `analysis/logs/run_observables_leppt_mt40.log` (2026-10-01).

### Fit quality
All six fits status 0 / covQual 3; Asimov closure PASS (max |POI − 1| = 0.0000);
`--statonly` vs conditioned covariance ≤ 0.0001 (grand), 0.0002–0.0003 (μ),
0.0005–0.0006 (e). The **μ-only fit is identical to 09-29 to every printed
digit** (muon inputs bit-identical) — the control.

### Cross sections, lab binning (nb, stat ± syst)

| fit | σ(W⁺) | σ(W⁻) | σ(W) | σ(Z) | ρ |
|---|---|---|---|---|---|
| grand | 59.84 ± 1.03 ± 1.85 | 47.16 ± 0.91 ± 1.47 | **107.00 ± 1.38 ± 3.31** | **9.937 ± 0.407 ± 0.322** | +0.532 |
| μ-only | 59.66 ± 1.28 ± 1.82 | 48.62 ± 1.16 ± 1.49 | 108.28 ± 1.73 ± 3.30 | 9.413 ± 0.492 ± 0.284 | +0.389 |
| e-only | 60.04 ± 1.77 ± 2.58 | 43.74 ± 1.55 ± 1.91 | 103.78 ± 2.35 ± 4.43 | 10.689 ± 0.708 ± 0.709 | +0.584 |
| *09-29 grand* | *60.72* | *47.89* | *108.60* | *10.199* | |
| *09-29 e-only* | *63.54* | *46.27* | *109.81* | *11.923* | |

Grand effective r: W⁺ 1.2015, W⁻ 1.1900, W 1.1964 (09-29: 1.214). fb binning:
σ_W grand 92.42 ± 1.29 ± 2.90, μ 94.14, e 87.07 nb (theory 76.69); σ_Z 10.015 /
9.443 / 10.681 nb. The e-only systematic grew (σ_W 3.64 → 4.43, σ_Z 0.467 →
0.709 nb): the TnP Reco error is ~3× the pp table's.

### e/μ (statistical errors; lumi cancels)

| | 2026-09-30 | 09-29 (pp Reco) |
|---|---|---|
| σ(W⁺) | 1.0065 ± 0.0366 (total ± 0.064) | 1.0652 ± 0.0388 |
| σ(W⁻) | 0.8996 ± 0.0384 (total ± 0.062) | 0.9517 ± 0.0407 |
| σ(W) lab | **0.9585 ± 0.0266** (total ± 0.057: −1.6σ stat, −0.7σ total) | 1.0142 ± 0.0281 |
| σ(W) fb | 0.9248 ± 0.0288 | 0.9789 ± 0.0304 |
| σ(Z) lab | **1.1356 ± 0.0958** (+1.4σ stat; ≈ +1.2σ with the ≈5.9% non-lumi e-SF syst — estimate) | 1.2667 ± 0.1068 (+2.5σ) |
| σ(Z) fb | 1.1311 ± 0.0955 | 1.2610 ± 0.1064 |
| (Z e/μ)/(W e/μ)² | 1.236 | 1.231 |
| W⁺ / W⁻ double ratio | 1.12 ± 0.06 | 1.12 ± 0.06 |

Close to the pure per-electron rescaling: the W e/μ moved by 0.945 and the Z by
0.897 (the ⟨RECO⟩ ratios are 0.940 per W electron and 0.886 per Z pair; the rest
of the fit absorbs a little), and the invariant (Z e/μ)/(W e/μ)² stayed at 1.23. The Z→ee tension
is now ~1.4σ (stat) and the W sits ~1.6σ (stat) below 1 — no longer a single
outlier, but the Z-vs-W electron mismatch is unchanged, as it must be.
Lumi-free σ_W/σ_Z: μ 11.50 ± 0.65, e 9.71 ± 0.74, grand 10.77 ± 0.48 (theory 10.20).

Per lab bin y0…y11 (stat ±0.10–0.17) — the charge split, unchanged in shape:
- W⁺: 0.935 1.112 1.002 1.106 0.921 1.066 1.053 1.081 1.021 1.043 0.858 0.923
- W⁻: **0.653 0.814 0.700 0.874** 0.887 1.138 0.961 1.010 1.013 1.082 0.939 **0.702**

### Nuisance parameters (post-fit value ± constraint, lab | fb)

Grand fit:
- r_Z 1.1332 ± 0.0591 | 1.1421 ± 0.0597 (09-29: 1.1630 | 1.1696)
- **eleSF −0.175 ± 0.680 | −0.501 ± 0.668** (09-29: **+1.170** | +0.820) — the
  tension it absorbed is gone; its width 0.68 = the grand fit measuring the
  electron normalization against the muons through the shared r's
- muSF −0.045 | −0.076; lumi +0.017 | +0.022; alphaS −0.062 | −0.116
- nPDF −0.128 ± 0.938 | −0.065 ± **1.498**; qcdScale −0.211 | +0.028
- qcd_rate: μ⁺ −0.201 | −0.157, μ⁻ +0.103 | +0.035, e⁺ −0.565 | −0.807, **e⁻ −1.822 | −2.225**
- qcdFF_ele: W⁺ −0.141 ± 0.489 | −0.758 ± 0.615, W⁻ +0.087 ± 0.493 | −0.084 ± 0.497
- SR QCD multipliers (lab): μ⁺ 0.919, μ⁻ 0.975, e⁺ 0.901, e⁻ 0.767

e-only: r_Z **1.2190 ± 0.1142** | 1.2181 ± 0.1123 (09-29: 1.3597); eleSF
−0.159 ± 1.101 | −0.236 ± 1.029; qcd_rate e⁺ −0.475 | −0.356, e⁻ −1.180 |
−1.609; nPDF −0.065 ± 0.957 | −0.173 ± **1.498**; qcdScale −0.311 | −0.198;
qcdFF W⁺ −0.145 ± 0.490 | −0.577 ± 0.580, W⁻ +0.296 ± 0.491 | +0.121 ± 0.463;
multipliers 0.902 / 0.830 (lab), 0.919 / 0.783 (fb).

**nPDF in the fb binning has a post-fit width of 1.498 in both the grand and the
e-only fit** (09-29: 1.007 / 1.027; lab 0.94 / 0.96; μ-only fb 0.970) — a
non-parabolic direction of the electron fb likelihood (the Hesse curvature below
the constraint's), not a failure: closure, the statonly cross-check and the
impacts (nPDF 6.9% / 4.3% of σ_r) are normal. The count of one-sided template
bins is unchanged (24), so it is the likelihood, not the templates, that changed.

### Goodness of fit (inclusive postfit stacks, lab, BC χ²/ndf)

| fit | μ W⁺ | μ W⁻ | μ W | e W⁺ | e W⁻ | e W |
|---|---|---|---|---|---|---|
| grand | 62.3/28 | 44.6/28 | 55.7/15 | 53.5/28 | 67.1/28 | 71.7/15 |
| per flavour | 62.3/31 | 43.2/31 | 55.2/18 | 53.1/29 | 67.4/29 | 71.4/16 |

Unchanged in substance (a flat factor does not reshape); grand e W⁻
data/postfit 0.955 (09-29: 0.969) — the shared r's now over-predict the e⁻.
μ-vs-e shape (stat χ²): A_ch 13.1/12, R_FB 6.2/6 (W⁺ 8.0/6, W⁻ 1.2/6),
dσ/dη W⁺ **5.9/12** (09-29 9.5), W⁻ **23.2/12, p = 0.026** (09-29 17.5), W 16.1/12
— the e/μ normalization offset now sits in W⁻.

Primary grand-fit observables vs 09-29: A_ch moved ≤ 0.001, R_FB ≤ 0.010
(≤ 0.14 stat σ); e-only ≤ 0.002 / ≤ 0.026 (≤ 0.22σ) — a flat per-electron
factor cancels in the ratios, as it must.

### Impacts (mean over the 24 W r's of impact/σ_r; on r_Z), lab
- grand: lumi 44.1% (57.6%), **eleSF 9.1% (22.1%)** (09-29: 5.4% / 13.7%), nPDF 6.4%, qcdScale 4.8%, muSF 4.5% (9.5%), qcd_sB/qcd_rate ≤ 2.7% each
- μ-only: unchanged (lumi 36.4%, nPDF 5.6%, muSF 4.7%, ...)
- e-only: lumi 27.4% (33.1%), **eleSF 24.8% (61.0%)** (09-29: 12.0% / 33.6%), qcd_sB_ele 9.8/9.5% (29.4/28.2% on r_Z), qcd_rate_ele 3.9/4.2%, nPDF 3.4%

### Profiled (σ_W, σ_Z) contours
All drawn; scan minimum vs covariance best fit ≤ 0.014 nb (gate passes).
Ellipse vs profiled 68% contour: grand −3.9/+3.2% (lab), −4.1/+3.3% (fb); μ
−4.0/+3.6% (both); e −13.1/+3.7% (lab), −10.5/+4.2% (fb) (09-29: e −5.0/+4.7%,
−5.1/+4.7%) — the larger, asymmetric e-SF error makes the e profile less
Gaussian. **New on this fit: failed grid points** — combine wrote no row for
16 nodes of the grand lab scan and 22 of the e-only lab scan, each set on a
straight line (σ_W + 1.8–1.9 × σ_Z ≈ const; the 09-29 scans had none). The
grand line lies outside the 99% region (ΔNLL 5.3–8.8); the e-only line is
NOT beyond 3σ, as first written here: it passes through ΔNLL ≈ 2 (inside the
95% region), and two missing nodes, (σ_W, σ_Z) = (111.40, 13.33) and
(115.41, 11.25), sit on the 95% contour. The drawn contour bridges each hole
from its neighbours; the surface is smooth across every hole (ΔNLL second
differences ≈ 0.05, as everywhere else), so the effect is well below a grid
step. Cause, checked 2026-10-01: a profile fit that did not converge — combine
runs ONE Migrad (strategy 0, tolerance 0.1) per node, no fallback configured,
and drops a failed node without a message. NOT a parameter at a range limit,
which was the first guess here (the bounded a_i and qcd_s* stay ≥ 0.67 from
their limits, Gaussian-propagated to the holes), and not a feature of the
likelihood (none in the stat-only, μ-only or fb scans). Which configuration
trips the minimizer is not recorded (no fit status or profiled nuisances
saved); the test is a rerun with `--cminFallbackAlgo Minuit2,Migrad,1:0.1 -v 1`.
But one
missing node was a CORNER of the e-only window, which cut the convex hull, and
`TGraph2D` fills bins outside the hull with 0 → a fake island (a 2nd contour
piece, "+455%" in the Gaussianity line, outside the plotted frame). Fixed in
`plotting/xsec_contour.C` (`SetMarginBinsContent(1e3)`); every contour PNG/CSV
regenerated byte-identical, only the log lines changed.

### Summary
The pO Reco SF removed the single Z→ee outlier: σ_Z e/μ 1.27 → 1.14 (+2.5σ →
+1.4σ stat), the grand fit no longer pulls `eleSF` (+1.17σ → −0.18σ). The price
is the expected one: σ_W e/μ 1.01 → 0.96 (−1.6σ stat, −0.7σ total), because
(Z e/μ)/(W e/μ)² ≈ 1.23 is invariant under any per-electron factor. Still open,
untouched by this iteration: the electron charge split (W⁻ edge bins y0/y2/y11
at 0.65–0.70 of the muon, e⁻ QCD pull −1.8σ) — the next step.

## 7. Backups of this iteration (taken 2026-10-01)

| path | content |
|---|---|
| `skim/rootfile_pre_recoegm/` | the 14 electron skim files (`WToElecNu_*`, `ZToEE_*`), with the LHE Up/Down twins |
| `skim/rootfile_pre_recoegm/observables/` | `{fidyields,charge_asym_fid,FBratio_fid}_{comb,simfit_mu,simfit_ele}_leppt_mt40.root` |
| `skim/logs/pre_recoegm/` | `skim_{Wel,Zee}_*.log`, `lhe_updown_{Wel,Zee}.log`, `compare_recotnp.log` (the 09-30 gate vs 09-29) |
| `skim/output/pre_recoegm/` | `WToElecNu_pO_PFMet_*.txt` cutflows |
| `correction/rootfile/pre_recoegm/` | `qcd_abcd_{ele,mu}.root` |
| `correction/logs/pre_recoegm/` | `qcd_abcd_{ele,mu}.log` |
| `plotting/plots/pre_recoegm/{Elec,mu}/` | all Combine inputs + `_systs.txt` sidecars |
| `plotting/plots/pre_recoegm/{comb,flavfit}/` | every observables plot and CSV of this fit |
| `plotting/logs/pre_recoegm/` | `mtandmet_*`, `dileptonpeak_*`, `syst_shapes_*` logs |
| `analysis/logs/pre_recoegm/` | `run_observables_leppt_mt40.log` |
| fork `test/pO_fit_out_leppt_mt40_pre_recoegm/` | the downloaded fit trees `simfit/`, `simfit_mu/`, `simfit_ele/` (summaries, cards as fitted, fitDiagnostics, contours, impacts, logs) |

Not re-run, hence not backed up: the Reco tag-and-probe outputs
(`correction/rootfile/reco_sf_ele/`, `correction/plots/reco_sf_ele/`), unchanged.
