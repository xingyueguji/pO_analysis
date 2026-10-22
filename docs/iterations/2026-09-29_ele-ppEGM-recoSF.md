# Iteration 2026-09-29 — electron WP90 ID+iso, pp EGM Reco SF, fake-factor electron QCD

> **Paths (2026-10-20):** every backup below now lives in `backups/pre_<tag>/<area>/`
> (`skim/rootfile_pre_<tag>/` → `backups/pre_<tag>/skim/rootfile/`, `<area>/pre_<tag>/` →
> `backups/pre_<tag>/<area>/`, fork `test/pO_fit_out<sfx>_pre_<tag>/` →
> `test/backups/pre_<tag>/pO_fit_out<sfx>/`); the SF-study outputs moved to
> `correction/electron_sf_check/`. Full list: `backups/MOVES_2026-10-20.txt`.

**Status:** fitted on lxplus on 2026-09-29 (grand + μ-only + e-only fits, lab + fb,
Asimov closure, impacts). **Superseded on 2026-09-30** by the iteration that
replaces the pp EGM electron Reco SF with our pO tag-and-probe value
(`2026-09-30_ele-pO-TnP-recoSF.md`).

Every number below was **read back from the fit/extraction outputs and job logs
on 2026-09-30**, not copied from notes. The full outputs are preserved in the
`*pre_recotnp*` backups listed at the end.

---

## 1. What defined this iteration

| item | setting |
|---|---|
| MC / data | Aug-20 embedded MC (POWHEG + Angantyr), July-29 data, L = 46.5 nb⁻¹ |
| electron selection | `eleMVAIdWP90` (ECAL crack vetoed on η_SC) + iso-pass = `eleMVAIsoWP90 && Δβ relIso < 0.3`; pT > 25, \|η\| < 2.4 (on `eleEta`), trigger match ΔR < 0.4 |
| muon selection | TightID, Δβ relIso < 0.15 (unchanged since 2026-09-14/21) |
| W discriminant | lepton pT at m_T > 40 (`leppt_mt40`), fit window 24–100 GeV |
| QCD | in-fit ABCD (`QCD_MODE=abcd`); electron SR template fake-factor weighted |
| fit | `./run_pO_fits.sh both all --asimov` + `./run_pO_impacts.sh --fit all` |
| κ in the cards | μ 1.09, e 1.05, lumi 1.03 (from the extractor call in `extract_simfit.log`) |
| shape nuisances | nPDF, qcdScale, alphaS, muSF, eleSF, qcdFF_ele_Wp, qcdFF_ele_Wm |

### Electron scale factors (`skim/electron_sf.h`, `kApplyEleRecoSF = true`)

| factor | source | ⟨factor⟩ W⁺ / W⁻ | per Z→ee pair |
|---|---|---|---|
| RECO | pp EGM 2025Prompt `Reco20to75` / `RecoAbove75` (0.9138–0.9532 at 20–75 GeV) | 0.92861 / 0.92899 | 0.86486 |
| ID+iso | EGM `wp90iso` | 0.95385 / 0.95463 | 0.91089 |
| trigger | our MB study, inclusive 1.0010 (+0.0035 −0.0040), 2084/2139 in data | 1.00096 | 1.00013 |
| **total ⟨SF⟩** | | **0.88653 / 0.88763** (DY-in-W 0.88934) | **0.78785** |

Z legs below 20 GeV clamped to the 20–45 GeV Reco cell: 108618 of 3970704 Reco
lookups in the DY Z→ee job.

`eleSF` nuisance on the summed weights: W⁺ reco ±0.684%, ID+iso ±0.961%, trigger
+0.352/−0.400% → **+1.231/−1.246%**; Z→ee pairs reco +1.216/−1.209%, ID+iso
+2.380/−2.346%, trigger +0.014/−0.018% → **+2.673/−2.640%**.

Muon SFs (unchanged): ⟨SF⟩ W⁺ 0.98792, W⁻ 0.98756, Z→μμ pairs 0.98247.

### Electron QCD (`correction/qcd_abcd.C`, log `correction/logs/qcd_abcd_ele.log`)

| | e⁺ | e⁻ |
|---|---|---|
| FF slope s (68% profile likelihood) | 0.2299 [0.1429, 0.3152] (F ×1.258 / 10 GeV) | 0.2506 [0.1688, 0.3312] (×1.285) |
| SR template total A0ff (FF) | 476.6 | 504.7 |
| flat A0 = B0·C40/D0 | 423.3 | 444.7 |
| slope ±1σ twins | 506.5 / 451.7 (+6.3 / −5.2%) | 534.9 / 479.9 (+6.0 / −4.9%) |
| T_mT | 0.4917 | 0.4961 |
| window residual → κ_e | 4.5% | 4.1% → **1.045 → card 1.05** |

## 2. Prefit (Combine inputs, `leppt_mt40`)

| | data | DY / signal MC | other MC | ratio |
|---|---|---|---|---|
| Z→μμ peak | 364 | 339.17 | 0.54 | data/MC **1.0715 ± 0.0562** |
| Z→ee peak | 225 | 171.67 | 0.20 | data/MC **1.3091 ± 0.0873** (e/μ 1.222) |

W, inclusive over the fit window, r-like = (data − z − ztau − qcd_abcd)/(signal + wtau):

| | data | signal | z | ztau | wtau | qcd_abcd | r-like | data/total |
|---|---|---|---|---|---|---|---|---|
| μ⁺ | 2730 | 2043.9 | 92.7 | 4.4 | 50.1 | 129.5 | 1.1955 | 1.1764 |
| μ⁻ | 2268 | 1635.0 | 84.4 | 4.3 | 42.9 | 125.2 | 1.2242 | 1.1988 |
| μ | 4998 | 3678.9 | 177.1 | 8.7 | 93.0 | 254.7 | **1.2083** | 1.1865 |
| e⁺ | 2299 | 1409.3 | 34.9 | 3.0 | 32.3 | 476.6 | 1.2379 | 1.1753 |
| e⁻ | 1791 | 1129.6 | 33.6 | 2.9 | 27.0 | 504.7 | 1.0806 | 1.0549 |
| e | 4090 | 2538.8 | 68.5 | 5.9 | 59.3 | 981.3 | **1.1679** | 1.1194 |

Selected data: W→eν candidates 6182 (m_T > 40: 4158), Z→ee 225; W→μν 6055
(5050), Z→μμ 364. `eleSF` on the W signal templates +1.29/−1.31%.

## 3. Fit results

### Fit quality
All six fits (grand / μ-only / e-only × lab / fb): status 0, covQual 3; Asimov
closure PASS (max |POI − 1| = 0.0000); `--statonly` companion vs conditioned
covariance: max |rErr ratio − 1| = 0.0002 (grand), 0.0002/0.0003 (μ),
0.0005/0.0006 (e); every profiled contour drawn (scan minimum within 0.009 nb of
the covariance best fit).

### Cross sections, lab binning (nb, stat ± syst)

| fit | σ(W⁺) | σ(W⁻) | σ(W) | σ(Z) | ρ(σ_W, σ_Z) |
|---|---|---|---|---|---|
| grand (μ+e) | 60.72 ± 1.05 ± 1.85 | 47.89 ± 0.93 ± 1.47 | **108.60 ± 1.41 ± 3.31** | **10.199 ± 0.418 ± 0.314** | +0.516 |
| μ-only | 59.66 ± 1.28 ± 1.82 | 48.62 ± 1.16 ± 1.49 | 108.28 ± 1.73 ± 3.30 | 9.413 ± 0.492 ± 0.284 | +0.389 |
| e-only | 63.54 ± 1.87 ± 2.13 | 46.27 ± 1.64 ± 1.59 | 109.81 ± 2.49 ± 3.64 | 11.923 ± 0.789 ± 0.467 | +0.389 |
| theory (r = 1) | 49.80 | 39.63 | 89.43 | 8.769 | |

Grand-fit effective r: W⁺ 1.219, W⁻ 1.209, W 1.214. fb binning: σ_W grand
93.70 ± 1.32 ± 2.87, μ 94.14, e 92.16 nb (theory 76.69); σ_Z 10.257 / 9.443 /
11.908 nb.

### e/μ (statistical errors; lumi cancels)

| | lab | fb |
|---|---|---|
| σ(W⁺) | 1.0652 ± 0.0388 | |
| σ(W⁻) | 0.9517 ± 0.0407 | |
| σ(W) | **1.0142 ± 0.0281** | 0.9789 ± 0.0304 |
| σ(Z) | **1.2667 ± 0.1068** (+2.5σ) | 1.2610 ± 0.1064 |

Per lab bin y0…y11 (stat errors ±0.10–0.17):
- W⁺: 1.004 1.154 1.063 1.172 0.973 1.132 1.117 1.143 1.082 1.107 0.892 0.994
- W⁻: **0.704 0.845 0.740 0.926** 0.933 1.207 1.019 1.067 1.074 1.150 0.975 **0.759**

### Nuisance parameters (post-fit value ± constraint, lab | fb)

Grand fit:
- r_Z 1.1630 ± 0.0596 | 1.1696 ± 0.0600
- **eleSF +1.170 ± 0.835 | +0.820 ± 0.850**; muSF −0.306 | −0.270; lumi +0.016 | +0.021
- nPDF −0.180 | −0.228; qcdScale −0.267 | −0.022; alphaS −0.064 | −0.125
- qcd_rate: μ⁺ −0.304 | −0.246, μ⁻ −0.016 | −0.062, e⁺ −0.359 | −0.641, **e⁻ −1.678 | −2.114**
- qcdFF_ele: W⁺ −0.213 ± 0.485 | −0.833 ± 0.602, W⁻ +0.043 ± 0.477 | −0.127 ± 0.498
- SR QCD multipliers (lab): μ⁺ 0.897, μ⁻ 0.950, e⁺ 0.925, e⁻ 0.781

μ-only: r_Z 1.0734 ± 0.0648 | 1.0768 ± 0.0650; every |θ| < 0.25 (qcd_rate_mu
−0.177/−0.166, nPDF −0.125, qcdScale +0.116, muSF −0.099); multipliers 0.930 / 0.934.

e-only: r_Z **1.3597 ± 0.1045** | 1.3579 ± 0.1047; qcd_rate_ele e⁺ −0.469 | −0.344,
e⁻ −1.157 | −1.582; nPDF −0.103 | −0.274; qcdScale −0.354 | −0.241; eleSF
+0.160 | +0.024; qcdFF_ele W⁺ −0.271 ± 0.519 | −0.748 ± 0.607, W⁻ +0.171 ± 0.470
| +0.018 ± 0.457; multipliers 0.886 / 0.819 (lab), 0.904 / 0.773 (fb).

### Goodness of fit (inclusive postfit stacks, lab, Baker-Cousins χ²/ndf)

| fit | μ W⁺ | μ W⁻ | μ W | e W⁺ | e W⁻ | e W |
|---|---|---|---|---|---|---|
| grand | 62.7/28 | 43.5/28 | 55.2/15 | 55.0/28 | 65.8/28 | 71.1/15 |
| per flavour | 62.3/31 | 43.2/31 | 55.2/18 | 52.5/29 | 66.4/29 | 70.2/16 |

μ-vs-e shape agreement (stat χ²): A_ch 13.1/12, R_FB 6.4/6 (W⁺ 8.2/6, W⁻
1.2/6), dσ/dη W⁺ 9.5/12, W⁻ 17.5/12, W 13.9/12.

### Impacts (mean over the 24 W r's of impact/σ_r; in brackets on r_Z), lab
- grand: lumi 44.2% (58.1%), nPDF 6.3% (4.5%), **eleSF 5.4% (13.7%)**, qcdScale 5.0% (4.0%), muSF 4.1% (8.1%), qcd_sB/qcd_rate ≤ 2.6% each (≤ 6.8% max)
- μ-only: lumi 36.4% (49.4%), nPDF 5.6%, muSF 4.7% (8.1%), qcdScale 4.0%, qcd_sB/qcd_rate 2.5–3.0%
- e-only: lumi 28.3% (38.8%), **eleSF 12.0% (33.6%)**, qcd_sB_ele 5.7/5.6% (23.4/22.3% on r_Z), qcd_rate_ele 4.2/4.3% (max 15.8%, r_Wm_y11), nPDF 3.9%, qcdFF_ele 2.2/2.8%

## 4. Open issues this iteration left

1. **Z→ee normalization:** σ_Z e/μ = 1.27 ± 0.11. The pO Reco tag-and-probe
   (`correction/reco_sf_ele.C`) measured SF_reco = 0.988 +0.018 −0.023 (stat)
   ± 0.002 against the applied pp 0.929 for the same probes — the pp factor
   over-corrects each pO electron by ~6%. → addressed by the 2026-09-30 iteration.
2. **Electron charge split:** σ(W⁺) e/μ 1.065 vs σ(W⁻) 0.952 (double ratio
   1.12 ± 0.06); e/μ(W⁻) low in y0–y3 and y11, where the e⁻ QCD fraction is
   largest — suspected y-split of the y-inclusive QCD template.

## 5. Backups of this iteration (taken 2026-09-30)

| path | content |
|---|---|
| `skim/rootfile_pre_recotnp/` | the 14 electron skim files (`WToElecNu_*`, `ZToEE_*`) |
| `skim/rootfile_pre_recotnp/observables/` | `{fidyields,charge_asym_fid,FBratio_fid}_{comb,simfit_mu,simfit_ele}_leppt_mt40.root` |
| `skim/logs/pre_recotnp/` | `skim_{Wel,Zee}_*.log`, `lhe_updown_{Wel,Zee}.log` |
| `skim/output/pre_recotnp/` | `WToElecNu_pO_PFMet_*.txt` cutflows |
| `correction/rootfile/pre_recotnp/` | `qcd_abcd_{ele,mu}.root`, `reco_sf_ele/skim_{Data,DY}.root` |
| `correction/logs/pre_recotnp/` | `qcd_abcd_{ele,mu}.log`, `reco_sf_ele_*.log`, `njet_WZ_{Wel,Zee}.log` |
| `correction/plots/pre_recotnp/` | `reco_sf_ele.csv` |
| `plotting/plots/pre_recotnp/{Elec,mu}/` | all Combine inputs + `_systs.txt` sidecars |
| `plotting/plots/pre_recotnp/{comb,flavfit}/` | every observables plot and CSV of this fit |
| `plotting/logs/pre_recotnp/` | `mtandmet_*`, `dileptonpeak_*`, `syst_shapes_*` logs |
| `analysis/logs/pre_recotnp/` | `run_observables_leppt_mt40.log` |
| fork `test/pO_fit_out_leppt_mt40_pre_recotnp/` | the downloaded fit trees `simfit/`, `simfit_mu/`, `simfit_ele/` (summaries, cards, fitDiagnostics, contours, impacts, logs) |
