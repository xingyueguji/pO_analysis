# pO_analysis — end-to-end runbook

W and Z boson cross-section measurement in proton–oxygen (pO) collisions.

Four physics channels — W→μν, W→eν, Z→μμ, Z→ee — each producing fitted signal
yields, charge asymmetries vs rapidity, and forward/backward (F/B) ratios in
|y_CM|. Tau channels (W→τν, Z→ττ) are treated as backgrounds.

This is the **single runbook for the whole procedure**, across **two repos**:

| repo | role | needs |
|------|------|-------|
| `pO_analysis` (this repo) | skim → MC norm → plotting → **structured Combine inputs** → final observables | ROOT 6 (+ACLiC); tested 6.32 |
| `HiggsAnalysis-CombinedLimit` fork, branch `zheng/po-analysis` | the Combine **fit**: datacards, FitDiagnostics, fitted-yield extraction, postfit plots | **`cmsenv`** (combine + text2workspace.py) |

Only the fit (Module 4) needs `cmsenv`; everything else is plain ROOT. See
[CLAUDE.md](CLAUDE.md) for the architecture (what each macro does, why histogram
naming is load-bearing, inter-channel asymmetries). The pipeline:

```
skim → ngen → ABCD QCD → structured Combine inputs (mtandmet/dileptonpeak)
     → fit (fork: run_pO_fits.sh)  → <chan>_fitted_yields.root
     → charge_asym.C / FBratio.C   → observables.C (per-channel + merged overlay)
```

## TL;DR (full chain)

```bash
# ---- pO_analysis (plain ROOT) ----
cd skim        && ./run_all.sh all && ./run_ngen.sh                          # 1,2 skims + N_gen
cd ../correction && root -l -b -q 'qcd_abcd.C+' && root -l -b -q 'qcd_abcd.C+(true)'   # 3a ABCD QCD (mu,ele)
cd ../plotting && for a in 'mtandmet.C+(false)' 'mtandmet.C+(true)' \
                          'dileptonpeak.C+(false)' 'dileptonpeak.C+(true)' \
                          'plotRpOtheory.C+'; do root -l -q -b "$a"; done       # 3b inputs + theory

# ---- fork (cmsenv) -- locally, or push to lxplus (see Module 4) ----
cd ../../HiggsAnalysis-CombinedLimit/test && cmsenv && ./run_pO_fits.sh both all   # 4 fit

# ---- pO_analysis (plain ROOT) ----
cd ../../pO_analysis/analysis                                                    # 5 observables
for c in mu ele; do F=../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/$c/summary/${c}_fitted_yields.root
  root -l -b -q "charge_asym.C+(\"$F\",\"../skim/rootfile/charge_asym_fit_$c.root\")"
  root -l -b -q "FBratio.C+(\"$F\",\"../skim/rootfile/FBratio_fit_$c.root\")"; done
cd ../plotting && root -l -b -q 'observables.C+(false)' && root -l -b -q 'observables.C+(true)' \
  && root -l -q -b -e 'gROOT->LoadMacro("observables.C+"); observables_overlay();'
```

## Layout

| Directory          | Stage                                  | Entry point |
| ------------------ | -------------------------------------- | ----------- |
| `merge_rootfile/`  | (one-time) discover/hadd EOS ntuples   | `make_filelist.sh`, `hadd_from_list.sh` |
| `skim/`            | selection → per-sample histos; N_gen   | `run_all.sh`, `run_ngen.sh` |
| `correction/`      | ABCD QCD, isolation WP, Data/MC checks | `qcd_abcd.C`, … |
| `plotting/`        | data/MC overlays + **Combine inputs**  | `mtandmet.C`, `dileptonpeak.C`, `plotRpOtheory.C` |
| `analysis/`        | charge asymmetry, F/B ratio            | `charge_asym.C`, `FBratio.C` |
| (fork) `test/`     | Combine fit pipeline                   | `run_pO_fits.sh`, `sync_lxplus.sh` |

## Prerequisites

- ROOT 6 with ACLiC (`.C+`).
- Input ntuples: per-sample May-26 files; path single-sourced in
  `skim/skim_common.h` (`kDefaultDataFile` + `ResolveMCSample`), currently the
  **local** `~/pO_2026_May_26/…`. Repoint there to read off EOS/lxplus.
- For the fit: `cmsenv` (combine + text2workspace.py) and the fork checked out to
  `zheng/po-analysis`.

---

# Workflow

## Module 1 — input data (`merge_rootfile/`, one-time per data version)

Per-sample ROOT files (the May-26 production): `Data_May_26.root` + 9 MC. The
skim reads them directly — `merge_rootfile/` only matters when (re)building from
EOS. Discover + merge:

```bash
cd merge_rootfile/
./make_filelist.sh        # discover EOS files -> per-sample .txt lists
./hadd_from_list.sh       # hadd into the per-sample ROOT files
```

## Module 2 — skim (`skim/`)

8-step W cutflow / Z selection → rapidity-binned MT & MET histograms per sample.

```bash
cd skim/
grep '^DATA_FILE=' run_all.sh           # 2.1 sanity-check the input path
./run_all.sh Wmu Data                   # 2.2 smoke-test one channel×sample
./run_all.sh all                        # 2.3 everything (Zmm Zee Wmu Wel × 7 samples)
```
CLI: `./run_all.sh <Zmm|Zee|Wmu|Wel|all> [samples…]` (samples ⊂
`Data DY Wp Wm DYtau Wptau Wmtau`).

**Check:** `skim/logs/*.log` (per-job OK/FAIL), `skim/output/*.txt` (cutflow),
`skim/rootfile/*_hist.root`. W files hold `h_{mt,met}_W{p,m}_y0..11` (+ `_FB`)
and the ABCD planes `h_iso_{met,mt}_{mu,ele}{Plus,Minus}`; Z files hold `hMass*`
+ Z kinematics + recoil histos.

## Module 2b — MC normalization (`skim/run_ngen.sh` + `skim/mc_norm.h`) — DONE

Every MC sample is scaled (downstream, not in the skim) by
`k_s = A·σ·L/N_gen` so templates sit at their **absolute** pO yield. A=16
(Oxygen A-scaling), L=46.5 nb⁻¹, σ read from the POWHEG weights (`⟨w⟩=σ`:
W⁺ 6.376, W⁻ 5.464, DY 1.175 nb).

```bash
cd skim/
./run_ngen.sh        # Σ gen-weight over ALL events -> rootfile/ngen.root + output/ngen.txt
```

`pONorm::MCScale("Wp_mu")` (in `mc_norm.h`) reads `ngen.root` and is **wired into**
`mtandmet.C` + `dileptonpeak.C`, which set `ps.normBkgToData=false` → stacks drawn
ABSOLUTE (no area norm). Escape hatch: if absolute MC is ~16× off, set `kA_O=1.0`.

## Module 3 — plotting + structured Combine inputs (`plotting/`, `correction/`)

```bash
# 3a. Data-driven ABCD QCD templates (low-MET background). FROM correction/.
cd correction/
root -l -b -q 'qcd_abcd.C+'                 # muon     -> rootfile/qcd_abcd_mu.root
root -l -b -q 'qcd_abcd.C+(true)'           # electron -> rootfile/qcd_abcd_ele.root

# 3b. Data/MC MT+MET plots AND the structured Combine inputs. FROM plotting/.
cd ../plotting/
root -l -b -q 'mtandmet.C+(false)'          # muon  -> plots/combine_input_W.root
root -l -b -q 'mtandmet.C+(true)'           # ele   -> plots/Elec/combine_input_W.root
root -l -b -q 'dileptonpeak.C+(false)'      # Zmm   -> plots/combine_input_Z.root
root -l -b -q 'dileptonpeak.C+(true)'       # Zee   -> plots/Elec/combine_input_Z.root

# 3c. Theory predictions for R_FB (producer; reads pQCDLightIon/ + filelist_theory.txt).
root -l -b -q 'plotRpOtheory.C+'            # -> RpO_rootfile/RpO_FB_graphs.root
```

`combine_input_W.root` is **structured**: one TDirectory per fit region
(`Wp_lab_y0..11`, `Wm_lab_y*`, `Wp_fb_y*`, `Wm_fb_y*`, `Wp_incl`, `Wm_incl`,
`W_incl`), each with the 6 **absolute** templates `data_obs/signal/z/ztau/wtau/qcd`
(MET discriminant; per-y ABCD QCD). `combine_input_Z.root` has a `Z_incl/` dir
(`data_obs/signal/w/wtau/ztau`, mass peak).

Optional/cosmetic: `mtandmet_overlay.C`, `plotZcurve.C`. Scratch (ignore):
`test.C`, `test111.C`, `Z_MC_overlay.C`.

## Module 4 — Combine fit (fork `zheng/po-analysis`, needs `cmsenv`)

One driver does everything per channel. Full details:
`HiggsAnalysis-CombinedLimit/test/README_pO_fits.md`.

```bash
cd HiggsAnalysis-CombinedLimit/test
cmsenv
./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|all] [--dry-run] [--no-postfit] [--draw-only]
```

It (1) finds this repo's structured inputs (env `PO_PLOTS` / `--plots-dir`, else
autodetect), (2) generates **all** datacards (`make_pO_datacards.sh`: 48 per-(charge,y)
lab+FB — each a **two-channel card fitted simultaneously with `Z_incl`** —
per-charge incl, `W_incl`, `Z_incl`, simultaneous `WZ`), (3) runs
`text2workspace` + `combine -M FitDiagnostics` per region into a clean tree
`pO_fit_out/<chan>/{datacards,fits/<region>,postfit,summary}`, (4) extracts
fitted yields (`extract_pO_yields.C`), (5) draws postfit plots (`draw_postfit_pO.C`,
same cosmetics as `mtandmet.C`).

**Fit model (two-parameter, 2026-07-01):** two MC scales per fit — the POI
**`r` = all W-related MC** (W `signal` + `wtau`; plus the `w`/`wtau` backgrounds
under the Z peak in simultaneous cards) and **`dy_norm` = all DY-related MC**
(`z` + `ztau`; plus the Z signal in simultaneous cards). Composition WITHIN each
group stays locked by the absolute `k_s` templates. Standalone `Z_incl` card:
roles flip — the POI `r` IS the DY scale (`signal`+`ztau`), the W backgrounds
get a free `w_norm`. Data-driven ABCD `qcd` → its own free `qcd_norm`.
**Simultaneous cards** (`WZ` and every per-bin W card) have two fit channels
(W region + `Z_incl`): `r` scales the W-related in both, the shared `dy_norm`
scales the DY-related in both — the high-purity Z peak pins it (this replaces
the old shared `eff_lumi`). Systematics deliberately deferred (stat-only fits).

**Outputs** (`pO_fit_out/<chan>/summary/`): `<chan>_W_yields.csv`,
`<chan>_summary.csv`, and `<chan>_fitted_yields.root` — single-bin
`h_mt_W{p,m}_y0..11` (lab) + `..._FB` histos with Sumw2 = fit error, i.e. the
exact names `charge_asym.C`/`FBratio.C` read.

Channel = `mu|ele|both`; mode = `perbin` (48 per-(charge,y) W regions, each
simultaneous with `Z_incl`), `incl` (`Wp_incl Wm_incl W_incl Z_incl`,
standalone), `combined` (the `WZ` fit only), or `all`. `--dry-run` builds
datacards without `cmsenv`; `--no-postfit` skips plots; `--draw-only` redraws
the postfit plots from an existing fit run (only `root` needed — for cosmetic
`draw_postfit_pO.C` changes; requires the `fits/` tree, so redraw where the
fits ran and `sync_lxplus.sh download --postfit`).
Per-region failures (e.g. a tail FB bin with an all-zero `qcd` template) are
logged and skipped, not fatal — check `pO_fit_out/<chan>/fits/<region>/fit.log`.
Sanity-check `dy_norm ≈ 1` in `<chan>_summary.csv` after the fits (`w_norm ≈ 1`
for `Z_incl`). Each postfit plot also carries on-plot fit-quality: `χ²/ndf (p)`
(Poisson GoF), the fitted `r`, the fitted `DY norm`/`W norm`/`QCD norm`
(whichever float in that fit), and a red `status/covQ` flag if the fit didn't
converge cleanly (see the fork README "Diagnosing fit quality").

### Running the fit on lxplus (split workflow)

Build inputs locally (Modules 1–3), fit on lxplus, run observables locally. The
helper `test/sync_lxplus.sh` wraps every transfer over ONE SSH connection
(single auth; `kinit zheng@CERN.CH` first for none):

```bash
cd HiggsAnalysis-CombinedLimit/test
./sync_lxplus.sh upload                  # 4 inputs + scripts -> lxplus
ssh zheng@lxplus.cern.ch                 # then: cmsenv; cd $FORK_LX/test
#   PO_PLOTS=/afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots ./run_pO_fits.sh both all
./sync_lxplus.sh download                # summary/ (fitted yields + CSVs) <- lxplus
./sync_lxplus.sh download --postfit      # also the postfit plots
```
lxplus paths (override via env): `ANA_LX=/afs/cern.ch/user/z/zheng/pO_analysis`,
`FORK_LX=/afs/cern.ch/user/z/zheng/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit`.

If a downloaded `<chan>_fitted_yields.root` is ever empty but the CSV is fine,
rebuild it without re-running the fit:
`root -l -b -q 'my_script/make_yields_from_csv.C("<chan>_W_yields.csv","<chan>_fitted_yields.root")'`.

## Module 5 — final observables (`analysis/` + `plotting/`)

Run `charge_asym.C` / `FBratio.C` on the **fitted-yields** file (they take the
input as their first arg — fitted yields replace raw, no edits), then plot.

```bash
cd analysis/
for c in mu ele; do
  F=../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/$c/summary/${c}_fitted_yields.root
  root -l -b -q "charge_asym.C+(\"$F\",\"../skim/rootfile/charge_asym_fit_$c.root\")"
  root -l -b -q "FBratio.C+(\"$F\",\"../skim/rootfile/FBratio_fit_$c.root\")"
done

## one line version: 
cd /Users/zhenghuang/pO_analysis/analysis && for c in mu ele; do F=../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/$c/summary/${c}_fitted_yields.root; root -l -b -q "charge_asym.C+(\"$F\",\"../skim/rootfile/charge_asym_fit_$c.root\")"; root -l -b -q "FBratio.C+(\"$F\",\"../skim/rootfile/FBratio_fit_$c.root\")"; done

cd ../plotting/
root -l -b -q 'observables.C+(false)'       # muon individual   -> plots/{charge_asym,FBratio}/
root -l -b -q 'observables.C+(true)'        # electron individual -> plots/Elec/...
# merged muon+electron overlay -> plots/merged/ :
root -l -q -b -e 'gROOT->LoadMacro("observables.C+"); observables_overlay();'

## one line version
cd /Users/zhenghuang/pO_analysis/plotting && root -l -b -q 'observables.C+(false)' && root -l -b -q 'observables.C+(true)' && root -l -b -q -e 'gROOT->LoadMacro("observables.C+"); observables_overlay();'

```

`observables.C` overlays the data with **all four** nPDF theory bands
(EPPS21/nCTEQ15HQ/nNNPDF3.0/TUJU21nlo, drawn as filled bands only — no central
line / error bars). The merged overlay shows muon (black circles) + electron
(red squares) on the same axes; the sum-channel theory is weighted by the
combined μ+e fitted yields. Theory file optional (missing → data-only).
Fiducial W cross sections (W⁺, W⁻, W inclusive; μ and e overlaid) straight from
the summary CSVs:
```bash
root -l -b -q 'xsec_fiducial.C+'        # W+/W-/W incl (mu+e) -> plots/xsec/W_fiducial.{png,pdf}
#   sigma_fid = N_fit / L (L = pONorm::kLumi_invnb = 46.5 nb^-1).
# differential vs CM-frame lepton pseudorapidity, d(sigma_fid)/d(eta), per channel
# (uses the per-bin <chan>_W_yields.csv; eta bins taken from g_chargeAsym_mt):
root -l -q -b -e 'gROOT->LoadMacro("xsec_fiducial.C+"); xsec_fiducial_diff(false); xsec_fiducial_diff(true);'
#   -> plots/xsec/W_dsigma_deta_{mu,ele}.{png,pdf}
```
NB this is the fiducial σ **before** the lepton-efficiency correction (= σ_fid×ε),
which is why μ and e differ — the gap is the channel efficiency, and they should
converge once ε is applied (universality). Stat (fit) uncertainty only.
For arbitrary regions, read the absolute yields directly from the CSVs.

## Module 6 — corrections & studies (`correction/`)

Run from `correction/`; outputs in `correction/plots/`, `correction/rootfile/`.

```bash
cd correction/
root -l -b -q 'qcd_abcd.C+'                       # ABCD QCD (muon)  [also Module 3a]
root -l -b -q 'qcd_abcd.C+(true)'                 # ABCD QCD (electron)
root -l -b -q 'isolation_mu_tight.C+("<data.root>")' # muon iso study (TightID, Δβ relIso) [current]
root -l -b -q 'isolation_ele.C+("<data.root>")'   # electron iso/ID ROC study (Δβ relIso)
root -l -b -q 'PlotsIsoROC.C+(false)'             # / PlotIsoROC_ele.C -> ROC plots
root -l -b -q 'plot_iso_summary.C+'               # per-ID iso summary (reads both studies)
# isolation.C = legacy muon multi-cone/multi-def scan (uncorrected iso, kept for reference)
root -l -b -q 'dataMC_kinematics.C+("Zmm")'       # Data/MC Z kinematics (Zmm/Zee)
root -l -b -q 'recoil_raw.C+'                      # raw hadronic recoil (recoil_raw_ele.C for e)
```
`qcd_sideband_fit_and_extrapolate.C` is the **superseded** Rayleigh QCD method
(kept for reference; ABCD replaced it).

---

## Notes / health of the pipeline

- **MC normalization is wired (absolute).** `k_s = A·σ·L/N_gen` applied in
  `mtandmet.C` + `dileptonpeak.C`; templates and Combine inputs are absolute (no
  area norm). The Z peak validates it: Z `signal` lands within ~4% of data.
- **Combine inputs renamed.** `combine_input_W.root` / `combine_input_Z.root`
  (structured, TDir per region) replace the old `combine_input_inclusive.root` /
  `combine_input_dilepton.root`. The old fork pipeline
  (`run_fit.sh`, `make_combine_input*.C`, `testdatacard_*.txt`,
  `draw_postfit_{inclusive,Zmumu,Zee}.C`) is **superseded** by `run_pO_fits.sh`.
- **Sumw2 / weights.** All skim histos `Sumw2()`'d; per-event gen weight applied
  to MC; `analysis_helpers.h` Sumw2-aware (`YieldInRange`→`IntegralAndError`,
  `AsymErr`/`RatioErr` propagate σ²). Fitted-yield histos carry the fit error in
  Sumw2, so charge_asym/FBratio errors are the fit uncertainties.
- **Corrections WIP.** Recoil, lepton SFs, momentum scale/smearing not all
  applied — don't assume MC in `skim/rootfile/` is fully corrected.
- **Acceptance / efficiency.** The fitted yields are reconstructed, in-acceptance
  yields; the asymmetry / F/B are fiducial, lepton-level observables. Acceptance
  is only needed to extrapolate to total cross-sections or boson-level theory.
- **Pre-existing inter-channel asymmetries** (intentional; see CLAUDE.md):
  DY-veto pT 15 (μ) vs 10 (e) GeV; isolation 0.15 (μ) vs 0.095 (e) — both
  re-confirmed optimal under the Δβ-corrected relIso (2026-07-06, MuonPOG
  convention, all channels; needs re-skim). All electron ID gates are
  `eleMVAIdWP95` since 2026-07-02 (`skim_Zee` included).
- **No tests/CI.** Validate by re-running + inspecting logs/plots. Commit
  messages are often "xx" — read diffs, not `git log`.
