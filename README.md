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
     → fit (fork: run_pO_fits.sh; DEFAULT = simfit, the μ+e GRAND SIMULTANEOUS
            fit → simfit/summary/comb_fitted_yields.root [+ covariance];
            legacy per-flavour modes → <chan>_fitted_yields.root)
     → charge_asym.C / FBratio.C   → observables.C (comb = primary;
            legacy per-channel + merged overlay kept for comparison)
```

## TL;DR (full chain)

```bash
# ---- pO_analysis (plain ROOT) ----
cd skim        && ./run_all.sh all && ./run_ngen.sh                          # 1,2 skims + N_gen
cd ../correction && ./run_qcd_abcd.sh                                        # 3a ABCD QCD (mu+ele, logged)
cd ../plotting && for a in 'mtandmet.C+(false)' 'mtandmet.C+(true)' \
                          'dileptonpeak.C+(false)' 'dileptonpeak.C+(true)' \
                          'plotRpOtheory.C+'; do root -l -q -b "$a"; done       # 3b inputs + theory

# ---- fork (cmsenv) -- locally, or push to lxplus (see Module 4) ----
cd ../../HiggsAnalysis-CombinedLimit/test && cmsenv && ./run_pO_fits.sh --asimov   # 4 fit (simfit = DEFAULT; --asimov adds the closure fit)
#   ./run_pO_fits.sh both all   # legacy per-flavour per-bin fits + simfit

# ---- pO_analysis (plain ROOT) ----
cd ../../pO_analysis/analysis && ./run_observables.sh                            # 5 observables (met)
#   ./run_observables.sh leppt | leppt_mt40 | all    -> per-discriminant folders (Module 5)
```

## Layout

| Directory          | Stage                                  | Entry point |
| ------------------ | -------------------------------------- | ----------- |
| `merge_rootfile/`  | (one-time) discover/hadd EOS ntuples   | `make_filelist.sh`, `hadd_from_list.sh` |
| `skim/`            | selection → per-sample histos; N_gen   | `run_all.sh`, `run_ngen.sh` |
| `correction/`      | ABCD QCD, isolation WP, Data/MC checks | `run_qcd_abcd.sh` → `qcd_abcd.C`, … |
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
#     ALWAYS go through the wrapper -- it keeps the log (see the note below).
cd correction/
./run_qcd_abcd.sh                           # both channels -> rootfile/qcd_abcd_{mu,ele}.root
#   ./run_qcd_abcd.sh mu                    # one channel only

# 3b. Data/MC MT+MET plots AND the structured Combine inputs. FROM plotting/.
cd ../plotting/
root -l -b -q 'mtandmet.C+(false)'          # muon  -> plots/combine_input_W.root
root -l -b -q 'mtandmet.C+(true)'           # ele   -> plots/Elec/combine_input_W.root
root -l -b -q 'dileptonpeak.C+(false)'      # Zmm   -> plots/combine_input_Z.root
root -l -b -q 'dileptonpeak.C+(true)'       # Zee   -> plots/Elec/combine_input_Z.root

# 3c. Theory predictions for R_FB (producer; reads pQCDLightIon/ + filelist_theory.txt).
root -l -b -q 'plotRpOtheory.C+'            # -> RpO_rootfile/RpO_FB_graphs.root
```

**Why the wrapper and not `root -l -b -q 'qcd_abcd.C+'` directly:** the macro's
console output IS a deliverable, not chatter. It carries the region composition
and closure, the factorisation tests (`T` in sub-slices of the low-y band), the
anti-iso window scan, the two-plane transport with its r-scan, the multijet
fraction of every fitted selection, and the assembled systematic budget → the
lnN κ used in the datacards — plus, since 2026-08-23, the **in-fit ABCD block**
(the m_T-plane B/C40/D counts + A40 prediction exported as `abcd_counts_*`, the
reduced κ for `QCD_MODE=abcd`) and the **per-pT-bin fake-factor diagnostic**
(F(pT) tables + `ff_*` plots). Those are the numbers quoted in
[AN_qcd_background.tex](AN_qcd_background.tex), and ROOT prints them to stdout
only. `run_qcd_abcd.sh` pre-builds once (so the two channels cannot race on the
ACLiC artifacts), tees each channel to `correction/logs/qcd_abcd_<chan>.log`
(~290 lines) and echoes the transfer factors and the κ's to the terminal.
Takes ~2 s per channel, so regenerate freely — but never run the macro bare and
lose the report.

`combine_input_W.root` is **structured**: one TDirectory per fit region
(`Wp_lab_y0..11`, `Wm_lab_y*`, `Wp_fb_y*`, `Wm_fb_y*`, `Wp_incl`, `Wm_incl`,
`W_incl`), each with the 6 **absolute** templates `data_obs/signal/z/ztau/wtau/qcd`
(MET discriminant; per-y ABCD QCD). `combine_input_Z.root` has a `Z_incl/` dir
(`data_obs/signal/w/wtau/ztau`, mass peak). **`combine_input_W_leppt_mt40.root`
additionally carries the in-fit-ABCD objects (2026-08-23):** a 7th SR template
`qcd_abcd` (same shape, total = B0·C40/D0) and 6 CR dirs `{Wp,Wm}_CR{B,C,D}`
of 1-bin templates — consumed only by `QCD_MODE=abcd` cards; run 3a BEFORE 3b
or the writer warns and skips them.

Optional/cosmetic: `mtandmet_overlay.C`, `plotZcurve.C`. Scratch (ignore):
`test.C`, `test111.C`, `Z_MC_overlay.C`.

## Module 4 — Combine fit (fork `zheng/po-analysis`, needs `cmsenv`)

One driver does everything per channel. Full details:
`HiggsAnalysis-CombinedLimit/test/README_pO_fits.md`.

```bash
cd HiggsAnalysis-CombinedLimit/test
cmsenv
./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|simfit|all] [--disc met|leppt|leppt_mt40] [--dry-run] [--no-postfit] [--draw-only] [--asimov]
```

### The grand simultaneous fit — `simfit` (2026-08-04, the DEFAULT)

`./run_pO_fits.sh` with no arguments runs **one likelihood per binning variant
(lab, fb)**: all 48 W channels ({μ,e} × {W⁺,W⁻} × y0..11) **plus both
Z-inclusive peaks**, with 2N+1 = **25 POIs**:

- **`r_<C>_y<i>`** (24) — the W signal strength of rapidity bin *i*, charge *C*,
  scaling that bin's W-related MC (`signal` + `wtau`) in the muon AND electron
  channels (**μ/e shared**: lepton universality, with the relative μ/e
  acceptance×efficiency taken from MC — note lepton SFs are not applied yet).
- **`r_Z`** (the "+1") — ONE global scale on all DY-related MC: `z` + `ztau` in
  every W channel and the DY signal (+`ztau`) under both Z peaks. The DY
  rapidity dependence across W bins comes fixed from MC; only this global
  normalization floats, pinned jointly by the two Z peaks.
- QCD (data-driven ABCD templates) — three modes via `QCD_MODE`:
  **`lnN` (default)**: one log-normal nuisance per (flavour, charge),
  `qcd_rate_{mu,ele}_{Wp,Wm}`, κ μ 1.15 / e 1.20 at the ABCD prediction
  (2026-08-17); **`free`**: the pre-2026-08-17 48 per-channel `qcd_norm`
  rateParams; **`abcd` (2026-08-23, `--disc leppt_mt40` only)**: the IN-FIT
  ABCD — 12 counting CR channels (`<F>_<C>_CR{B,C,D}`, m_T plane) with free
  scales `qcd_s{B,C,D}_<F>_<C>`, the SR `qcd_abcd` template scaled by the
  formula rateParam `(sB·sC/sD)`, EWK subtraction riding the POIs (CRB DY →
  r_Z, CRB W → per-y `w_y*` mapped to the r's; `QCD_WCR=frozen` freezes it),
  and the reduced residual κ μ 1.09 / e 1.15 (`QCD_ABCD_LNN_MU/ELE`).
- `w`/`wtau` under the Z peaks — **frozen at absolute MC** (0.03–0.06 events
  under 372/252-event peaks; decision 2026-08-04).

lab and fb are the same events rebinned, so they are fitted **separately** (two
workspaces, same POI names; charge asymmetry ← lab, R_FB ← fb). This replaces
the legacy scheme's statistical flaw — 48 per-bin fits per flavour each
re-using the same Z data with the induced correlations ignored — with one
correct likelihood, and it produces the full covariance of the r's:
`extract_pO_simfit.C` stores it as `h_cov_yield[_FB]` and
`charge_asym.C`/`FBratio.C` automatically include the cross terms when present.

Implementation: `make_pO_simfit_cards.sh` writes one 50-channel datacard + the
`multiSignalModel` map file per variant (plain files — `--dry-run` needs no
`cmsenv`); then `text2workspace.py -P ...:multiSignalModel` and `combine -M
FitDiagnostics --skipBOnlyFit` per variant. `--asimov` adds a prefit-Asimov
closure fit (`-t -1`): every fitted POI must come back at 1 — checked and
reported (PASS/FAIL) by the extraction. Outputs under
`pO_fit_out<suffix>/simfit/`: `summary/comb_W_yields.csv`, `comb_summary.csv`
(all POIs, fit status/covQual, covariance-propagated W⁺/W⁻/W inclusive sums,
Asimov closure rows) and **`comb_fitted_yields.root`** — the `h_yield_*` names
Module 5 reads (yields = r × the μ+e-summed prefit template integral) plus the
covariance matrices. Postfit plots for every channel of the grand fit land in
`simfit/postfit/` (info box: that bin's `r_<C>_y<i>`, the global `r_Z` shown as
"DY norm", the channel's `qcd_norm`).

### Legacy per-flavour pipeline (`perbin|incl|combined` — the superseded scheme, kept runnable for comparison)

It (1) finds this repo's structured inputs (env `PO_PLOTS` / `--plots-dir`, else
autodetect), (2) generates **all** datacards (`make_pO_datacards.sh`: 48 per-(charge,y)
lab+FB — each a **two-channel card fitted simultaneously with `Z_incl`** —
per-charge incl, `W_incl`, `Z_incl`, simultaneous `WZ`), (3) runs
`text2workspace` + `combine -M FitDiagnostics` per region into a clean tree
`pO_fit_out/<chan>/{datacards,fits/<region>,postfit,summary}`, (4) extracts
fitted yields (`extract_pO_yields.C`), (5) draws postfit plots (`draw_postfit_pO.C`,
same cosmetics as `mtandmet.C`). Mode `all` = this whole legacy pipeline PLUS
the simfit (when channel = `both`).

### W discriminant variants (`--disc`, 2026-07-30)

The W fit can run on three discriminants; `--disc` selects which (default
`met`, so all existing usage is unchanged):

| `--disc` | discriminant | selection | input file (per channel) | output tree |
|---|---|---|---|---|
| `met` | PF MET shape | plain W selection | `combine_input_W.root` | `pO_fit_out/` |
| `leppt` | lepton pT | plain W selection | `combine_input_W_leppt.root` | `pO_fit_out_leppt/` |
| `leppt_mt40` | lepton pT | pT>25 && m_T>40 | `combine_input_W_leppt_mt40.root` | `pO_fit_out_leppt_mt40/` |

Datacards and the fit model are IDENTICAL for all three (same region names,
same two-parameter scheme, `qcd_norm` FREE for all — 2026-07-30 decision);
only the input file, the output tree and the postfit x-title change.

> **WARNING — carry the SAME `--disc` through the ENTIRE workflow of a variant
> run.** The out-trees carry no marker of which discriminant produced them —
> the coupling is only the directory suffix — so mixing steps corrupts or
> mislabels silently:
> - **fit**: `./run_pO_fits.sh both all --disc leppt_mt40`
> - **redraw**: `--draw-only` MUST repeat the same `--disc` (it selects the
>   out-tree AND the axis title; without it, lepton-pT plots get relabeled
>   "PF MET (GeV)" with no error). Same rule if you ever use `--out`.
> - **download**: `sync_lxplus.sh download` needs NO flag — it sweeps all
>   three out-trees automatically, skipping absent ones.
> - **observables (Module 5)**: run `analysis/run_observables.sh <disc>` —
>   it carries the tag through every step automatically (reads the matching
>   `pO_fit_out<suffix>/` tree, writes disc-tagged graph files and per-disc
>   plot folders, stamps the discriminant on every plot). The histogram
>   names inside the yields files are identical across variants
>   (`h_yield_*`), so if you ever drive the macros by hand instead, the tree
>   you point at is the ONLY thing distinguishing a MET result from a pT
>   result — the driver exists so you never have to get that right manually.
> - **inputs**: `mtandmet.C` writes all three files in one run; a variant
>   whose skim histograms are missing is skipped AND its stale file deleted,
>   so a missing `combine_input_W_leppt*.root` means "re-run Module 3", never
>   "use the old one".
>
> Physics note for the pT variants: with no low-MET region in the
> discriminant, the fit constrains `qcd_norm` only weakly — check the fitted
> `QCD norm` on the postfit plots / in `<chan>_summary.csv` before trusting
> the composition (if it wanders far from 1, consider constraining it to the
> ABCD prediction instead of leaving it free).

**Legacy fit model (two-parameter, 2026-07-01):** two MC scales per fit — the POI
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

Channel = `mu|ele|both`; mode = `simfit` (DEFAULT, see above; ignores the
channel argument — always μ+e), `perbin` (48 per-(charge,y) W regions, each
simultaneous with `Z_incl`), `incl` (`Wp_incl Wm_incl W_incl Z_incl`,
standalone), `combined` (the `WZ` fit only), or `all` (legacy + simfit).
`--dry-run` builds
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
./sync_lxplus.sh upload                  # inputs (4 required + pT variants if built) + scripts
ssh zheng@lxplus.cern.ch                 # then: cmsenv; cd $FORK_LX/test
#   # simfit (DEFAULT) + Asimov closure:
#   PO_PLOTS=/afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots ./run_pO_fits.sh --asimov
#   # legacy per-flavour pipeline + simfit together:
#   PO_PLOTS=... ./run_pO_fits.sh both all
#   # discriminant variants (SEE THE --disc WARNING above -- keep the flag
#   # consistent for every later step of that variant's workflow):
#   PO_PLOTS=... ./run_pO_fits.sh --asimov --disc leppt
#   PO_PLOTS=... ./run_pO_fits.sh --asimov --disc leppt_mt40
./sync_lxplus.sh download                # summary/ <- lxplus, ALL out-trees (met + variants)
./sync_lxplus.sh download --postfit      # also the postfit plots
```
lxplus paths (override via env): `ANA_LX=/afs/cern.ch/user/z/zheng/pO_analysis`,
`FORK_LX=/afs/cern.ch/user/z/zheng/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit`.

If a downloaded `<chan>_fitted_yields.root` is ever empty but the CSV is fine,
rebuild it without re-running the fit:
`root -l -b -q 'my_script/make_yields_from_csv.C("<chan>_W_yields.csv","<chan>_fitted_yields.root")'`.

## Module 5 — final observables (`analysis/` + `plotting/`)

**One command per discriminant (2026-08-03):** `analysis/run_observables.sh`
runs the whole chain, carrying the disc tag through every filename and output
folder so the three variants coexist without overwriting each other. Since
2026-08-04 it has two conditional blocks, each run only when its fit outputs
exist: the **PRIMARY simfit chain** (`charge_asym.C` + `FBratio.C` on
`simfit/summary/comb_fitted_yields.root` — errors include the fit covariance —
then `observables_comb`), and the **legacy per-flavour chain** (both channels'
fitted yields, per-channel / overlay / fiducial-σ plots, kept for comparison):

```bash
cd analysis/
./run_observables.sh              # met (default) -- the PF-MET-shape fit
./run_observables.sh leppt        # lepton-pT variant   (its fit must exist)
./run_observables.sh leppt_mt40   # lepton-pT, mT>40 variant
./run_observables.sh all          # every variant whose out-tree exists (others SKIP)
```

Outputs, per `<disc>` = `met` | `leppt` | `leppt_mt40`:

| output | path |
|---|---|
| **PRIMARY: simfit (μ+e comb)** graphs | `skim/rootfile/{charge_asym,FBratio}_fit_comb_<disc>.root` |
| **PRIMARY: simfit (μ+e comb)** plots | `plotting/plots/comb/{charge_asym,FBratio}/<disc>/` |
| **PRIMARY: simfit** fiducial σ: post-fit vs reco-MC vs **gen-MC** | `plotting/plots/comb/xsec/<disc>/` |
| **PRIMARY: simfit** y-inclusive postfit stacks (μ/e × W⁺/W⁻/W) | `plotting/plots/comb/postfit_incl/<disc>/` |

For the generator-level overlay on the comb σ plots, produce the (one-time,
discriminant-independent) gen histograms first: `cd skim && root -l -b -q
'gen_xsec.C+'` → `skim/rootfile/gen_xsec.root` (missing file ⇒ the overlay is
skipped with a note, everything else unaffected).
| legacy graph files (charge asym, R_FB) | `skim/rootfile/{charge_asym,FBratio}_fit_{mu,ele}_<disc>.root` |
| legacy per-channel plots | `plotting/plots[/Elec]/{charge_asym,FBratio}/<disc>/` |
| legacy merged μ+e overlay | `plotting/plots/merged/<disc>/` |
| legacy fiducial σ (incl + dσ/dη) | `plotting/plots/xsec/<disc>/` |

NB with the simfit's μ/e-shared `r`'s, the per-flavour observables are 100%
correlated with the comb ones (they differ only through MC template ratios) —
the comb plots are *the* result; per-flavour plots are meaningful as
independent measurements only from the legacy per-bin fits.

Every plot also carries the discriminant as a header line ("PF MET fit" /
"lep p_T fit" / "lep p_T (m_T>40) fit"), so a saved PNG self-identifies.
The disc→path mapping is single-sourced in `plotting/disc_variants.h`
(unknown tags are rejected, never silently misfiled). Pre-2026-08-03
*untagged* graph files (`charge_asym_fit_<chan>.root`) are accepted as a
met-only fallback when the tagged file is absent; old *flat* plot outputs
(`plots/charge_asym/*.png`, `plots/merged/*.png`, `plots/xsec/W_*.png`) are
stale leftovers — current outputs live in the per-disc subfolders.

Manual equivalents (what the driver runs; `disc` defaults to `"met"`
everywhere, shown explicit here — sub in `leppt`/`leppt_mt40` AND the matching
`pO_fit_out_leppt[_mt40]` tree for a variant):

```bash
cd analysis/
for c in mu ele; do
  F=../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/$c/summary/${c}_fitted_yields.root
  root -l -b -q "charge_asym.C+(\"$F\",\"../skim/rootfile/charge_asym_fit_${c}_met.root\")"
  root -l -b -q "FBratio.C+(\"$F\",\"../skim/rootfile/FBratio_fit_${c}_met.root\")"
done

cd ../plotting/
root -l -b -q 'observables.C+(false, "met")'   # muon     -> plots/{charge_asym,FBratio}/met/
root -l -b -q 'observables.C+(true, "met")'    # electron -> plots/Elec/{charge_asym,FBratio}/met/
# merged muon+electron overlay -> plots/merged/met/ :
root -l -q -b -e 'gROOT->LoadMacro("observables.C+"); observables_overlay("met");'
# fiducial cross sections -> plots/xsec/met/ :
root -l -b -q 'xsec_fiducial.C+("met")'        # W+/W-/W incl (mu+e), sigma_fid = N_fit / L
root -l -q -b -e 'gROOT->LoadMacro("xsec_fiducial.C+"); xsec_fiducial_diff(false, "met"); xsec_fiducial_diff(true, "met");'
```

`observables.C` overlays the data with **all four** nPDF theory bands
(EPPS21/nCTEQ15HQ/nNNPDF3.0/TUJU21nlo, drawn as filled bands only — no central
line / error bars). The merged overlay shows muon (black circles) + electron
(red squares) on the same axes; the sum-channel theory is weighted by the
combined μ+e fitted yields. Theory file optional (missing → data-only).
`xsec_fiducial.C` reads the summary CSVs (σ_fid = N_fit / L, L =
`pONorm::kLumi_invnb` = 46.5 nb⁻¹); `xsec_fiducial_diff` adds dσ_fid/dη per
channel (per-bin `<chan>_W_yields.csv`; η bins from `g_chargeAsym`).
NB this is the fiducial σ **before** the lepton-efficiency correction (= σ_fid×ε),
which is why μ and e differ — the gap is the channel efficiency, and they should
converge once ε is applied (universality). Stat (fit) uncertainty only.
For arbitrary regions, read the absolute yields directly from the CSVs.

## Module 6 — corrections & studies (`correction/`)

Run from `correction/`; outputs in `correction/plots/`, `correction/rootfile/`.

```bash
cd correction/
./run_qcd_abcd.sh [mu|ele|both]                    # ABCD QCD, logged  [also Module 3a]
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
