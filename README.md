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
            flavfit = the same fit per flavour → simfit_{mu,ele}/summary/)
     → fiducial_yields.C (r × σ_gen — EVERY observable is built from these,
            never from the raw fitted counts)
     → charge_asym.C / FBratio.C   → observables.C (comb = primary;
            μ-only vs e-only overlays from the flavfit trees, plots/flavfit/)
```

## TL;DR (full chain)

```bash
# ---- pO_analysis (plain ROOT) ----
cd correction  && ./run_trig_eff_mb.sh mu                                    # 0 muon trigger SF (the muon MC skim reads its rootfile since 2026-09-14)
cd ../skim     && ./run_all.sh all && ./run_lhe_updown.sh && ./run_ngen.sh   # 1,2 skims (muon MC weighted by ID x ISO x TRIG SFs) + nPDF Up/Down + N_gen
cd ../correction && ./run_qcd_abcd.sh                                        # 3a ABCD QCD (mu+ele, logged)
cd ../plotting && ./run_combine_inputs.sh                                    # 3b Combine inputs (mu+ele, logged)
#   theory graphs are one-time + yield-independent: root -l -q -b 'plotRpOtheory.C+'

# ---- fork (cmsenv) -- locally, or push to lxplus (see Module 4) ----
cd ../../HiggsAnalysis-CombinedLimit/test && cmsenv && ./run_pO_fits.sh --asimov   # 4 fit (simfit = DEFAULT; --asimov adds the closure fit)
#   ./run_pO_fits.sh both flavfit --asimov   # the mu-only + e-only simultaneous fits
#   ./run_pO_fits.sh both all --asimov       # grand + mu-only + e-only in one go

# ---- pO_analysis (plain ROOT) ----
cd ../../pO_analysis/analysis && ./run_observables.sh                            # 5 observables (met)
#   ./run_observables.sh leppt_mt40 | all    -> per-discriminant folders (Module 5; + mu-vs-e overlays when flavfit ran)
```

## Layout

| Directory          | Stage                                  | Entry point |
| ------------------ | -------------------------------------- | ----------- |
| `merge_rootfile/`  | (one-time) discover/hadd EOS ntuples   | `make_filelist.sh`, `hadd_from_list.sh` |
| `skim/`            | selection → per-sample histos; N_gen   | `run_all.sh`, `run_ngen.sh` |
| `correction/`      | ABCD QCD, isolation WP, Data/MC checks | `run_qcd_abcd.sh` → `qcd_abcd.C`, … |
| `plotting/`        | data/MC overlays + **Combine inputs**  | `mtandmet.C`, `dileptonpeak.C`, `plotRpOtheory.C` |
| `analysis/`        | r × σ_gen yields, charge asymmetry, F/B ratio | `fiducial_yields.C`, `charge_asym.C`, `FBratio.C` |
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
./run_lhe_updown.sh                     # 2.4 nPDF/qcdScale/alphaS Up/Down templates from the LHE weights (MC files; needs lhe_env.sh;
                                        #     members area-normalized before combining since 2026-09-14 -- see skim/lhe_updown.py)
```
CLI: `./run_all.sh <Zmm|Zee|Wmu|Wel|all> [samples…]` (samples ⊂
`Data DY Wp Wm DYtau Wptau Wmtau`).

**Check:** `skim/logs/*.log` (per-job OK/FAIL), `skim/output/*.txt` (cutflow),
`skim/rootfile/*_hist.root`. W files hold `h_{mt,met}_W{p,m}_y0..11` (+ `_FB`)
and the ABCD planes `h_iso_{met,mt}_{mu,ele}{Plus,Minus}`; Z files hold `hMass*`
+ Z kinematics + recoil histos.

**2.4 — LHE-weight templates (2026-09-07).** MC files additionally carry, for
every fit-template histogram (`h_met_*`, `h_leppt_mt40_*`, `hMass`), the TH2D
member twins `<h>_epps21` (107 = the LHAPDF set `EPPS21nlo_CT18Anlo_O16`),
`<h>_scale` (9) and `<h>_alphas` (5), filled by the skim with
`w·ttbar_w[i]/ttbar_w[0]` (layout: `skim/lhe_index.h`). `./run_lhe_updown.sh
[Wmu|Wel|Zmm|Zee|all]` then writes three Up/Down pairs per template into the
same files: `<h>_nPDFUp/Down` (LHAPDF's `PDFSet.uncertainty()` per bin),
`<h>_qcdScaleUp/Down` (μR/μF envelope over all 9 points, per-bin max/min) and
`<h>_alphaSUp/Down` (the α_s 0.119/0.117 member templates); log:
`skim/logs/lhe_updown_<chan>.log`; it must be re-run after every re-skim. Environment: `source skim/lhe_env.sh`
(LHAPDF 6.5.6 built under `~/local/lhapdf` for the PyROOT python; recipe in the
file). Re-skim regression: `./compare_reskim.sh <backup-dir>` checks every
pre-existing histogram for bit-identity (`compare_hists.C`).

**2.5 — muon efficiency SFs (2026-09-14).** The muon MC skims (`Wmu`, `Zmm`)
multiply the event weight by ID × ISO × trigger scale factors
(`skim/muon_sf.h`; data untouched): the Muon POG pp-2025 `TightID` and
`TightPFIso` SFs from the committed reduced JSON
`skim/sf/muon_sf_2025_TightID_PFIso_schemaV2.json` (regenerate from the POG
file with `python3 skim/sf/extract_muon_sf.py <file>`; the Z→μμ iso cut was
harmonized to the W's 0.15 the same day so ONE iso SF serves both) and the
MB-derived trigger SF, inclusive in rapidity, from
`correction/rootfile/trig_eff_mb_mu.root` (run `correction/run_trig_eff_mb.sh mu`
FIRST — a missing input is FATAL, never a silent SF = 1;
`PO_MUON_SF=off ./run_all.sh Wmu` disables it for checks). Each MC job logs an
`[SF]` block (inputs + provenance, the inclusive trigger SF and the per-y
consistency check, ⟨SF⟩ per source, Up/Down totals). The fit templates get
per-source twins `<h>_{muID,muIso,muTrig}Up/Down` (one source at ±1σ;
diagnostics) and the combined `<h>_muSFUp/Down` (the three added in
quadrature per bin), which is carried into the muon Combine inputs and cards
as ONE nuisance, exactly like nPDF/qcdScale/alphaS (+0.35/−0.38% on the W
templates, ±0.52% on the Z peak). Electron SFs: not yet.

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
[docs/AN_qcd_background.tex](docs/AN_qcd_background.tex), and ROOT prints them to stdout
only. `run_qcd_abcd.sh` pre-builds once (so the two channels cannot race on the
ACLiC artifacts), tees each channel to `correction/logs/qcd_abcd_<chan>.log`
(~290 lines) and echoes the transfer factors and the κ's to the terminal.
Takes ~2 s per channel, so regenerate freely — but never run the macro bare and
lose the report.

`combine_input_W.root` is **structured**: one TDirectory per fit region
(`Wp_lab_y0..11`, `Wm_lab_y*`, `Wp_fb_y*`, `Wm_fb_y*`, `Wp_incl`, `Wm_incl`,
`W_incl`), each with the 6 **absolute** templates `data_obs/signal/z/ztau/wtau/qcd`
(MET discriminant; per-y ABCD QCD). `combine_input_Z.root` has a `Z_incl/` dir
(`data_obs/signal/w/wtau/ztau`, mass peak). **Since 2026-09-07 every MC
process of every region also carries the LHE shape systematics
`<process>_{nPDF,qcdScale,alphaS}Up/Down`** (from the skim's Up/Down twins,
Module 2.4; their absence for plain `leppt` is what retired that variant on
2026-09-21 — see Module 4) **and, in the muon inputs,
the combined SF systematic `<process>_muSFUp/Down` (Module 2.5,
2026-09-14)**, and a sidecar
`<input>_systs.txt` next to each file lists what was written — the fork's
card generator reads it for its `shape` rows. Diagnostics of those templates
(per-region Up/Down-over-nominal plots, per-region overlays of the nominal
with all 106 EPPS21 member templates in `members/`, integral-shift summaries
vs rapidity, the inclusive-consistency tables and the LHAPDF-vs-Hessian
closure):
`./run_syst_shapes.sh [met|leppt_mt40|all]` → `plots[/Elec]/syst_shapes/<disc>/`
+ `logs/syst_shapes_<disc>.log`. **`combine_input_W_leppt_mt40.root`
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
./run_pO_fits.sh [mu|ele|both] [simfit|flavfit|all] [--disc met|leppt_mt40] [--dry-run] [--no-postfit] [--draw-only] [--asimov] [--no-statonly] [--no-contour] [--extract-only]
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
- **LHE shape systematics (2026-09-07, `LHE_SYST=auto|off|list`)** — three
  `shape` nuisances `nPDF` / `qcdScale` / `alphaS` from the
  `<process>_<syst>Up/Down` templates in the Combine inputs (Module 2.4 → 3),
  flag `1` on the MC processes, `-` on the data-driven `qcd` and the CR
  channels; one θ each for the whole card, group `lhe` (stat-only comparison:
  `--freezeNuisanceGroups lhe`). The generator reads the inputs' `_systs.txt`
  sidecars; `LHE_SYST=off` reproduces the earlier cards. Pulls + constraints
  land in `comb_summary.csv` as `<name>_theta` rows and enter the Asimov
  closure. With shape nuisances, `postfit_incl.C` needs the FitDiagnostics
  files (`sync_lxplus.sh download` now pulls them).
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
`extract_pO_simfit.C` stores it as `h_cov_yield[_FB]` (+ the POI-space
`h_cov_poi[_FB]`), `fiducial_yields.C` carries it over to the r × σ_gen yields,
and `charge_asym.C`/`FBratio.C` automatically include the cross terms when present.

Implementation: `make_pO_simfit_cards.sh` writes one 50-channel datacard + the
`multiSignalModel` map file per variant (plain files — `--dry-run` needs no
`cmsenv`); then `text2workspace.py -P ...:multiSignalModel` and `combine -M
FitDiagnostics --skipBOnlyFit` per variant. `--asimov` adds a prefit-Asimov
closure fit (`-t -1`): every fitted POI must come back at 1 — checked and
reported (PASS/FAIL) by the extraction. Outputs under
`pO_fit_out<suffix>/simfit/`: `summary/comb_W_yields.csv`, `comb_summary.csv`
(all POIs, fit status/covQual, covariance-propagated W⁺/W⁻/W inclusive sums,
Asimov closure rows) and **`comb_fitted_yields.root`** — the fit's record:
`h_yield_*` = r × the μ+e-summed prefit template integral (RECO-level counts)
plus the covariance matrices (incl. the 25×25 POI covariance `h_cov_poi*`).
Module 5 never builds an observable on those counts: `fiducial_yields.C` first
turns the r's into r × σ_gen (see Module 5). Postfit plots for every channel of the grand fit land in
`simfit/postfit/` (info box: that bin's `r_<C>_y<i>`, the global `r_Z` shown as
"DY norm", the channel's `qcd_norm`).

### The per-flavour fits — `flavfit` (2026-09-22)

`./run_pO_fits.sh both flavfit` runs the **same simultaneous fit once per lepton
flavour**: a μ-only likelihood (the 24 muon W channels + the μμ peak + the 6
muon ABCD control regions) and an e-only one, per binning variant. Everything
else is the grand fit's: 25 POIs (`r_<C>_y<i>` + `r_Z`, now measured by that
flavour alone — its own Z peak pins its own `r_Z`), every nuisance that acts
on the flavour (`lumi`, its two QCD rows, the LHE shapes, `muSF` in the muon
fit only), the three passes (nominal + `--statonly` = the stat error +
`--contour`), the `--asimov` closure and the same extraction. The card
generator writes the flavour halves of the grand card (`SIMFIT_FLAVS`;
verified identical per column, with Combine's own `multiSignalModel`), and the
grand card is byte-identical to before. Outputs `pO_fit_out<suffix>/simfit_mu/`
and `simfit_ele/` with the grand fit's layout and summary files
`simfit_<flav>_{W_yields.csv,summary.csv,fitted_yields.root}`; `mu flavfit`
runs one flavour; `both all` runs the grand fit and both flavour fits in one go.
Purpose: compare μ with e (Module 5 overlays them) — e.g. while the electron
SFs are not applied — which the grand fit cannot, since it forces one r on both.

**Removed 2026-09-22 (user decision): the legacy per-flavour per-bin pipeline**
(modes `perbin|incl|combined`: 48 separate two-channel cards per flavour, the
two-parameter `r` + `dy_norm` + free `qcd_norm` model with no systematics,
re-fitting the same Z data in every card; scripts `make_pO_datacards.sh`,
`extract_pO_yields.C`, `make_yields_from_csv.C`). Those modes now exit with a
pointer to `flavfit`; the scripts are in the fork's git history.

### W discriminant variants (`--disc`, 2026-07-30)

The W fit runs on two discriminants; `--disc` selects which:

| `--disc` | discriminant | selection | input file (per channel) | output tree | role |
|---|---|---|---|---|---|
| `met` | PF MET shape | plain W selection | `combine_input_W.root` | `pO_fit_out/` | backup |
| `leppt_mt40` | lepton pT | pT>25 && m_T>40 | `combine_input_W_leppt_mt40.root` | `pO_fit_out_leppt_mt40/` | **PRIMARY** |

Datacards and the fit model are IDENTICAL for both (same region names, same
POI scheme); only the input file, the output tree, the postfit x-title and the
QCD treatment change (`QCD_MODE` defaults to `abcd` for `leppt_mt40`, `lnN`
for `met` — see Module 4).

**RETIRED 2026-09-21 — a third variant `leppt`** (lepton pT, plain W selection,
`combine_input_W_leppt.root` → `pO_fit_out_leppt/`). It was dropped as a
discriminant on 2026-08-16; `pO_fit_out_leppt/` never had a simfit, and the
skim stores LHE/SF systematic twins only for `h_met_*` and `h_leppt_mt40_*`,
so its Combine input carried **no shape systematics** and could not be fitted
by the current card generator. `mtandmet.C` no longer writes that file or its
`plots[/Elec]/leppt/` stacks, and deletes any pre-retirement copy on each run;
`disc_variants.h` and `run_observables.sh` reject the tag with a message
naming the replacement. To bring it back, restore the two-entry variant table
at the top of `mtandmet.C` (`kVarNom`/`kVarMt40`) and the `leppt` rows in
those two consumers. NB `skim.C` still fills `h_leppt_W{p,m}_y*` (the
no-m_T-cut per-y histos) — they are now unread, but dropping them needs a
full re-skim, so they were left in place.

> **WARNING — carry the SAME `--disc` through the ENTIRE workflow of a variant
> run.** The out-trees carry no marker of which discriminant produced them —
> the coupling is only the directory suffix — so mixing steps corrupts or
> mislabels silently:
> - **fit**: `./run_pO_fits.sh both all --disc leppt_mt40`
> - **redraw**: `--draw-only` MUST repeat the same `--disc` (it selects the
>   out-tree AND the axis title; without it, lepton-pT plots get relabeled
>   "PF MET (GeV)" with no error). Same rule if you ever use `--out`.
> - **download**: `sync_lxplus.sh download` needs NO flag — it sweeps every
>   out-tree automatically, skipping absent ones.
> - **observables (Module 5)**: run `analysis/run_observables.sh <disc>` —
>   it carries the tag through every step automatically (reads the matching
>   `pO_fit_out<suffix>/` tree, writes disc-tagged graph files and per-disc
>   plot folders, stamps the discriminant on every plot). The histogram
>   names inside the yields files are identical across variants
>   (`h_yield_*`), so if you ever drive the macros by hand instead, the tree
>   you point at is the ONLY thing distinguishing a MET result from a pT
>   result — the driver exists so you never have to get that right manually.
> - **inputs**: `mtandmet.C` writes both files in one run; a variant whose
>   skim histograms are missing is skipped AND its stale file deleted, so a
>   missing `combine_input_W*.root` means "re-run Module 3", never "use the
>   old one".
>
> Physics note for the pT variants: with no low-MET region in the
> discriminant the QCD normalization is constrained by the in-fit ABCD control
> regions (`QCD_MODE=abcd`, the leppt_mt40 default) — check the fitted QCD
> multipliers and `qcd_rate_*` pulls in `<fit>_summary.csv` / on the postfit
> plots before trusting the composition.

Channel = `mu|ele|both` (the flavour(s) of `flavfit`; `simfit` is always μ+e);
mode = `simfit` (DEFAULT), `flavfit`, or `all` (grand + per-flavour).
`--dry-run` builds datacards without `cmsenv`; `--no-postfit` skips plots;
`--draw-only` redraws the postfit plots from an existing fit run (only `root`
needed — for cosmetic `draw_postfit_pO.C` changes; requires the `fits/` tree,
so redraw where the fits ran and `sync_lxplus.sh download --postfit`). A
failed pass is logged and skipped, not fatal — check
`pO_fit_out<suffix>/<fit>/fits/simfit_<B>/fit.log` and
`summary/extract_<fit>.log`. Each postfit plot carries on-plot fit-quality:
`χ²/ndf (p)` (Poisson GoF), the bin's `r`, `DY norm` (= `r_Z`), the QCD
parameter, and a red `status/covQ` flag if the fit didn't converge cleanly
(see the fork README "Diagnosing fit quality").

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
#   # the mu-only + e-only fits (flavfit), or all three fits in one go:
#   PO_PLOTS=... ./run_pO_fits.sh both flavfit --asimov
#   PO_PLOTS=... ./run_pO_fits.sh both all --asimov
#   # discriminant variants (SEE THE --disc WARNING above -- keep the flag
#   # consistent for every later step of that variant's workflow).
#   # NB --disc defaults to leppt_mt40, so the MET backup needs it explicitly:
#   PO_PLOTS=... ./run_pO_fits.sh --asimov --disc leppt_mt40
#   PO_PLOTS=... ./run_pO_fits.sh --asimov --disc met
./sync_lxplus.sh download                # <- lxplus, ALL out-trees (met + variants), each fit (simfit + simfit_{mu,ele})
./sync_lxplus.sh download --postfit      # also the postfit plots
```
lxplus paths (override via env): `ANA_LX=/afs/cern.ch/user/z/zheng/pO_analysis`,
`FORK_LX=/afs/cern.ch/user/z/zheng/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit`.

A downloaded fit can be re-extracted locally without re-fitting (only `root`):
`./run_pO_fits.sh [both simfit | both flavfit] --extract-only --disc <disc>`.

## Module 5 — final observables (`analysis/` + `plotting/`)

**One command per discriminant (2026-08-03):** `analysis/run_observables.sh`
runs the whole chain, carrying the disc tag through every filename and output
folder so the variants coexist without overwriting each other. It has two
conditional blocks, each run only when its fit outputs exist: the **PRIMARY
simfit chain** (`fiducial_yields.C` on the grand fit → r × σ_gen with the r
covariance, `charge_asym.C` + `FBratio.C` on THAT — errors include the fit
covariance — then `observables_comb`, the cross sections, the (σ_W, σ_Z)
contour and the inclusive postfit stacks), and since 2026-09-22 the **μ-vs-e chain** of the
per-flavour fits (`flavfit`, when `simfit_{mu,ele}/` exist), which overlays
the μ-only and e-only results — each with stat bars + syst boxes, exactly as
the grand fit's plots:

```bash
cd analysis/
./run_observables.sh              # met (default) -- the PF-MET-shape fit (backup)
./run_observables.sh leppt_mt40   # lepton-pT, mT>40 -- the PRIMARY fit
./run_observables.sh all          # every variant whose out-tree exists (others SKIP)
```

Outputs, per `<disc>` = `met` | `leppt_mt40`:

| output | path |
|---|---|
| **every fit**: fiducial yields (r × σ_gen + r covariance) | `skim/rootfile/fidyields_<fit>_<disc>.root` (`<fit>` = `comb`, `simfit_mu`, `simfit_ele`) |
| **every fit**: fiducial A_ch / R_FB graphs | `skim/rootfile/{charge_asym,FBratio}_fid_<fit>_<disc>.root` |
| **PRIMARY: simfit (μ+e comb)** plots | `plotting/plots/comb/{charge_asym,FBratio}/<disc>/` |
| **PRIMARY: simfit** fiducial σ: post-fit vs reco-MC vs **gen-MC** | `plotting/plots/comb/xsec/<disc>/` |
| **PRIMARY: simfit** y-inclusive postfit stacks (μ/e × W⁺/W⁻/W) | `plotting/plots/comb/postfit_incl/<disc>/` |
| **μ vs e** A_ch, R_FB overlays (+ `mu_vs_e_chi2.csv`) | `plotting/plots/flavfit/{charge_asym,FBratio}/<disc>/` |
| **μ vs e** σ (W⁺/W⁻/W, with the grand fit), dσ/dη per charge, (σ_W, σ_Z) contour overlay + each fit's own contour, `xsec_flavfit[_ratio].csv` | `plotting/plots/flavfit/xsec/<disc>/` |
| **μ vs e** each flavour fit's y-inclusive postfit stacks | `plotting/plots/flavfit/postfit_incl/<disc>/` |

For the generator-level overlay on the comb σ plots, produce the (one-time,
discriminant-independent) gen histograms first: `cd skim && root -l -b -q
'gen_xsec.C+'` → `skim/rootfile/gen_xsec.root` (missing file ⇒ the overlay is
skipped with a note, everything else unaffected).

**Every observable is built from r × σ_gen, never from raw counts** (user rule,
2026-09-22: there is no dedicated efficiency/acceptance correction, and
r × σ_gen-fid is what applies it, from MC). The σ's have been r × σ_gen since
2026-08-12; the charge asymmetry and R_FB were built from the count-based
fitted yields r × S until 2026-09-22, on the argument that A×ε cancels in a
ratio. That holds for A_ch (W⁺ and W⁻ share a bin; the switch moved it by
≤ 0.008) but NOT for R_FB, which divides two different |η_lab| regions
(η_lab = η_CM + 0.35): the electron ECAL crack sits in F at |η_CM| ≈ 1.2 and
in B at ≈ 1.9, so the count-based comb R_FB was off by ~20% there (1.026 →
0.848 and 0.908 → 1.087, 2.5–3 stat σ, on the 2026-09-21 fit). The count-based
r × S remain the fit's record (the CSVs, `*_fitted_yields.root`) and appear
only in the postfit stacks and the `xsec_fiducial_diag` A×ε diagnostic.

**Reading the μ-vs-e plots.** Every flavour result is quoted in ONE fiducial
volume (the pooled μ+e gen σ, bare lepton pT > 25 GeV, |η_lab| < 2.4), so
σ_e/σ_μ = r_e/r_μ exactly. The A_ch and R_FB overlays use the same
**acceptance-corrected** yields r × σ_gen (`analysis/fiducial_yields.C`); with
the count-based r × S each flavour's own A×ε would enter R_FB, and on
identical physics the μ and e R_FB would differ by up to 60%. The two fits' statistical errors are independent; the
systematic boxes are largely COMMON (lumi 3% dominates σ and moves both
flavours together; it cancels in the ratios), so the printed σ_e/σ_μ carries
the stat error and the per-bin χ² is quoted against stat and total errors.
(With the grand fit's μ/e-shared `r`'s the flavours cannot differ at all —
the comb plots are the result; the flavfit plots are the consistency check.)

Every plot also carries the discriminant as a header line ("PF MET fit" /
"lep p_T (m_T>40) fit"), so a saved PNG self-identifies. The disc→path
mapping is single-sourced in `plotting/disc_variants.h` (unknown tags are
rejected, never silently misfiled). The graph files are
`<stem>_fid_<fit>_<disc>.root` since 2026-09-22 (`pODisc::GraphFile`, no
fallback): the count-based `*_fit_*` ones — incl. the pre-2026-08-03
*untagged* `charge_asym_fit_<chan>.root`, once a met-only fallback — are no
longer read. The matching *untagged* plot outputs
(`plots/charge_asym/chargeAsym_mt.png`, `plots/FBratio/RFB_mt_*.png`,
`plots/merged/*_overlay.png`, `plots/xsec/W_*.png`, + the `Elec` twins —
30 files) were **deleted 2026-09-21**: no live code path could overwrite
them (every writer targets `<disc>/`), and their `_mt` names came from the
discriminant naming retired in the 2026-07-30 audit, so they advertised a
quantity the fit never used. Current outputs live only in the per-disc
subfolders.

Manual equivalents (what the driver runs, for `leppt_mt40`; the loop is the
fiducial step of BOTH chains, the `comb` pass is also what `observables_comb`
reads; `S=../../HiggsAnalysis-CombinedLimit/test/pO_fit_out_leppt_mt40`):

```bash
cd analysis/
for T in simfit_mu simfit_ele comb; do
  D=$T; [ $T = comb ] && D=simfit
  root -l -b -q "fiducial_yields.C+(\"$S/$D/summary/${T}_W_yields.csv\",\"$S/$D/summary/${T}_fitted_yields.root\",\"../skim/rootfile/fidyields_${T}_leppt_mt40.root\")"
  root -l -b -q "charge_asym.C+(\"../skim/rootfile/fidyields_${T}_leppt_mt40.root\",\"../skim/rootfile/charge_asym_fid_${T}_leppt_mt40.root\")"
  root -l -b -q "FBratio.C+(\"../skim/rootfile/fidyields_${T}_leppt_mt40.root\",\"../skim/rootfile/FBratio_fid_${T}_leppt_mt40.root\")"
done
cd ../plotting/
root -l -b -q -e 'gROOT->LoadMacro("observables.C+");   observables_flav("leppt_mt40");'        # A_ch, R_FB overlays
root -l -b -q -e 'gROOT->LoadMacro("xsec_fiducial.C+");  xsec_fiducial_flav("leppt_mt40");'     # sigma + dsigma/deta overlays
root -l -b -q -e 'gROOT->LoadMacro("xsec_contour.C+");   xsec_contour_WZ_flav("leppt_mt40","lab");'           # contour overlay
root -l -b -q -e 'gROOT->LoadMacro("xsec_contour.C+");   xsec_contour_WZ_fit("leppt_mt40","lab","simfit_mu");' # one fit's own contour
root -l -b -q -e 'gROOT->LoadMacro("postfit_incl.C+");   postfit_incl_fit("leppt_mt40","simfit_mu");'
```

`observables.C` overlays the data with **all four** nPDF theory bands
(EPPS21/nCTEQ15HQ/nNNPDF3.0/TUJU21nlo, drawn as filled bands only — no central
line / error bars); theory file optional (missing → data-only). The per-fit
series use the conventions of `plotting/fit_variants.h` (μ blue circle, e red
square, combined black diamond). **Removed 2026-09-22 with the legacy
per-flavour fits:** `observables(isElec, disc)`, `observables_overlay`, the
N_fit/L `xsec_fiducial(disc, muCsv, eleCsv)` and `xsec_fiducial_diff`
(`observables(disc)` / `xsec_fiducial(disc)` now run the comb + μ-vs-e views).
Their last outputs (`plots[/Elec]/{charge_asym,FBratio}/`, `plots/merged/`,
`plots/xsec/`, `skim/rootfile/{charge_asym,FBratio}_fit_{mu,ele}*.root`) are
orphaned — nothing rewrites them — as are the count-based comb graphs
`skim/rootfile/{charge_asym,FBratio}_fit_comb_<disc>.root` superseded by the
`_fid_` ones on 2026-09-22.

## Module 6 — corrections & studies (`correction/`)

Run from `correction/`; outputs in `correction/plots/`, `correction/rootfile/`.

```bash
cd correction/
./run_qcd_abcd.sh [mu|ele|both]                    # ABCD QCD, logged  [also Module 3a]
./run_trig_eff_mb.sh [mu|ele|both]                 # trigger turn-on, data vs W MC, MB-triggered denominator, logged
./run_charge_flip.sh [mu|ele|both]                 # lepton charge-flip rate from W MC, logged
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
  `AsymErr`/`RatioErr` propagate σ²). The fiducial-yield histos
  (`fidyields_*`) carry the fit error (r covariance × σ_gen) in Sumw2 plus the
  covariance matrices, so charge_asym/FBratio errors are the fit uncertainties.
- **Corrections WIP.** Recoil, lepton SFs, momentum scale/smearing not all
  applied — don't assume MC in `skim/rootfile/` is fully corrected.
- **Acceptance / efficiency.** The fitted yields r × S are RECONSTRUCTED counts
  (S carries the MC A×ε), so no observable is built on them: every σ, A_ch and
  R_FB uses r × σ_gen-fid (bare lepton pT > 25 GeV, |η_lab| < 2.4), which applies
  the MC A×ε bin by bin — there is no dedicated efficiency/acceptance
  correction (user rule 2026-09-22). A×ε nearly cancels in A_ch but NOT in
  R_FB (F and B are different detector regions). Going beyond the fiducial
  volume (total σ, boson-level theory) would additionally need the acceptance.
- **Pre-existing inter-channel asymmetries** (intentional; see CLAUDE.md):
  DY-veto pT 15 (μ) vs 10 (e) GeV; isolation 0.15 (μ) vs 0.095 (e) — both
  re-confirmed optimal under the Δβ-corrected relIso (2026-07-06, MuonPOG
  convention, all channels; needs re-skim). All electron ID gates are
  `eleMVAIdWP95` since 2026-07-02 (`skim_Zee` included).
- **No tests/CI.** Validate by re-running + inspecting logs/plots. Commit
  messages are often "xx" — read diffs, not `git log`.
