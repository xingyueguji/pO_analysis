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
- [skim/legacy/](skim/legacy/) — the original four per-channel macros, preserved verbatim for diffing/fallback. The unified files are intended to be byte-identical in physics output to these.
- [skim/count_ngen.C](skim/count_ngen.C) + [skim/run_ngen.sh](skim/run_ngen.sh) — compute N_gen (Σ gen weight over **all** events, no selection) per MC sample → `skim/rootfile/ngen.root` + `skim/output/ngen.txt`. The cross-section-normalization denominator. See "MC normalization" below.
- [skim/mc_norm.h](skim/mc_norm.h) — single source of truth for the absolute MC→data scale `k_s = σ·L/N_gen` (`pONorm::MCScale`); consumes `ngen.root`. See "MC normalization" below.

Isolation studies (ROC curves, working-point tuning) are per-flavour and live in
`correction/` (moved out of `skim/`):
- [correction/isolation.C](correction/isolation.C) — muon
- [correction/isolation_ele.C](correction/isolation_ele.C) — electron
- plotted by [correction/PlotsIsoROC.C](correction/PlotsIsoROC.C) / [correction/PlotIsoROC_ele.C](correction/PlotIsoROC_ele.C). See the `correction/` section below.

**Selection (8-step cutflow in the W macros)**:
1. ≥1 PF lepton with pT > 25 GeV
2. pO event-quality filters (auto-detected from branches)
3. Trigger: `HLT_PAL3Mu12` (muon) or `HLT_PAL3Ele12` (electron)
4. DY veto: reject events with an OS dilepton pair (pT>15, Tight, relIso<0.15, m>30 GeV)
5. ≥1 Tight ID lepton
6. Leading lepton pT > 25 GeV
7. Leading lepton PF relIso < 0.15
8. Trigger matching on leading lepton

**Output**: per-sample ROOT files in `skim/rootfile/` (W: `{channel}_pO_PFMet_{sample}_hist.root`;
Z: `ZToMuMu_pO2025_*` / `ZToEE_pO2025_*`). W files contain 12 rapidity-binned MT
histograms (`h_mt_Wp_y0..y11`, `h_mt_Wm_y0..y11`), MET distributions, isolation
histos. **Both `skim_Wmu` and `skim_Wel`** additionally store the ABCD inputs for
the QCD/low-MET background: per-charge 2D histos `h_iso_met_{mu,ele}{Plus,Minus}`
(relIso × PF MET) and `h_iso_mt_{mu,ele}{Plus,Minus}` (relIso × m_T), filled for
every event passing the full W selection **except** the isolation cut, so the
relIso axis spans iso-pass (signal) and the anti-iso sideband. Binning: relIso
0–1.0 in 0.01 steps; MET 0–120/2 GeV; m_T 0–200/2.5 GeV. Regions are chosen
downstream by projecting these (no re-skim) in `correction/qcd_abcd.C`. Z files contain dilepton mass (`hMass`, `hMass_extended`, `hMass_vipul`)
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

- [plotting/mtandmet.C](plotting/mtandmet.C) — MT and MET data/MC plots. **MET stacks (muon AND electron) include the data-driven ABCD QCD** from `correction/rootfile/qcd_abcd_{mu,ele}.root` (run `correction/qcd_abcd.C[+(true)]` first): added to the inclusive MET stack (rigorous) and split across the per-y bins by the per-bin low-MET excess (data−EWK at MET<30) proxy, keeping the inclusive template shape (per-y sums back to inclusive). Also written as `qcd` into `combine_input_inclusive.root`. Channel auto-selected by `isElec`. (For electrons QCD is the *dominant* low-MET component — see `qcd_abcd.C`.)
- [plotting/mtandmet_overlay.C](plotting/mtandmet_overlay.C) — MC stack overlays
- [plotting/dileptonpeak.C](plotting/dileptonpeak.C) — Z peak plots
- (QCD sideband fit and isolation ROC curves moved to `correction/` — see below)
- [plotting/observables.C](plotting/observables.C), [plotting/plotZcurve.C](plotting/plotZcurve.C), [plotting/plotRpOtheory.C](plotting/plotRpOtheory.C) — final-observable plots, theory comparisons
- [plotting/CMS_lumi.C](plotting/CMS_lumi.C), [plotting/plotting_helper.C](plotting/plotting_helper.C) — style / helpers

### Stage 3 — `analysis/`

Final observable extraction from the rapidity-binned histograms:

- [analysis/analysis_helpers.h](analysis/analysis_helpers.h) — shared header. `pOAnalysis::YieldInRange`, `AsymErr`, `RatioErr`, and the `kPORapidityShift = 0.3466` constant (single source of truth for the pO→CM boost).
- [analysis/charge_asym.C](analysis/charge_asym.C) — `A = (N+ − N−) / (N+ + N−)` vs rapidity (12 bins, y ∈ [−2.4, 2.4]).
- [analysis/FBratio.C](analysis/FBratio.C) — forward/backward ratio in |y_CM|. Lab-frame `yEdges` chosen symmetric around Δy so the CM-frame bins are symmetric around 0 (required for F/B pairing).

### `correction/` — corrections & studies

Orthogonal to the main skim→plotting→analysis flow: code that derives/validates
corrections and background estimates. **All macros are run from `correction/`**
and write outputs there (`correction/plots/`, `correction/rootfile/`). The two
that need `plotting_helper.C` `#include "../plotting/plotting_helper.C"`; the
isolation study is self-contained (writes/reads `correction/rootfile/`); the QCD
fit reads the *main* W skim output at `../skim/rootfile/`.

- [correction/qcd_abcd.C](correction/qcd_abcd.C) — **current** data-driven QCD/low-MET background via the ABCD method, **both muon and electron** (`qcd_abcd.C+` = muon, `qcd_abcd.C+(true)` = electron). Plane = relIso × PF MET (also stores relIso × m_T). Regions: iso-pass (`relIso < isoCut`) / iso-fail (anti-iso window `[isoFailLo, isoFailHi)`) × MET high/low (`metCut`); all four counts are QCD-only (`data − Σ EWK·k_s`, absolute MC from `mc_norm.h`). Reports the transfer factor `T = B/D`, the signal-region QCD `A = B·C/D`, and writes the iso-pass QCD MET/m_T **template** (anti-iso shape × T) to `correction/rootfile/qcd_abcd_{mu,ele}.root` for Combine, plus closure overlays (data vs EWK ± ABCD QCD). Channel diffs: electron uses `isoCut=0.095`, anti-iso `[0.20,1.0)` (muon `0.15`, `[0.30,1.0)`) and the `Wp_ele/Wm_ele/DYee` MCScale labels. Boundaries are projections → retune `ABCDConfig` with no re-skim. **Result:** muon QCD is small (T≈0.17; ≈3% of the high-MET signal region). Electron QCD was originally huge with the loose `eleCutIdWP95` (T≈0.30; ≈24–26%), which **motivated the ID switch to `eleMVAIdWP95`** (see the isolation/ID study + FIXMEs); after the switch the electron QCD dropped to **T≈0.33–0.38, ≈10–11%** of the high-MET signal region (iso-pass QCD ≈10037→3634, W signal only −6.5%, S/B 0.34→0.87). Closure is good in both channels.
- [correction/qcd_sideband_fit_and_extrapolate.C](correction/qcd_sideband_fit_and_extrapolate.C) — **superseded** by `qcd_abcd.C`. QCD shape from anti-iso sideband, Rayleigh-like fit + linear shape-parameter extrapolation to signal iso; did not behave at pO statistics. Kept for reference.
- [correction/isolation.C](correction/isolation.C) (muon), [correction/isolation_ele.C](correction/isolation_ele.C) (electron) — isolation/ID ROC study. Signal = OS Z-window pairs; backgrounds = SS pairs (≈empty for μ) and single-lepton MET<5 (QCD-like). **Pass the current data file** (default path is stale EOS) and run with ACLiC: `root -l -q 'isolation.C+("/Users/zhenghuang/pO_2026_May_26/Data_May_26.root", false)'` (false=tight muon ID; needed `TParameter.h` include). Muon scans continuous relIso × 3 cones × 3 iso-defs (+ a `ggbranch` variant = the skim's `RelIsoPF`); **electron `isolation_ele.C` scans the 8 ID WPs × 4 MVA-Iso WPs AND now a continuous-relIso scan per ID** (added 2026-06-24 — the discrete MVA-Iso-WP path did NOT match the skim's continuous relIso cut; it also prints a per-ID optimality table). **Conclusions:** μ relIso<0.15 is the Youden-J(QCD) optimum (ΔβPU adds ~nothing); e the *ID* was the problem (`CutWP95` worst, AUC_MET 0.870 vs 0.911 MVA) → switched `skim_Wel` to `eleMVAIdWP95`, relIso<0.095 kept.
- [correction/PlotsIsoROC.C](correction/PlotsIsoROC.C) (`PlotsIsoROC.C+(false)` for the tight set), [correction/PlotIsoROC_ele.C](correction/PlotIsoROC_ele.C) — isolation ROC curves → `correction/plotsROC[_ele]/`; plus `plotsROC_ele/roc_MET_continuous_allID.png` (one-off overlay of the 8 IDs' continuous-relIso QCD-ROC).
- [correction/plot_iso_summary.C](correction/plot_iso_summary.C) — **per-fixed-ID summary in one consistent cosmetic**: a 2-pad canvas (LEFT = continuous efficiency vs relIso cut for signal/QCD/SS with the operating cut marked + ε_sig/ε_QCD; RIGHT = the corresponding ROC + AUC with the operating-point star). Electron → one per ID (`plotsROC_ele/summary_<ID>.png`, reads the continuous scan from `isolation_ele.C`); muon → single tight-ID summary (`plotsROC/summary_muon.png`, branch relIso = skim RelIsoPF). `root -l -q 'plot_iso_summary.C+'`.
- [correction/dataMC_kinematics.C](correction/dataMC_kinematics.C) — Data vs signal-MC (DY) overlay + ratio pad for the Z kinematics (`h_Zpt/h_Zeta/h_Zy/h_Zphi`, `h_lepPt/h_lepEta/h_lepPhi`, `hMass`), shape-normalized. Used for the Data/MC boson-pT check. `root -l -q 'dataMC_kinematics.C+("Zmm")'`.
- [correction/recoil_raw.C](correction/recoil_raw.C) — raw look at the Z→μμ hadronic recoil (`u_par`/`u_perp`) for the MET recoil correction: data vs DY MC, inclusive and in q_T slices (projected from the 2D histos), **plus a printed entry-count table per q_T slice** to judge statistics/binning before fitting. Slice edges = editable `kQtEdges` array (projections → no re-skim to change). NO fit yet; the planned `recoil_fit.C` (double-Gaussian per q_T bin → μ(q_T)/σ(q_T), AN Sec 6) is the next step. `root -l -q 'recoil_raw.C+'`.
- [correction/recoil_raw_ele.C](correction/recoil_raw_ele.C) — electron-channel twin of `recoil_raw.C`: same raw-recoil look for Z→ee (reads `ZToEE_pO2025_*`, writes `plots/recoil_Zee/`). `root -l -q 'recoil_raw_ele.C+'`.

The shared ratio-pad helper `SaveDataMCRatio(...)` lives in [plotting/plotting_helper.C](plotting/plotting_helper.C).

### `merge_rootfile/` — input prep

Utilities to scan EOS and `hadd` raw ntuple files into the **per-sample** ROOT
files that `skim/run_all.sh` reads (now the May-26 production — see "Input data";
earlier it was a single consolidated `pO_2025.root`):
- [merge_rootfile/make_filelist.sh](merge_rootfile/make_filelist.sh) — discover files on EOS
- [merge_rootfile/hadd_from_list.sh](merge_rootfile/hadd_from_list.sh) — merge them
- `MC_*.txt`, `DATA_pass_*.txt`, `version_*.txt` — per-sample filelists

## Input data

**Per-sample** ROOT files (the May-26 production), not a single consolidated file.
Single source of truth: `kDefaultDataFile` (data) and `ResolveMCSample` (MC) in
`skim/skim_common.h` (`run_all.sh` greps the data path; override per-run via the
`DATA_FILE` env var).

Currently pointing at the **local** copy under `~/pO_2026_May_26/` — written as
absolute paths (`/Users/zhenghuang/pO_2026_May_26/…`) because ROOT's
`TFile::Open` doesn't reliably expand `~` — switched from EOS 2026-06-23. To read
off EOS/lxplus again, repoint `kDefaultDataFile` + `inputBase` (in
`ResolveMCSample`) to
`root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2026_May_26/`.
Files: `Data_May_26.root` + 9 MC (`MC_DY{mu,ee}_May26.root`,
`MC_W{p,m}_{mu,ele}_May26.root`, `MC_DYtau_May26.root`, `MC_W{p,m}_tau_May26.root`).
Three low-mass DY files (`MC_DY{mu,ee,tau}_low_May26.root`) are present but
**unwired** (skipped — see "MC normalization").

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

Files on `zheng/po-analysis`:
- `test/run_fit.sh` — driver
- `test/my_script/make_combine_input.C`, `make_combine_input_Z.C` — read this repo's
  histograms from `/afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots/`, produce
  Combine inputs (`data_obs`, signal, backgrounds)
- `test/my_script/testdatacard_inclusive.txt` — W (inclusive) datacard
- `test/my_script/testdatacard_Zmumu.txt`, `testdatacard_Zee.txt` — Z datacards
- `test/my_script/draw_postfit_{inclusive,Zmumu,Zee}.C` — postfit plots

Each channel is fit separately. Processes are signal + W/W-tau/Z-tau backgrounds.
`rateParam`s (range 0–20) on W/Z/tau rates let the fit float them — the POIs
are effectively the cross-sections. Shape-only datacards (no analytic templates).
Fit method: `combine -M FitDiagnostics`.

## Current state of the work

Active (in flight, recent commits):
- **Electron isolation + background** — being added in `DrawWToElecNu_PFMet.C` /
  `mtandmet.C` / `isolation_ele.C`. Working point not finalized; PU correction
  for electron isolation is TODO.
- **Tau channels as background** — sample enums (`kWptau`, `kWmtau`, `kDYtau`)
  and `Wptau/Wmtau/DYtau` filelists exist; selection code largely there but
  validation pending. Commit `447a2b5` ("before adding tau") marks the boundary.
- **Corrections still WIP** — recoil corrections, lepton scale factors,
  momentum scale/smearing are not all applied yet. Don't assume MC is
  fully corrected when reading skim output. (The per-event *generator*
  weight IS now applied — see "Generator event weight" below — but these
  reco-level corrections are separate and still missing.)
- **ABCD QCD background (muon + electron, done 2026-06-23)** — replaces the
  Rayleigh shape-extrapolation. Both `skim_Wmu`/`skim_Wel` store the relIso × MET
  (and × m_T) 2D per charge (`h_iso_{met,mt}_{mu,ele}{Plus,Minus}`);
  `correction/qcd_abcd.C[+(true)]` does `N_A = N_B·N_C/N_D` on QCD-only
  (EWK-subtracted) counts and emits the iso-pass QCD template, wired into the
  `mtandmet.C` MET stacks for both channels. Electron QCD is ~10× the muon and
  is the dominant low-MET background there. Remaining: use these templates in
  the Combine fork datacards.

Future plans / on the TODO list (not yet started):
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
  cut the electron ABCD QCD ~64% (≈10037→3634 iso-pass; signal −6.5%). **`skim_Zee`
  still uses `eleCutIdWP95`** (`skim.C:1598,1899` comments) — apply the same swap there
  if desired (not yet done).
- `isolation_ele.C` / `skim_Wel` — electron relIso still lacks the pileup (ΔβPU)
  correction; re-confirm the 0.095 optimum once it's added (muon study showed PU gives
  ~no AUC gain, electrons untested).

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
  `analysis/FBratio.C` read by name (`h_mt_Wp_y{0..11}`, `h_mt_Wm_y{0..11}`),
  and so do the Combine input scripts. Rename in skim → must rename in
  analysis AND in the Combine fork.
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
overshoot = uncalibrated **MET/recoil** (NOT efficiency). The Z mass peak (no MET)
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
- **Isolation cut**: 0.15 in `skim_Wmu`, 0.095 in `skim_Wel` (with the FIXME
  about PU correction) — **intentional, still being tuned**.
- **DY-veto gate type**: continuous PF relIso (muon) vs integer
  `eleCutIdWP95` / `eleMVAIsoWP95` (electron) — pre-existing.
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
