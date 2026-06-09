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
histos. Z files contain dilepton mass (`hMass`, `hMass_extended`, `hMass_vipul`)
plus kinematics of the dilepton system and its leptons (`h_Zpt/h_Zeta/h_Zy/h_Zphi`,
`h_lepPt/h_lepEta/h_lepPhi` — both legs; `h_Zy` is lab-frame rapidity), filled for
iso-selected OS pairs in the Z peak [60,120] GeV (consumed by
`correction/dataMC_kinematics.C`). **`skim_Zmm` additionally** stores the hadronic
recoil for the MET correction (Sec 6 of AN2017_058): `h_uPar`/`h_uPerp` (1D
inclusive) and `h_uPar_qT`/`h_uPerp_qT` (2D, recoil component vs q_T), where
`u = −MET − q_T` and `q_T` is the dimuon pT. u-axis = AN's 2 GeV/c binning; q_T
axis fine so the fit binning can be chosen later by projecting (no re-skim). PF
tree is **soft-gated** (loud `[WARN]`, empty recoil histos if `pftree` absent).
Electron recoil (`skim_Zee`) not added yet.

Per-job logs in `skim/logs/`, cutflow text in `skim/output/`.

### Stage 2 — `plotting/`

Data/MC overlay, background estimation, plots.

- [plotting/mtandmet.C](plotting/mtandmet.C) — MT and MET data/MC plots
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

- [correction/qcd_sideband_fit_and_extrapolate.C](correction/qcd_sideband_fit_and_extrapolate.C) — QCD shape from anti-iso sideband, Rayleigh-like fit (moved from `plotting/`).
- [correction/isolation.C](correction/isolation.C), [correction/isolation_ele.C](correction/isolation_ele.C) — isolation working-point study inputs (moved from `skim/`).
- [correction/PlotsIsoROC.C](correction/PlotsIsoROC.C), [correction/PlotIsoROC_ele.C](correction/PlotIsoROC_ele.C) — isolation ROC curves (moved from `plotting/`).
- [correction/dataMC_kinematics.C](correction/dataMC_kinematics.C) — Data vs signal-MC (DY) overlay + ratio pad for the Z kinematics (`h_Zpt/h_Zeta/h_Zy/h_Zphi`, `h_lepPt/h_lepEta/h_lepPhi`, `hMass`), shape-normalized. Used for the Data/MC boson-pT check. `root -l -q 'dataMC_kinematics.C+("Zmm")'`.
- [correction/recoil_raw.C](correction/recoil_raw.C) — raw look at the Z→μμ hadronic recoil (`u_par`/`u_perp`) for the MET recoil correction: data vs DY MC, inclusive and in q_T slices (projected from the 2D histos), **plus a printed entry-count table per q_T slice** to judge statistics/binning before fitting. Slice edges = editable `kQtEdges` array (projections → no re-skim to change). NO fit yet; the planned `recoil_fit.C` (double-Gaussian per q_T bin → μ(q_T)/σ(q_T), AN Sec 6) is the next step. `root -l -q 'recoil_raw.C+'`.

The shared ratio-pad helper `SaveDataMCRatio(...)` lives in [plotting/plotting_helper.C](plotting/plotting_helper.C).

### `merge_rootfile/` — input prep

Utilities to scan EOS and `hadd` per-sample ROOT files into the consolidated
`pO_2025.root` that `skim/run_all.sh` reads:
- [merge_rootfile/make_filelist.sh](merge_rootfile/make_filelist.sh) — discover files on EOS
- [merge_rootfile/hadd_from_list.sh](merge_rootfile/hadd_from_list.sh) — merge them
- `MC_*.txt`, `DATA_pass_*.txt`, `version_*.txt` — per-sample filelists

## Input data

Centralized file on EOS. Single source of truth: `kDefaultDataFile` in
`skim/skim_common.h` (`run_all.sh` greps it; override per-run via the `DATA_FILE`
env var). MC file paths live in the same header (`ResolveMCSample`).

```
root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root
```

Tree: `ggHiNtuplizer/EventTree` (~2.2M events). Also contains `hltanalysis/HltTree`
and PF candidate / muon / electron collections.

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

Future plans / on the TODO list (not yet started):
- **ABCD method for QCD background** — the current
  `qcd_sideband_fit_and_extrapolate.C` (Rayleigh-like shape from low-MT /
  anti-iso sideband) doesn't work well at the available statistics. Plan is
  to replace with an ABCD estimate using two ~uncorrelated handles
  (typically iso × MT): `N_A = N_B · N_C / N_D`.
- **Tau cross-section measurement** — currently τ is only a background to W/Z.
  Adding W→τν / Z→ττ as their own signal channels is on the wishlist (tagged
  as "probably just an attempt").

Known FIXMEs (mostly in the electron path):
- `DrawDielectronPeak.C` — ECAL coverage comment says 2.4 but should be 2.5
  with the 1.4442–1.566 transition gap excluded
- `DrawDielectronPeak.C` — "WP95 is loose" — tight-ID cut needs revisiting
- `isolation_ele.C` — electron isolation lacks pileup correction
- `DrawWToElecNu_PFMet.C` — "Need to redo electron isolation study, not using
  working point"

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
