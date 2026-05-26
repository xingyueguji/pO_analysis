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

Isolation studies (ROC curves, working-point tuning) are still per-flavour and were not folded into `skim.C`:
- [skim/isolation.C](skim/isolation.C) — muon
- [skim/isolation_ele.C](skim/isolation_ele.C) — electron

**Selection (8-step cutflow in the W macros)**:
1. ≥1 PF lepton with pT > 25 GeV
2. pO event-quality filters (auto-detected from branches)
3. Trigger: `HLT_PAL3Mu12` (muon) or `HLT_PAL3Ele12` (electron)
4. DY veto: reject events with an OS dilepton pair (pT>15, Tight, relIso<0.15, m>30 GeV)
5. ≥1 Tight ID lepton
6. Leading lepton pT > 25 GeV
7. Leading lepton PF relIso < 0.15
8. Trigger matching on leading lepton

**Output**: `skim/rootfile/{channel}_pO_PFMet_{sample}_hist.root` containing
12 rapidity-binned MT histograms (`h_mt_Wp_y0..y11`, `h_mt_Wm_y0..y11`),
MET distributions, isolation histos, dilepton mass for Z channels.

Per-job logs in `skim/logs/`, cutflow text in `skim/output/`.

### Stage 2 — `plotting/`

Data/MC overlay, background estimation, plots.

- [plotting/mtandmet.C](plotting/mtandmet.C) — MT and MET data/MC plots
- [plotting/mtandmet_overlay.C](plotting/mtandmet_overlay.C) — MC stack overlays
- [plotting/qcd_sideband_fit_and_extrapolate.C](plotting/qcd_sideband_fit_and_extrapolate.C) — QCD shape from low-MT / anti-iso sideband, Rayleigh-like fit, extrapolate to signal region
- [plotting/PlotsIsoROC.C](plotting/PlotsIsoROC.C), [plotting/PlotIsoROC_ele.C](plotting/PlotIsoROC_ele.C) — isolation ROC curves
- [plotting/dileptonpeak.C](plotting/dileptonpeak.C) — Z peak plots
- [plotting/observables.C](plotting/observables.C), [plotting/plotZcurve.C](plotting/plotZcurve.C), [plotting/plotRpOtheory.C](plotting/plotRpOtheory.C) — final-observable plots, theory comparisons
- [plotting/CMS_lumi.C](plotting/CMS_lumi.C), [plotting/plotting_helper.C](plotting/plotting_helper.C) — style / helpers

### Stage 3 — `analysis/`

Final observable extraction from the rapidity-binned histograms:

- [analysis/analysis_helpers.h](analysis/analysis_helpers.h) — shared header. `pOAnalysis::YieldInRange`, `AsymErr`, `RatioErr`, and the `kPORapidityShift = 0.3466` constant (single source of truth for the pO→CM boost).
- [analysis/charge_asym.C](analysis/charge_asym.C) — `A = (N+ − N−) / (N+ + N−)` vs rapidity (12 bins, y ∈ [−2.4, 2.4]).
- [analysis/FBratio.C](analysis/FBratio.C) — forward/backward ratio in |y_CM|. Lab-frame `yEdges` chosen symmetric around Δy so the CM-frame bins are symmetric around 0 (required for F/B pairing).

### `merge_rootfile/` — input prep

Utilities to scan EOS and `hadd` per-sample ROOT files into the consolidated
`pO_2025.root` that `skim/run_all.sh` reads:
- [merge_rootfile/make_filelist.sh](merge_rootfile/make_filelist.sh) — discover files on EOS
- [merge_rootfile/hadd_from_list.sh](merge_rootfile/hadd_from_list.sh) — merge them
- `MC_*.txt`, `DATA_pass_*.txt`, `version_*.txt` — per-sample filelists

## Input data

Centralized file on EOS, hardcoded in `skim/run_all.sh`:

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
  fully corrected when reading skim output.

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
- **`analysis/charge_asym.C`, `analysis/FBratio.C`** — only call
  `Integral()` on already-Sumw2'd inputs, so the *yields* are fine. But
  the error formulas in `analysis/analysis_helpers.h` (`AsymErr`,
  `RatioErr`) are hand-derived Poisson on raw counts. **Today this is
  correct** because skim Fill calls are unweighted. **But the moment any
  event weight enters the skim Fill calls** (lepton SFs, recoil correction,
  pileup weight, etc.), those analysis errors will be silently wrong —
  they don't look at the histogram's stored σ². The defensive fix is to
  switch `YieldInRange` to return both integral and error via
  `TH1::IntegralAndError`, and have `AsymErr` / `RatioErr` take errors
  as inputs instead of recomputing them. Not yet done — flag this when
  the user adds the first event weight.

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
