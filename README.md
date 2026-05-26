# pO_analysis

W and Z boson cross-section measurement in proton–oxygen (pO) collisions.

Four physics channels — W→μν, W→eν, Z→μμ, Z→ee — each producing cross-sections,
charge asymmetries vs rapidity, and forward/backward ratios in |y| (CM frame).
Tau channels (W→τν, Z→ττ) are treated as backgrounds; a possible standalone
tau measurement is on the wishlist.

The analysis is a three-stage ROOT-C++ pipeline (no python). Histograms produced
here feed a separate HiggsAnalysis-CombinedLimit fork for the final fit. See
[CLAUDE.md](CLAUDE.md) for the architectural overview (what each macro does,
why histogram naming is load-bearing, the inter-channel asymmetries, etc.).
**This README is the runbook**: in what order to run things, what each step
reads, what it writes, and how to check success.

## Layout

| Directory          | Stage                            | Entry point                              |
| ------------------ | -------------------------------- | ---------------------------------------- |
| `merge_rootfile/`  | Build consolidated input file    | `make_filelist.sh` + `hadd_from_list.sh` |
| `skim/`            | Selection → per-sample histos    | `run_all.sh <channel> [samples]`         |
| `plotting/`        | Data/MC overlays, QCD bkg, plots | several macros — see Module 3 below      |
| `analysis/`        | Charge asymmetry, F/B ratio      | `charge_asym.C`, `FBratio.C`             |

`skim/` is unified post-refactor: shared utilities in `skim/skim_common.h`,
four channel functions plus a dispatcher in `skim/skim.C`. The pre-refactor
per-channel macros (`DrawW*`, `DrawDi*`) are preserved in `skim/legacy/` for
reference / fallback diffing.

`analysis/` uses a shared `analysis/analysis_helpers.h` for `YieldInRange`,
`AsymErr`, `RatioErr`, and the `kPORapidityShift = 0.3466` constant.

## Prerequisites

- ROOT 6 with ACLiC (`.C+` compilation)
- xrootd / EOS access for the input ntuple (CMS lxplus or equivalent)
- For the fit stage: CMSSW environment with `combine` and `text2workspace.py`
  on `PATH`, and your fork of HiggsAnalysis-CombinedLimit checked out to
  the `zheng/po-analysis` branch (see CLAUDE.md "Downstream fit")

---

# Workflow

Five modules, run roughly in order. Module 1 is one-time per data version;
modules 2–5 are the standard reanalyse-everything path.

## Module 1 — `merge_rootfile/`  (build the consolidated input ntuple)

Scope: gather per-sample ROOT files on EOS and `hadd` them into a single
`pO_2025.root` that the skim reads. Run this once per data/MC version.

```bash
cd merge_rootfile/

# 1.1  Discover files on EOS for each sample. Writes / updates the per-sample
#      filelists (MC_Wp_mu.txt, MC_DY_ele_Z.txt, DATA_pass_2.txt, etc.).
./make_filelist.sh

# 1.2  Inspect and edit the .txt lists if you want to subset / pin versions.
#      (The version_*.txt files are historical snapshots — don't edit those.)

# 1.3  Merge into the consolidated ROOT file. End result is published at
#      root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root
./hadd_from_list.sh
```

After this, you should not need to touch `merge_rootfile/` until the input
data changes.

## Module 2 — `skim/`  (event selection → per-sample histograms)

Scope: read the consolidated ntuple, apply the 8-step W cutflow (or the Z
selection), produce rapidity-binned MT / MET histograms per sample.

```bash
cd skim/

# 2.1  Sanity: confirm DATA_FILE in run_all.sh points at the right file.
grep '^DATA_FILE=' run_all.sh

# 2.2  Run a single channel × sample combo to smoke-test the build:
./run_all.sh Wmu Data
# This calls: root -l -q -b 'skim.C+(kWmu, "<DATA_FILE>", kData)'
# On first run, ACLiC will compile skim.C → skim_C.so (a few seconds).

# 2.3  Run a whole channel (all 7 samples):
./run_all.sh Wmu

# 2.4  Run everything (4 channels × 7 samples = 28 jobs, serial):
./run_all.sh all
```

CLI: `./run_all.sh <Zmm|Zee|Wmu|Wel|all> [samples...]` where samples is a
subset of `Data DY Wp Wm DYtau Wptau Wmtau` (default: all 7).

**Check success:**
- `skim/logs/skim_<channel>_<sample>.log` — per-job log; the driver prints
  `OK` or `FAIL` per job and exits non-zero if any failed.
- `skim/output/*_{Data,MC}.txt` — human-readable cutflow tables (W channels).
- `skim/rootfile/*_hist.root` — the histogram outputs. Names:
  - W: `WToMuNu_pO_PFMet_hist.root` / `..._<sample>_hist.root`
  - W (ele): `WToElecNu_pO_PFMet_hist.root` / `..._<sample>_hist.root`
  - Z: `ZToMuMu_pO2025_{Data,MC}_hist.root` (+ per-sample for MC)
  - Z (ee): `ZToEE_pO2025_{Data,MC}_hist.root` (+ per-sample for MC)

  Each W file contains 12 rapidity-binned MT/MET histograms per charge
  (`h_mt_Wp_y0..y11`, `h_mt_Wm_y0..y11`, ditto MET, plus the `_FB` variants
  and the iso-binned `h_met_iso_*_bin{0..2}` used by the QCD sideband fit).
  Each Z file contains `hMass`, `hMass_extended`, `hMass_vipul`.

**Optional isolation studies** (only when retuning a working point):
```bash
root -l -q 'isolation.C+'        # muon ROC inputs
root -l -q 'isolation_ele.C+'    # electron ROC inputs
```
These write extra ROOT files consumed by `plotting/PlotsIsoROC.C` and
`plotting/PlotIsoROC_ele.C`.

## Module 3 — `plotting/`  (data/MC, background estimation, intermediate plots)

Scope: read the skim outputs, build data/MC overlays, estimate QCD background,
prepare Combine inputs, draw publication plots. Some macros depend on others;
the recommended order is below.

All commands assume `cd plotting/`. Outputs land under `plotting/plots/` and
`plotting/plots/Elec/` (the Combine fork reads these via
`/afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots/`).

```bash
cd plotting/

# 3.1  QCD background shape (low-iso sideband fit, Rayleigh-like). Standalone.
#      Currently doesn't work well at the available statistics — replacing
#      with an ABCD method is on the TODO list.
root -l -q 'qcd_sideband_fit_and_extrapolate.C+'

# 3.2  Isolation ROC curves (only run if you re-ran skim/isolation*.C).
root -l -q 'PlotsIsoROC.C+(0)'           # muon
root -l -q 'PlotIsoROC_ele.C+'           # electron

# 3.3  Data/MC overlays for MT and MET in each rapidity bin.
#      ALSO writes plots/combine_input_inclusive.root (or plots/Elec/...),
#      which the Combine fork consumes.
root -l -q 'mtandmet.C+(0)'              # muon channel  -> plots/
root -l -q 'mtandmet.C+(1)'              # electron      -> plots/Elec/

# 3.4  Dilepton (Z) peak plots. ALSO writes plots/combine_input_dilepton.root.
root -l -q 'dileptonpeak.C+(0)'          # Z->mumu
root -l -q 'dileptonpeak.C+(1)'          # Z->ee

# 3.5  (Optional, cosmetic) muon/electron channel overlay.
root -l -q 'mtandmet_overlay.C+'

# 3.6  (Optional, cosmetic) Z mass-peak curves.
root -l -q 'plotZcurve.C+'

# 3.7  Theory predictions for R_FB (reads pQCDLightIon/ + filelist_theory.txt).
#      MUST run before observables.C in Module 5, because observables.C reads
#      its output RpO_rootfile/RpO_FB_graphs.root.
root -l -q 'plotRpOtheory.C+'

# Now switch to Module 4 (analysis/), then come back for 3.8 below.

# 3.8  Final-observable plots: charge asymmetry and R_FB, overlaying data
#      (from analysis/) with theory (from 3.7). Run AFTER Module 4.
root -l -q 'observables.C+'
```

**Files that are scratch / stale and should NOT be in the normal run order:**
`test.C`, `test111.C`, `Z_MC_overlay.C`. (Safe to delete eventually.)

**`plotting_helper.C` and `CMS_lumi.{C,h}`** are headers / helpers `#include`d
by the macros above. Don't run them directly.

## Module 4 — `analysis/`  (extract final observables)

Scope: read the W-channel rapidity-binned MT histograms from the skim, compute
charge asymmetry and forward/backward ratio per rapidity bin. Outputs feed
plotting/observables.C and (downstream) the cross-section fit.

```bash
cd analysis/

# 4.1  Charge asymmetry: A(y) = (N+ - N-) / (N+ + N-) per rapidity bin.
#      Reads ../skim/rootfile/WToMuNu_pO_PFMet_hist.root (and equivalents).
#      Writes ../skim/rootfile/charge_asym.root with a TGraphErrors per channel.
root -l -q 'charge_asym.C+'

# 4.2  Forward/backward ratio R_FB(|y_CM|). Same input, same output dir.
#      Writes ../skim/rootfile/FBratio.root.
root -l -q 'FBratio.C+'
```

Function signatures take optional kwargs (input/output paths, MT vs MET,
integration range, NY). Default args match the standard workflow — read the
function prototypes in the .C files if you need to override.

Once these have run, go back to Module 3 step 3.8 to draw the final plots.

## Module 5 — Cross-section fit  (HiggsCombine fork, separate repo)

Scope: read the Combine input histograms produced by Module 3.3 / 3.4, run
the statistical fit per channel, extract cross-sections and post-fit
distributions. Lives in a **separate** repo with the fit code on a non-default
branch.

```bash
# 5.1  Switch to the Combine fork (or open a new shell there).
cd /Users/zhenghuang/HiggsAnalysis-CombinedLimit

# 5.2  Check out the working branch (main is unmodified stock v10.6.0).
git checkout zheng/po-analysis

# 5.3  Source the CMSSW environment (combine + text2workspace.py on PATH).
#      Typical lxplus pattern:
#   cd /path/to/CMSSW_X_Y_Z/src && cmsenv && cd -
#      ... your setup may differ.

# 5.4  Make sure plotting/plots/ from this repo is reachable by run_fit.sh.
#      The fit reads from /afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots/.

# 5.5  Run the fit driver.
bash test/run_fit.sh
```

`test/run_fit.sh` invokes `test/my_script/make_combine_input{,_Z}.C` to build
Combine inputs, runs `text2workspace.py` on each datacard
(`testdatacard_{inclusive,Zmumu,Zee}.txt`), then `combine -M FitDiagnostics`,
and finally `draw_postfit_*.C` for postfit plots. Three channels are fit
independently (W inclusive, Z→μμ, Z→ee).

---

## Quick rerun cheat sheet

If everything is set up and you just want to redo from the skim onward:

```bash
cd skim                                  && ./run_all.sh all
cd ../plotting                           \
    && root -l -q 'qcd_sideband_fit_and_extrapolate.C+' \
    && root -l -q 'mtandmet.C+(0)'  && root -l -q 'mtandmet.C+(1)'  \
    && root -l -q 'dileptonpeak.C+(0)' && root -l -q 'dileptonpeak.C+(1)' \
    && root -l -q 'plotRpOtheory.C+'
cd ../analysis                           \
    && root -l -q 'charge_asym.C+' && root -l -q 'FBratio.C+'
cd ../plotting                           && root -l -q 'observables.C+'
cd /Users/zhenghuang/HiggsAnalysis-CombinedLimit \
    && git checkout zheng/po-analysis    && bash test/run_fit.sh
```

---

## Notes / health of the pipeline

- **Sumw2 status (audited May 2026):** all 84 histograms in `skim/skim.C`
  call `Sumw2()` explicitly at construction, so per-bin errors are tracked
  through Scale / Clone / Add / Integral. `skim()` also calls
  `TH1::SetDefaultSumw2(kTRUE)` defensively so any *future* histogram is
  Sumw2'd by default. **Caveat:** `analysis/analysis_helpers.h::{AsymErr,
  RatioErr}` derive errors from raw counts (Poisson). Today that's correct
  because skim fills are unweighted; if event weights ever enter the skim
  Fill calls, those analysis errors will silently become wrong and need
  to switch to `TH1::IntegralAndError`.

- **Corrections WIP.** Recoil corrections, lepton scale factors, momentum
  scale/smearing are not all applied yet. Don't assume MC in `skim/rootfile/`
  is fully corrected when reading downstream.

- **Pre-existing inter-channel asymmetries** (preserved through the May 2026
  refactor — see CLAUDE.md): DY-veto pT cut differs between Wmu (15 GeV) and
  Wel (10 GeV); isolation cut differs (0.15 vs 0.095); DY-veto mass-window
  comments lie ("> 30" but actual is `80 < m < 110`); `tHi->GetEntry` is
  unguarded in Zmm but guarded in Zee. These may all be intentional but
  are worth a second look.

- **No tests, no CI.** Validation is by re-running the skim and inspecting
  cutflow logs + plots. Recent commit messages are placeholders ("xx") so
  `git log --oneline` is unreliable — read diffs.
