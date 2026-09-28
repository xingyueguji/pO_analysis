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
- [skim/skim_common.h](skim/skim_common.h) — shared utilities, enums (`SampleType`, `ChannelType`), kinematics, cutflow helpers, pO event filters, etc. Wrapped in `namespace pOSkim`. **Since 2026-09-14 also the electron ECAL-crack veto: `kEcalGapLo/Hi` (1.4442/1.566), `InEcalGap(scEta)` and `BuildEleIDNoGap(nEle, eleID, eleSCEta, out)` = the effective electron ID (MVA WP95 AND not in the crack) that every ID'd-electron definition uses — see selection step 5.**
- [skim/skim.C](skim/skim.C) — one file, four entry functions plus a dispatcher:
  - `skim_Wmu(file, sample)` — W→μν
  - `skim_Wel(file, sample)` — W→eν
  - `skim_Zmm(file, sample)` — Z→μμ
  - `skim_Zee(file, sample)` — Z→ee
  - `skim(channel, file, sample)` — dispatcher used by `run_all.sh`
- [skim/legacy/](skim/legacy/) — the original four per-channel macros, preserved verbatim for diffing/fallback. The unified files were byte-identical in physics output to these until 2026-07-06, when the relIso definition gained the Δβ PU correction (legacy keeps the uncorrected sum — expect small selection diffs).
- [skim/count_ngen.C](skim/count_ngen.C) + [skim/run_ngen.sh](skim/run_ngen.sh) — compute N_gen (Σ gen weight over **all** events, no selection) per MC sample → `skim/rootfile/ngen.root` + `skim/output/ngen.txt`. The cross-section-normalization denominator. See "MC normalization" below.
- [skim/gen_xsec.C](skim/gen_xsec.C) — **GENERATOR-level W σ per rapidity bin (2026-08-05; FIDUCIAL since 2026-08-12)**: loops ALL generated events of the 4 W MC files (no reco, no selection; needs the `vector<vector<int>>` dictionary like the skim), picks the gen charged lepton of the sample's flavour+charge (highest-pT, W-ancestor preferred; falls back with a WARN — the ntuple mother chain often omits the W), histograms its LAB η in `kYEdges`(+FB) with the POWHEG weight → `skim/rootfile/gen_xsec.root`. **Binned in `y = −η_lab`, the skim's convention** ("p-going (−Z) = forward", `skim.C:835,1497`) — **BUG FIXED 2026-08-12**: it previously filled raw +η, so gen bin i paired with the MIRRORED reco region y_i. Caught by the per-bin (A·ε) diagnostic (mirrored pairing made A·ε charge-dependent — μ W⁻ ran 1.27→0.73 across η — instead of the charge-symmetric detector response it must be; reversing the gen bins collapsed the W⁺/W⁻ spread from 0.54 to 0.06). Effect: per-bin σ was mirrored, inclusive σ shifted ≤2% (W⁻ met 40.15→39.68, leppt_mt40 44.86→43.89); counts/charge_asym/FBratio never read this file and were unaffected. **`h_gen_sig_{Wp,Wm}[_FB]` are FIDUCIAL since 2026-08-12**: gen (bare, post-FSR) lepton pT > 25 (`kFidPtMin` = the skim's nominal cut; |η_lab| < 2.4 implicit in the binning window; NO m_T cut — one fiducial serves all discs), so **σ_meas,i = r_i × σ_gen-fid,i directly** (= yield-based extraction with MC A·ε, algebraically identical; kA_O and kSigma cancel between r and σ_gen — see the r×σ_gen plan under "Future plans"). `h_gen_tot_{Wp,Wm}[_FB]` keep the old no-pT-cut definition — but **NOT a true total**: the ntuple gen collection `HiGenParticleAna/hi` is itself filtered at **pT > 5 GeV and |η| < 2.5** (measured 2026-08-12), so fid/tot = 0.675 (W⁺) / 0.755 (W⁻) is a pT>25/pT>5 ratio, not an acceptance. **σ_fid is unaffected** (the whole fiducial sits inside the filter and Σw runs over all events) — the filter is also why gen_xsec.C WARNs that 25–37% of events have "no gen lepton" (those are below 5 GeV or beyond |η| 2.5, outside the fiducial either way). Normalization is the weighted FRACTION × `mc_norm.h` σ: σ_i = kA_O·kSigma·Σw_i/Σw_all — **unit-proof, because the July-29 `weight` is σ in pb (⟨w⟩≈6376) while `kSigma_*` is nb**; the pipeline was always consistent since σ cancels in k_s = A·σ·L/Σw. Discriminant-independent one-time producer consumed by `xsec_fiducial_comb` (gen overlay). Gen storage window extends past |η|=2.4 (nonzero over/underflow → edge bins safe); per-flavour σ: fid(pT>25, |η|<2.4) W⁺ 49.8 / W⁻ 39.6 nb; in-window no-pT-cut 73.7/52.5; total 102.0/87.4. **Z GEN FIDUCIAL + a plain-text sidecar (2026-09-15):** `h_gen_sig_Z` (ONE bin) = σ(Z→ll) in the gen twin of the skim's dilepton selection — lead pT > 15, sub > 10, |η| < 2.4, 60 < m_ll < 120, highest-pT OS same-flavour gen pair, DY_mu and DY_ele pooled like the W's — so σ_Z = r_Z × h_gen_sig_Z is the Z axis of the (σ_W, σ_Z) plane (`plotting/xsec_contour.C`). **σ_Z^fid = 8.757 ± 0.006 nb per flavour**, acceptance 0.4658 of σ(DY→ll, m>50) = 16 × 1.175 = 18.80 nb; **μ 0.4703 vs e 0.4615 — the ~2% gap is bare post-FSR electrons losing more energy**, and it is the reason the two are quoted separately in the log. 46.7% of events have no OS gen pair at all: the same `HiGenParticleAna/hi` filter (pT > 5, |η| < 2.5) as the W's 25–37%, and equally harmless (the whole fiducial sits inside the filter, and Σw runs over ALL events). No Z-ancestor preference is attempted — the `hi` tree stores no boson, so `HasAncestor` is a structural no-op there (see `charge_flip.C`). Also writes **`skim/output/gen_xsec_fid.txt`** (48 W rows `lab|fb <charge> <ybin> <sigma_nb>` + `incl Z 0 <sigma_nb>`): the same numbers in a form the FORK's card generator reads with awk, because the reparametrized workspace that promotes σ_W to a POI must bake in the per-bin σ_gen,i. Uploaded by `sync_lxplus.sh upload` to `<ana>/skim/output/`. **EPPS21 MEMBER CROSS SECTIONS (2026-09-15c):** `h_gen_sig_{Wp,Wm}[_FB]_epps21` (TH2D, x = rapidity bin with the nominal's own axis, y = member) and `h_gen_sig_Z_epps21` (TH1D over members) = the gen fiducial σ for all 107 members of `EPPS21nlo_CT18Anlo_O16`, via `pOLhe::ComputeMemberWeights` on `ttbar_w` (the same member table the skim's reco twins use). **The normalization is the one thing that matters here:** the nominal is a weighted FRACTION of a FIXED total (Σw_i/Σw_all), but a PDF variation changes the total too, so dividing by that member's OWN Σw_all would divide exactly that change out and collapse every member onto the same σ. The members are therefore normalized to the **NOMINAL** total, σ_i(m) = kA_O·kSigma·Σw_i(m)/Σw_all(0), which makes σ_tot(m)/σ_tot(0) = Σw_all(m)/Σw_all(0) exactly and reproduces the nominal for m = 0 (**verified: max|member0 − nominal| = 0 per bin, W⁺, W⁻ and Z**). Accumulated in plain arrays, not TH2D::Fill — 7M events × 107 members is far too many fills. Member spreads: σ_total W⁺ 98.40–103.73 nb (nominal 102.02), W⁻ 84.71–88.67 (87.42), DY 18.19–19.10 (18.80); on the FIDUCIAL σ the 106 variations span W −3.2/+1.3% (lab) or −3.7/+1.6% (fb) and Z −3.1/+1.3%. NB these are individual-member EXTREMES, NOT the Hessian uncertainty (a quadrature over eigen-directions, ~4.6% at 90% CL — necessarily larger). Adds ~2× to the job's run time (reading 217 floats/event). Consumed by `xsec_contour.C`'s member scatter; it is also the first half of the long-deferred "gen twins" TODO.
- [skim/lhe_weights.C](skim/lhe_weights.C) + [skim/run_lhe_weights.sh](skim/run_lhe_weights.sh) — **decoder for `hiEvtAnalyzer/HiTree::ttbar_w` (2026-09-01), the per-event LHE systematic weights.** `ttbar_w` (`vector<float>`, **217 entries/event, identical layout in all 9 MC files**) is the FULL LHE `<rwgt>` block rescaled so element 0 == `weight` (`HiEvtAnalyzer.cc`: `ttbar_w[i] = weight/originalXWGTUP × LHEEventProduct::weights()[i]`; the name is historical, nothing to do with ttbar). The ntuple carries no `<initrwgt>` header; from the numbers alone three members are identical to the nominal (idx 0, 44, 110 — PDF centrals reproducing the generation PDF), splitting the vector into **regions A = 0–43, B = 44–109, C = 110–216** (the default figure shows exactly that, no names — user decision 2026-09-01). **IDENTITIES CONFIRMED 2026-09-02 against the producer's header generator, CMS POWHEG `make_rwl.py`** (user supplied the script and then the whole genproductions `bin/Powheg` archive a1a26254; **the producer confirmed the `"EPPS21" in Period` branch**). TWO mechanisms write the 217: **(1)** `make_rwl.py` → `pwg-rwl.dat` `<initrwgt>` read at generation = 3×3 scale grid (ids 1001–1009, `lhapdf=14600` = `run_pwg_condor.py`'s defaultPDF for an `EPPS21_*` ion) + a "hessian" PDF group (35 single proton centrals + CT18ANLO 59 + 4 CT18ANLO α_s) + a "replica" group (3 singles) = **110 weights = idx 0–109**; `lhapdf=X` swaps the LHAPDF set of BOTH beams (the card must have `lhans1 == lhans2 = 14600`). **(2)** `runcmsgrid_powheg.sh` AFTER generation: if the card has `nPDFerrSet`, it reruns `pwhg_main` 107 times with `rwl_add 1`, once per `nPDFerrSet = 1..107`, appending **ids 9001–9107 (weightgroup `EPPS21_variation`, combine=hessian) = idx 110–216**. The nuclear PDF itself is a POWHEG source patch (`patches/EPPS21/*`): for the beam with `ia ≥ 16` the LHAPDF proton PDF is multiplied by the EPPS21 R-factors of error set `nPDFerrSet` (`EPPS21.f` + the `EPPS21NLOR_16` grid) and isospin-averaged (Z/A mixing u↔d, ū↔d̄), so **f_O = R(set) × f_p^LHAPDF, and the two weight kinds vary the two factors SEPARATELY** — (1) varies f_p on both beams with R at set 1, (2) varies R with f_p fixed at CT18ANLO central. Hence idx 44 (`lhapdf=14600`) and idx 110 (`nPDFerrSet=1`) are both exactly the nominal configuration. Every testable feature matches: the α_s scans are monotonic (idx 11–18 NNPDF31 0.108→0.124, 23–28 NNPDF40, 30–34 CT18NNLO, 103–106 CT18ANLO — the last is linear at 0.47%/0.001), 38 == 39 (MSHT20 as118 ≡ as_smallrange m0), 107 == 9 (NNPDF31 replica mean ≡ mc_hessian central), 108 ≈ 22 (NNPDF40 replica vs hessian central), the scale grid obeys the additive/cross-term algebra with idx 1,2,3,6 = single-scale variations, and the EPPS21 nuclear/proton sub-blocks change character at 158/159. **The map (`kSets` in the macro, printed per index in the txt):** **0** id 1001 μR=μF=1 (nominal); **1–8** ids 1002–1009, (μR,μF) ∈ {1,2,0.5}² minus (1,1), inner loop over μF (1=(1,2) 2=(1,½) 3=(2,1) 4=(2,2) 5=(2,½) 6=(½,1) 7=(½,2) 8=(½,½)) — μF×2 +3.1% / μF×½ −4.7% / μR×2 −1.3% / μR×½ +1.7% (W⁺), envelope +4.0/−6.6% W⁺, +3.8/−6.5% W⁻, +3.6/−6.0% DY; per-event ratios have RMS 0.06–0.24 and flip sign (they RE-SHAPE; PDF members only re-scale, RMS 0.002–0.06); **9–43** 35 proton-PDF centrals (NNPDF3.0/3.1/4.0 + α_s scans, CT18NNLO(+α_s)/Z/A/X, MSHT20 ×3, PDF4LHC21, HERAPDF20 (= the +6.6% outlier at 42), ABMP16), +7.4/−4.2%; **44** CT18ANLO m0 (== the nominal LHAPDF set) + **45–102** CT18ANLO eigenvectors 1–58 on both beams, R central (29 pairs — THE proton-PDF uncertainty) + **103–106** CT18ANLO α_s 0.116/0.117/0.119/0.120 (±0.5/±0.9%); **107–109** replica-group centrals NNPDF31_nnlo_as_0118_mc / NNPDF40_nnlo_pdfas / NNPDF40_nnlo_pch_as_01180 (+3.0/+3.8/+4.9%; NOT an α_s scan as first guessed); **110** id 9001 `nPDFerrSet=1` = central EPPS21 R (== nominal; DY −0.008%) + **111–158** ids 9002–9049, EPPS21 R-factor sets 2–49 = the 24 nuclear eigen-direction pairs (92% opposite-signed) + **159–216** ids 9050–9107, sets 50–107 = the 29 CT18A-baseline eigen-direction pairs of R — **R only: the baseline PDF itself is NOT moved here** (that is what idx 45–102 does); the coherent EPPS21 baseline variation f_A,k = R_k × f_p,k is to first order the PRODUCT of the idx 44+k and idx 158+k per-event ratios. 9+35+59+4+3+107 = 217. **Nominal = CT18ANLO on both beams × EPPS21 R (set 1, O16, isospin-averaged) on the oxygen beam** (= EPPS21nlo_CT18Anlo_O16 central). Inclusive-σ Hessian uncertainties (90% CL as delivered; ÷1.645 for 68%): **EPPS21 nuclear sym 4.6% W⁺ / 3.8% W⁻ / 4.0% DY** (asym +3.3/−6.2% W⁺), EPPS21 baseline-direction R-only sets 2.3/2.0/2.1%, CT18ANLO eigenvectors (both beams) 2.5/2.4/2.4% (asym +2.5/−3.2% W⁺). On the **W⁺/W⁻ ratio** these collapse to nuclear 0.85%, EPPS21 R-baseline 0.29%, CT18ANLO 0.55%, scale envelope +0.15/−0.11% (the 35 alternative centrals move it by −0.4% on average, at most −0.8%); on the **charge asymmetry** (A₀ = 0.077) the absolute Hessian shift is 0.0042 (nuclear) / 0.0015 (proton), envelope −0.003/+0.001 — nPDF effects largely cancel in the ratio observables. Loops ALL events (no selection, like `count_ngen.C`; inputs from `ResolveMCSample`; `./run_lhe_weights.sh [mu|ele] [nmax]`, ~5 min for the three ~0.6 GB branches, log `skim/logs/lhe_weights_<flav>.log`) → `skim/output/lhe_weights.txt` (per-index table with set/member/weight-id/LHAPDF-id of every weight + block summaries incl. W⁺/W⁻ and A), `skim/output/lhe_weights_structure.png` (two-pad figure: σ_i/σ_0−1 per index with regions A/B/C + the per-event std. dev. — NB a plain unweighted std. dev., tail-dominated in region A by the near-zero-weight NLO sign-flip events: MAD-based spread 0.025 vs 0.157 at idx 5), `..._named.png` (same + the confirmed content per region and dashed sub-block boundaries) and `skim/rootfile/lhe_weights.root` (`h_wrel_/h_wrms_/h_sumw_<label>`); `./run_lhe_weights.sh draw` redraws both figures from the rootfile in seconds. **Analysis relevance:** σ_meas = r×σ_gen-fid depends on the PDF only through the MC A·ε and template shapes (kSigma is fixed, the inclusive shift cancels) — the natural next use is the per-bin A·ε variation, filling the skim histograms with `w·ttbar_w[i]/ttbar_w[0]`; idx 111–158 (EPPS21 R sets 2–49) are the nPDF systematic, idx 45–102 (CT18ANLO 1–58 on both beams) the proton-PDF one (idx 159–216 add only R's baseline dependence; combine as the product of the paired ratios for the full EPPS21 baseline variation), ids 1002–1009 the scale envelope. **Combination rule cross-checked against LHAPDF source (2026-09-02):** `PDFSet::uncertainty(values, cl=CL1SIGMA)` for `hessian` = the asymmetric Hessian over pairs (2i−1, 2i) with max(·,0) — errplus/errminus, rescaled to 68% by default via √(χ²q(68.27%)/χ²q(90%)) = 1/1.645 — identical to `lhe_weights.C::Hessian`'s up/dn ÷ 1.645, bin by bin. **NB its `errsymm` is NOT the symmetric-Hessian ½√Σ(v₂ᵢ₋₁−v₂ᵢ)² (that is `Hessian()`'s `sym`, the "sym" column of `lhe_weights.txt`): measured 2026-09-08 with LHAPDF 6.5.6 on a toy 107-vector, `errsymm = (errplus + errminus)/2` for an asymmetric `hessian` set** (2.1278 = (3.0398+1.2159)/2; the pipeline only reads errplus/errminus, so nothing depends on it); a colleague's pPb code calls exactly this on the 107 EPPS21 values (their `_nPDF` systematic). EPPS21 paper (arXiv:2112.12462 §4.2–4.3): 24 parameters, Δχ² ≈ 33 = 90% CL; baseline sets = the nuclear fit repeated with each CT18ANLO error set S_i^± as baseline (signs follow CT18A's), total = nuclear ⊕ baseline in quadrature (Eq. 39); **Eq. (40) defines the proton-error cross sections as σ(S_i^±) = f^p_{i−24,±} ⊗ σ̂ ⊗ f^A_{i,±} for i = 25–53 — the proton-beam PDF is the CT18A error set AND the nuclear PDF is the refit with that baseline**, so the EPPS21 block alone (our idx 159–216 = R_i only) is NOT the CT18A variation; the product with idx 44+m reproduces Eq. (40) exactly. With EPPS21.f numbering (even = S⁻, odd = S⁺) and the CTEQ convention (odd LHAPDF member = "+"), the coherent product pairing is CROSSED (idx 44+(2i−1) ↔ 158+2i, 44+2i ↔ 158+(2i−1)) → baseline term 3.0% (W⁺) instead of 3.7% — CTEQ sign convention still to be verified. Origin of every weight is documented: the relevant generator scripts are copied READ-ONLY to `skim/reference/genproductions_a1a26254/` (README there; `make_rwl.py`, `runcmsgrid_powheg.sh`, `run_pwg_condor.py`, `runGetSource_template.sh`, `patches/EPPS21/*`, `MetaData/npdflist_O_5f_run3.dat`), and `lhe_weights.txt` opens with those paths + the original zip, a PROVENANCE BY BLOCK section (idx range, the producing file:line and code entry, why the block exists) and per-row columns "written by (file:line)" + role, so every one of the 217 entries is traceable (user request 2026-09-02). Nothing left to ask the producer. **HOW TO APPLY THE PDF/nPDF SYSTEMATIC — agreed with the user 2026-09-07 (after checking the LHAPDF source, EPPS21 Eqs. 39–40 and the data; the user questioned the two-block construction several times — the decisive facts are below):** the EPPS21 uncertainty = **53 directions**: **24 nuclear pairs = idx (111,112)…(157,158)** ⊕ **29 baseline pairs whose members are the PRODUCT of the per-event ratios idx 44+m and idx 158+m (m = 1..58)** — Eq. (40) puts the CT18A member on the proton beam AND inside the nuclear PDF, block C alone is R-only (2.3% on σ_W⁺ vs 3.0–3.7% correct), and there is **NO separate CT18 term** (the 29 directions enter once — adding CT18ANLO independently would double count). Data evidence that block C is R-only: the event-by-event spread of its baseline members is **0.33–0.37× block B's** (W⁺/W⁻/DY; ≥0.7 expected if C carried the CT18A member on the O beam). Per bin: asymmetric Hessian (per pair the larger up / larger down move, quadrature over pairs — identical to LHAPDF `PDFSet::uncertainty`) then **÷1.645** (90%→68%); an envelope over Hessian members is WRONG (≈2× low). Scale = envelope of idx 1–8 (7-point drops idx 5 and 7), no 1.645; α_s = idx 103–106 separately (±0.5%/±0.001); idx 9–43 and 107–109 never enter an uncertainty. **Open:** the sign pairing inside each baseline pair (EPPS21.f even = S⁻/odd = S⁺ vs CTEQ odd member = "+" ⇒ probably CROSSED, idx 44+(2i−1) ↔ 158+2i; 3.0% vs 3.7% on W⁺) — verify from the CT18 paper/.info or ask the EPPS21 authors; conservative fallback = the larger. **APPLIED 2026-09-07 — SKIM STAGE (see the `lhe_index.h` and `lhe_updown.py` bullets below):** the reco fit templates carry per-member twins (RAW variation = nominal template with the extra weight, nothing divided out — user decision) and every MC skim file holds `<h>_nPDFUp/Down` from LHAPDF's `PDFSet.uncertainty()`, `<h>_qcdScaleUp/Down` (μR/μF envelope over all 9 points, per-bin max/min) and `<h>_alphaSUp/Down` (the 0.119/0.117 member templates). **Deferred by the same decision (TODO):** the gen-fiducial denominators in `gen_xsec.C` (the "acceptance-only" variation consistent with σ = r×σ_gen — the raw variation also shifts each template's normalization, degenerate with r_i), the closure/A·ε diagnostic (member 0 vs nominal is checked; idx 44/110 are not stored), the per-eigen-direction nuisances (53 U/D pairs = members ÷1.645, the correlation-preserving alternative to one collapsed `nPDF`), the baseline sign-pairing verification, and the downstream consumers (`mtandmet.C` Up/Down writer, fork `shape` rows).
- [skim/lhe_index.h](skim/lhe_index.h) — **SINGLE SOURCE for the `ttbar_w` layout AND the STORED-MEMBER FAMILIES (2026-09-07)**, namespace `pOLhe`. Part 1 = `kNW`, `kBlocks`/`BlockOf`, `kSets`/`SetOf`/`SetLabel`, `kScaleLabel`, `HessRes`/`Hessian()` moved VERBATIM out of `lhe_weights.C`'s anonymous namespace (the macro now `#include`s it; `lhe_weights.txt` regenerated byte-identical after the move). Part 2 = the three families the skim stores as TH2D twins of every FIT-TEMPLATE histogram (x axis copied from the nominal, y = member index, MC only): **`<h>_epps21` = the 107 members of the LHAPDF set `EPPS21nlo_CT18Anlo_O16`, member for member** — 0 nominal (`ttbar_w[0]`), 1–48 = idx 111–158 (the 24 nuclear pairs), 49–106 = the 29 baseline pairs built PER EVENT as `(ttbar_w[44+k]/w0)·(ttbar_w[158+k′]/w0)`, k = 1..58, k′ = the other member of k's pair under **`kBaselinePairing = kCrossed`** (the coherent pairing if EPPS21.f even=S⁻ / CTEQ odd="+" — still the open item; switching = edit the constant + a ~15 min MC re-skim); **`<h>_scale`** = idx 0–8 (the (μR,μF) grid); **`<h>_alphas`** = idx 0, 103–106. Per-event member weight = `w·ttbar_w[a]/ttbar_w[0] [·ttbar_w[b]/ttbar_w[0]]` (`ComputeMemberWeights`, once per event; vector null/≠217 → warn once + skip, `ttbar_w[0]==0` → skip; measured 0 skipped events in every MC file) — so **member 0 is BIT-IDENTICAL to the nominal histogram** (verified per cell incl. errors) and any SF later folded into `w` propagates. **FREAK-WEIGHT GUARD `kMaxMemberRatio = 10` (2026-09-15):** if max |ρ| over the 121 stored members exceeds 10, `ComputeMemberWeights` sets EVERY member weight to `w` (ρ → 1) and bumps the caller's counter — the event keeps its full nominal weight, only its variations die. **Why:** POWHEG generates UNWEIGHTED events, so `|ttbar_w[0]|` is ONE number per file (5538.8 in July-29) whose SIGN is the sign of B̄ = B + V + ∫R (0.69% negative); a variation weight is the same event re-evaluated, `ttbar_w[k] = ttbar_w[0]·B̄_k/B̄_0`, and **unweighting erased B̄_0**, so an event on a near-cancellation of B̄_0 gets a huge ratio **in every member at once** (they share the denominator). The proof that the denominator is at fault, not any member: on the worst such event **α_s ± 0.001 gives ρ = 11.5 / −10.2**, which no physical variation can do. Per-EVENT, not per-member — the cause is common, and varying some members but not others would distort the member-to-member correlations the envelope and the Hessian rely on. **Calibration:** |ρ| > 10 is a uniform ~0.002% population in all 8 MC files (max |ρ| per file 34–261); the guard fires **41–69 times per file (0.0035–0.0058%)** and shifts the inclusive sum of any of the 114 used members by **≤ 0.082% (mean 0.006%)**. Member 0 has ρ ≡ 1 by construction ⇒ **every nominal histogram stays bit-identical** (`compare_reskim.sh rootfile_pre_rhoguard`, 2026-09-15: every DIFF in every MC file is a `_epps21`/`_scale`/`_alphas` twin, 0 nominals, data files IDENTICAL; the same 44 of 96 templates move in all three families, which is the per-event signature). Counter reported in each job's `[INFO] LHE member twins booked:` line. The thing it removed: `July_29_MC_Wm_mu` **entry 70485** (w0 = −5538.8, μF×2 ρ = 228, nPDF baseline product 261; reco muon pT 36.5, η −0.23) supplied **45% of the qcdScale integral shift** of the leppt_mt40 W⁻ **y5 (fb) / y6 (lab)** signal template — see the `syst_shapes.C` bullet. **Diagnostic that found it, reusable with no re-skim:** `N_eff = content²/error²` per (bin, member) from the stored twins, and `frac = N_eff(member)/N_eff(nominal)` — a real variation reweights every event alike so ρ cancels and frac ≈ 1 (measured 0.91–1.00), a one-event artifact collapses it (0.03–0.40), nothing in between (the outlier enters the content linearly but Σw² quadratically). **Tooling NB:** `rdfentry_` is NOT a reliable global entry number under `EnableImplicitMT` — run single-threaded when you need entry numbers. API: `BookTwins(nominal, registry)` (name = nominal + `_epps21|_scale|_alphas`, the layout string stamped into the title, Sumw2), `FillTwins(twins, x, mw)`. **Why member templates at all:** LHAPDF's `PDFSet::uncertainty()` takes the 107 per-member VALUES of one bin (nonlinear — max/min per pair, then squares), so no per-event "up weight" exists and the 107 stored members are its irreducible input; 217 would be 110 too many (the comparison centrals, the raw CT18/R-only blocks — consumed by the products — and the two nominal duplicates never enter an uncertainty).
- [skim/lhe_updown.py](skim/lhe_updown.py) + [skim/run_lhe_updown.sh](skim/run_lhe_updown.sh) + [skim/lhe_env.sh](skim/lhe_env.sh) — **the theory Up/Down templates (2026-09-07): `nPDF` via LHAPDF's official `PDFSet.uncertainty()` on the `_epps21` twins, `qcdScale` and `alphaS` from the `_scale`/`_alphas` twins.** nPDF = the exact loop of the colleague's pPb `prepareSystVariation` `_nPDF` branch: per template and per x-bin `val = [member 0..106 bin contents]` → `err = lhapdf.getPDFSet("EPPS21nlo_CT18Anlo_O16").uncertainty(val)` (asymmetric Hessian over consecutive pairs; the set's `.info`: `ErrorType hessian`, `ErrorConfLevel 90`, rescaled by LHAPDF to 68.27% CL by default) → `<h>_nPDFUp = max(v₀+errplus, 0)`, `<h>_nPDFDown = max(v₀−errminus, 0)` (nominal bin errors kept, under/overflow untouched), written back INTO the MC skim file (UPDATE, overwrite — must be re-run after every re-skim, RECREATE wipes them). Prints per template the integral shifts and `max|member0 − nominal|` (must be 0). `./run_lhe_updown.sh [Wmu|Wel|Zmm|Zee|all]` → `skim/logs/lhe_updown_<chan>.log` (the record). PyROOT + `lhapdf` in ONE interpreter needs LHAPDF built against the PyROOT python: **LHAPDF 6.5.6 from source under `~/local/lhapdf`, `PYTHON=/opt/homebrew/bin/python3.12`** (= `root-config --python-version`; the `davidchall/hep` brew formula binds python 3.10 and cannot share a process with ROOT, and no pip wheel exists for macOS) + the O16 set untarred into `share/LHAPDF/` (the shipped `pdfsets.index` predates it) — recipe in `lhe_env.sh`, which exports `PATH/PYTHONPATH/LHAPDF_DATA_PATH/LHE_PYTHON`. First numbers (W⁺ MC, 68% CL, integrals): `h_leppt_mt40_Wp_y5` +2.5/−3.4%, `h_met_Wp_y0` +2.6/−2.1% (the all-events inclusive is nuclear 4.6% ⊕ baseline 3.0% at 90% ≈ 3.3% at 68%). **NB one collapsed `nPDF` nuisance is 100% correlated across bins/charges/channels**; the per-eigen-direction alternative (53 U/D pairs = members ÷1.645) needs no re-skim. **`qcdScale` (same day, user recipe): per bin Up = max / Down = min over the (μR,μF) members INCLUDING member 0** (so Up ≥ nominal ≥ Down bin-wise), over **ALL 9 points** (user decision; `--scale-points 7` gives the usual 7-point set that drops the (2,½)/(½,2) corners — those two corners ARE the envelope's extremes inclusively, +4.0/−6.6% vs +3.1/−4.7% without them), NO 1.645 — W⁺ signal `h_leppt_mt40_Wp_y5` +3.8/−6.7%, `h_met_Wp_y0` +4.1/−6.4% (integrals of the per-bin extremes; verified bin by bin: Up = max, Down = min, Down ≤ nominal ≤ Up in all 5280 bins of a W file). **`alphaS`: Up = the member-3 template (α_s 0.119), Down = the member-2 template (0.117)** taken as-is (user decision; ±0.001 = the PDF4LHC21 68% CL shift; `--alphas-up/--alphas-down` override) — ±0.47%. `--systs nPDF,qcdScale,alphaS` selects; all three are written by default, so one MC file gains 6 TH1Ds per template (576 per W file, 6 per Z file). **Consumed downstream since the same day:** `mtandmet.C`/`dileptonpeak.C` carry them into every Combine input as `<process>_<syst>Up/Down` + the `_systs.txt` sidecar (see "Structured inputs" under "Downstream fit"); next = the fork's three `shape` rows. **RECIPE CHANGED 2026-09-14 (user's final decision after the 09-07 fit showed the raw variation leaking the theory normalization into σ — see the "Current state" bullet): (i) `--norm reco` (DEFAULT, all three families): every member template is AREA-NORMALIZED to the nominal integral over the fit bins 1..nx BEFORE the family's combination (member k × I_0/I_k), i.e. the theory-cross-section change is divided out and only the shape variation reaches the fit; the acceptance×efficiency part (= normalizing to the GEN-level member integral instead, which needs the gen twins in `gen_xsec.C`) is deliberately deferred — `--norm none` = the pre-09-14 raw variation; (ii) `qcdScale` = envelope over the 6 non-antagonistic (μR,μF) points + nominal (`--scale-points 6`, drops (2,½)/(½,2); `8` = the old all-points choice); (iii) `alphaS` SYMMETRIZED (`--alphas-mode symm`): per bin err = (N̂[0.119] − N̂[0.117])/2, Up/Down = nominal ± err (no one-sided bins by construction; `pick` = the old two-member recipe). Per-template log lines now also print `I_k/I_0 min..max removed` (the normalization that was divided out; nPDF 0.937–1.035 in y11, scales 0.951–1.032, α_s 0.995–1.005). Resulting integral shifts of the stored Up/Down: nPDF ±0.4%, qcdScale ±0.3%, alphaS 0.00 (the residual of a per-bin combination of area-normalized members) vs the raw +2.5/−3.4% (+5.7/−8.2% in y11), +3.8/−6.7%, ±0.5%. Re-run over all 24 MC files + all four Combine inputs regenerated the same day (nominal templates bit-identical to the last-fitted copies, `compare_hists.C`: 427 IDENTICAL / 0 DIFF per W file; sidecars unchanged, so the fork needs NO change); `syst_shapes.C`'s closure recompute applies the same per-member factors (detected from the stored titles) and its `[INCL]` block prints a note that (b) is now the shape-only residual. NOT YET REFITTED (lxplus).**
- [skim/muon_sf.h](skim/muon_sf.h) + [skim/sf/](skim/sf/) — **MUON EFFICIENCY SCALE FACTORS folded into the MC event weight (2026-09-14, SF-application phase step 2; user decision: pp POG SFs as a tentative stand-in for pO tag-and-probe, trigger = our own MB-derived SF, electrons later).** Single source (namespace `pOSF`) for what is applied, from where, how it is looked up and how the variations are formed — `skim.C` only calls `W()`/`Z()` and fills twins. **ID** = `NUM_TightID_DEN_TrackerMuons` (our `muIDTight`), **ISO** = `NUM_TightPFIso_DEN_TightID` (PF relIso(Δβ, R=0.4) < 0.15 | TightID = exactly the W cut), from the Muon POG pp-2025 correctionlib file: signed-η (24 × 0.2) × pT [10,15,20,25,30,40,50,60,120,∞) bins, `nominal` + the XPOG `systup/systdown` (= nominal ± √(stat²+syst²), verified). The committed copy `skim/sf/muon_sf_2025_TightID_PFIso_schemaV2.json` (0.45 MB) is the user's 14 MB `ScaleFactors_Muon_ID_ISO_2025_schemaV2.json` (md5 5a1a6d1bc4611102f629a2cd2e1646cf, 28 corrections) reduced by `skim/sf/extract_muon_sf.py` to the two used corrections (nothing else — no diagnostic LoosePFIso table), a provenance paragraph appended to the file's description and echoed in every job log. The script is a one-off data-prep utility (Python stdlib, kept as Python — user 2026-09-25); nothing in the pipeline runs it, it is only re-run for a new POG file or a different WP. The header has its own ~100-line JSON reader (correctionlib is not installed for the ROOT python; NB the reason once given here, "ROOT 6.32 ships no generic parser", is wrong — nlohmann-json 3.11.3 comes with the Homebrew ROOT and ROOT 6.32 ships `RooFit/Detail/JSONInterface.h` — but the reader is tested and stays); lookups reproduce the Python reference values exactly (tested), out-of-range η/pT clamp to the edge bin and are counted (0 in every job). **TRIG** = the fired-and-matched SF of `correction/rootfile/trig_eff_mb_mu.root` (the `mt40` variant = the purer W sample; `nom` gives 0.9962), Clopper-Pearson on the data count vs MC. **APPLIED AS ONE INCLUSIVE NUMBER, 0.9971 (+0.0020 −0.0024)** = 2521/2547 data vs MC 0.9927 (`kTrigBinning = kTrigInclusive`) — **user decision 2026-09-21, after the 2026-09-15 per-|y| version; the rapidity structure is REAL and is still MEASURED, it is simply not corrected for.** Everything needed to justify that is kept and regenerated on every run: `muon_sf.h` still reads the 6-bin |y| table (`sf_absy_mt40` + `h_{den,num,bit}_absy_mt40_*`) and prints it in the `[SF]` block marked **`per-bin table (NOT applied)`** with its pulls, plus the signed-y folding cross-check; `correction/trig_eff_mb.C` still produces `eff_absy_<sel>[_charge]`, `eff_y_<sel>[_charge]`, the `sf_absy_/sf_y_` graphs, `trig_eff_<sel>.csv` and the `BINNING DECISION` likelihood-ratio block. What the record says: the per-|y| SF is **0.9817 / 0.9963 / 0.9908 / 1.0059 / 1.0002 / 1.0069** for |y| in 0–0.4 … 2.0–2.4, a **2.4% spread** with errors ±1.1/0.7/0.9/0.4/0.6/0.5% (379–462 data events per bin), and the calibrated likelihood-ratio test rejects flatness at **p = 0.0006 (mt40) / 0.0004 (nom)**, ~3σ, in both selection variants — so the choice is deliberate, NOT a claim that the effect is absent. Full test table (pT flat p = 0.71/0.80; 3 |y| bins NOT enough p = 0.0043/0.0132; folding justified p = 0.39/0.14; no pT×η interaction p = 0.68/0.68; charge-inclusive p = 0.13) in the `BINNING DECISION` block of `correction/logs/trig_eff_mb_mu.log` — see the `trig_eff_mb.C` bullet. **What is NOT the justification:** the Clopper-Pearson pull χ² also printed in the `[SF]` block (11.5/11, p = 0.41 signed; 11.5/5, p = 0.04 folded). Toys show it is **~16× under-powered** — under a true flat null with the real denominators it gives ⟨χ²⟩ = 3.15 for 5 dof and rejects at 0.3% instead of 5%, because CP intervals over-cover badly as ε → 1 (2 of the 6 bins have ε_data = 1 exactly) — so it cannot support a flatness claim; both the macro and the header say so in place. **Cost of not applying it:** a pure rapidity SHAPE effect, the inclusive normalization being unchanged by construction — W templates would move by −1.5% (|y| < 0.4) to +1.0% (|y| > 2.0), |y|- and charge-symmetric, so σ_inclusive and the charge asymmetry are untouched and only dσ/dη and R_FB move (see the effect table below). **Switching back is one constant** (`kTrigBinning = kTrigPerAbsY`) + a muon MC re-skim + the downstream re-run; `kTrigPerY` gives the signed-y table. Per-bin errors in the record are **SYMMETRIZED** (the larger CP side mirrored; halving the interval width would understate the ε = 1 bins, whose central value sits ON the interval edge) so no Up template would equal its nominal; the inclusive value keeps its raw asymmetric CP errors (2521/2547 is far from 1, both sides informative). The Z uses the per-lepton path-fired efficiencies in the two-leg formula (inclusive data 0.9914, MC 0.9927). **ETA CONVENTION (verified by unit test 2026-09-15, after the user flagged the risk):** the analysis bins histograms in **y = −η_lab** (p-going = forward), but `MuonSF::W()`/`Z()` take the **RAW DETECTOR η** (`muEta`) — the flip is a BINNING convention and must never reach the POG lookup, because the POG tables are in signed detector η and are genuinely asymmetric (**ID × ISO differs by up to 1.06% between +η and −η**), so feeding `y` would mirror the ID/ISO correction with no error message — the same class of bug that hit `gen_xsec.C` in 2026-08. The trigger term is immune either way (inclusive it is one number; binned it negates its own argument and folds to |y| = |η|, so it is symmetric by construction — confirmed: TRIG(+η) ≡ TRIG(−η) exactly). Unit test: `W(35, η)` reproduces `id.Eval(η)` and never `id.Eval(−η)`, at η = ±0.2, ±1.0, ±2.2. The parameters are named `etaLab` so the convention is visible at the call site; `skim.C` passes `muEta->at(i)` directly. **Application:** W — `SF = ID(η) × [passIso ? ISO(η) : 1] × TRIG` of the LEADING muon after step 8 (TRIG = the one inclusive 0.9971; `TRIG(|y|)` under `kTrigPerAbsY`), multiplying EVERY fill (all ABCD regions have a TightID'd, matched leading muon; the anti-iso sideband has no measured SF → 1); Z — both legs ID × [iso ? ISO : 1], times the per-EVENT trigger factor `[1−Π(1−ε_data,i)]/[1−Π(1−ε_MC,i)]` (`skim_Zmm` requires the bit, no matching; = 1 ± 1e-4); `hMass_vipul` (no iso requirement) gets ID × TRIG only. **The Z→μμ iso cut was harmonized to the W's 0.15 the same day** (it was 0.20 = the POG Medium WP, for which the pp file has no TightID-denominator table; checked before deciding: Tight|TightID vs Tight|MediumID agree to ≤ 0.2% per bin, so the denominator would not have mattered), so ONE Tight-iso SF serves the W and both Z legs — and `skim_Zee` was already at `skim_Wel`'s 0.095. Side observation worth a pO tag-and-probe: 0.20 → 0.15 removes 6.2% of the Z→μμ DATA (388 → 364) but only 2.9% of the DY MC, i.e. the data iso efficiency is BELOW the (unembedded) MC's, while the pp iso SF moves the MC the other way (+1.6% per pair) — the pp-SF-in-pO caveat made concrete. Not covered: the reco/tracking efficiency (no SF supplied), the DY-veto second lepton, the 1e-4 Z-trigger correlation with the W bins. Data: always weight 1, nothing booked. **The LHE member weights are rescaled by the SF too** (`pOLhe::ScaleMemberWeights`), so member 0 ≡ the nominal fill weight and nPDF/qcdScale/alphaS sit on the SF-weighted nominal — every nuisance varies one thing around the same nominal, which is how Combine composes them. **Effect (muon re-skims 2026-09-14, backup of the pre-SF files in `skim/rootfile_pre_musf/`; the W data file is bit-identical, 216 histograms — the Z data file changed only through the iso cut):** ⟨SF⟩ over the nominal W selection W⁺ 0.98802 (ID 0.9835, ISO 1.0075, TRIG 0.9971) / W⁻ 0.98765 / DY-in-W 0.98700 / τ samples 0.990–0.991; Z→μμ DY 0.98280 (two legs: ID 0.968, ISO 1.016, TRIG 0.99998) → `Z_incl/signal` 372.1 (0.20, no SF) → 355.4 (0.15, SF) vs data 388 → 364; per lab bin the W template totals scale by 0.979–0.999 (the ID×ISO η dependence only, now that the trigger factor is flat; y3 0.979, y4 0.999), charge-symmetric (W⁺ and W⁻ agree to ≤ 0.06%). `correction/njet_WZ.C Zmm` re-run: 364 data events = the skim's `hMass` integral (event-identical at the new cut). **Effect of the per-|y| trigger SF, measured 2026-09-15 and REVERTED 2026-09-21 (it is what the inclusive choice gives up; backups `skim/rootfile_pre_etasf/` = before it went in, `skim/rootfile_pre_incltrig/` = before it came out):** ⟨SF⟩ 0.98802 ⇄ **0.98794** (⟨TRIG⟩ 0.99708 ⇄ 0.99699) — i.e. the INCLUSIVE normalization is unchanged by construction and the correction is **pure rapidity shape**: per lab bin the W templates move by **−1.54% (y5/y6, |y| < 0.4) to +0.99% (y0/y11, |y| > 2.0)** going per-|y|, and the reverting re-skim reproduced exactly the inverse (measured ×1.01563 at y5/y6, ×0.99022 at y0/y11, matching SF_incl/SF_bin to 5 decimals), exactly |y|-symmetric (y0 ≡ y11, y5 ≡ y6) and **identical for W⁺ and W⁻ to all printed digits** (so the charge asymmetry is untouched either way; the F/B ratio is not, because F/B pairs in y_CM are shifted by +0.3466 from the |y_lab|-symmetric SF). Inclusive totals ×1.00005 (W⁺) / ×1.00039 (W⁻) on the revert. Verified on the 2026-09-21 revert re-skim (`compare_reskim.sh rootfile_pre_incltrig`): the W and Z **data** files bit-identical (216 / 14 histograms), the electron files not touched, and **every histogram that stayed IDENTICAL in an MC file is empty** (521 of 1272 in a W file = the cross-charge `h_*_W∓_*` set and its twins; checked explicitly — 0 non-empty identical), i.e. every filled MC histogram moved, as it must. Downstream re-run on both occasions: `run_lhe_updown.sh Wmu Zmm` (a re-skim RECREATEs and wipes the 576/6 LHE Up/Down twins — always re-run it; `max|member0 − nominal| = 0` confirms the member weights picked up the new SF), `run_qcd_abcd.sh mu` (A40 μ⁺ 132.4 / μ⁻ 127.0 and κ 1.09 all unchanged — the ABCD is y-inclusive, so a |y|-shape-only SF cannot move it), `mtandmet`/`dileptonpeak` mu (**sidecars unchanged ⇒ the fork needs NO change**; `data_obs` 4998 and `qcd_abcd` 259.4 unchanged), `run_syst_shapes.sh met leppt_mt40` (LHE closure 1.6e-8 / 7.6e-8, member 0 ≡ `signal` exactly, one-sided bins 29 in leppt_mt40 / 84 in met), Z peak data 364 / `Z_incl/signal` 355.4 unchanged. **NOT yet refitted (lxplus).** **Uncertainties: three INDEPENDENT sources, ONE nuisance in the fit (user decision 2026-09-14: "one Up and one Down with all three combined").** Every fit-template histogram gets the per-SOURCE TH1D twins `<h>_muIDUp/Down`, `<h>_muIsoUp/Down`, `<h>_muTrigUp/Down` (`pOSF::kMuonSFSourceNames`, `BookSFTwins`/`FillSFTwins`) = the TOTAL event factor with that ONE source at ±1σ and the others nominal, filled with the SF-free weight — NOT area-normalized (the normalization IS the effect) — plus the COMBINED pair **`<h>_muSFUp/Down`** built at the end of the job (`FinalizeSFTwins`): per bin Up = nominal + √Σ_s(Up_s − nominal)², Down = nominal − √Σ_s(nominal − Down_s)² — the exact 1σ of a product of independent factors (ID and ISO are factorized POG T&P results, the trigger comes from our MB sample; adding linearly would assume full correlation). (Verified against the three-source inputs: the per-region `muSF` templates equal the quadrature of the three old ones to 2e-16; only templates that SUM two samples — the legacy `W_incl` dirs, `w` under the Z peak — differ from the quadrature of the summed shifts, at ≤ 1.4e-4 relative, because a sum of per-sample quadratures is not the quadrature of the sum; the near-proportional shifts make it negligible.) Only `muSF` (`pOSF::kMuonSFSystNames`) is carried into the Combine inputs and the cards, where it is treated exactly like nPDF/qcdScale/alphaS: one `shape` row = one Gaussian-constrained θ, templates at θ = ±1, pull + constraint reported, Asimov closure at 0; the per-source twins stay in the skim files as diagnostics (768 histograms per W MC file = 96 templates × 8, 8 per Z MC file). One coherent parameter loses nothing here: all three factors are near-flat normalization factors (ID/ISO = the XPOG systup/systdown totals, systematic-dominated at pT > 25; muTrig = the inclusive CP error). Sizes (inclusive, on Σw·SF): sources muID ±0.055%, muIso ±0.285%, **muTrig +0.20/−0.24%** (the inclusive CP error) → **muSF +0.36/−0.39% on the W templates**, ±0.52% on the Z peak (two legs — the Z trigger factor is 1 ± 1e-4 in every binning). Measured per lab region 2026-09-21, the `muSF` shift is essentially **flat in rapidity, +0.32% … +0.40%**, the residual variation being the ID×ISO η dependence, not the trigger. **This is one more consequence of going inclusive:** with the per-|y| trigger SF the same numbers were muTrig ±0.704% and muSF ±0.762%, running ±0.54% (y2/y9) to ±1.19% (y5/y6) across rapidity — the growth was the price of the COHERENT correlation model, since the six |y| errors are statistically INDEPENDENT (disjoint event samples) and moving them all to ±1σ together adds them linearly rather than in quadrature, ≈√6 too large on the inclusive normalization. One inclusive number has no such ambiguity: it is one measurement with one error, coherent by definition, and it contributes exactly ZERO to the rapidity-shape observables (dσ/dη, R_FB, charge asymmetry), where a coherent rescaling cancels. Decorrelating (only meaningful under `kTrigPerAbsY`) = set `pOSF::kMuTrigCorr = "perbin"` and ship `muTrig` separately; the fork then splits it with `nuisance edit rename`. `pOSF::kMuTrigCorr` and the sidecar directive `#! muTrig corr` stay dormant while the combined `muSF` is shipped (they only matter if the three sources are ever listed separately again, `kMuonSFSystNames = kMuonSFSourceNames`). Every job prints the `[SF]` block (inputs + provenance, the y-table of trigger SF and efficiencies, ⟨SF⟩ per source, Up/Down totals, lookup counters) — `skim/logs/skim_{Wmu,Zmm}_<sample>.log` is the record. `PO_MUON_SF=off` fills MC without SFs (regression checks; the log says so). **Missing inputs are FATAL**, never a silent SF = 1: a fresh checkout must run `correction/run_trig_eff_mb.sh mu` before the muon MC skim. Electrons: none yet — the same pattern (a `pOSF`-style header, twins named e.g. `eID/eIso/eTrig`, listed by the electron sidecars only) is the next step after the refit.
- [skim/compare_hists.C](skim/compare_hists.C) + [skim/compare_reskim.sh](skim/compare_reskim.sh) — **bit-identity gate for re-skims** (2026-09-07): every TH1-derived key of a backup file vs the new file (axes, all cells incl. under/overflow, Sumw2, entries, exact equality) → IDENTICAL/DIFF/MISSING + the NEW keys; `./compare_reskim.sh rootfile_pre_lhe` loops a backup dir. Use it whenever a skim change is supposed to leave the existing histograms untouched.
- [skim/mc_norm.h](skim/mc_norm.h) — single source of truth for the absolute MC→data scale `k_s = σ·L/N_gen` (`pONorm::MCScale`); consumes `ngen.root`. See "MC normalization" below.

Isolation studies (ROC curves, working-point tuning) are per-flavour and live in
`correction/` (moved out of `skim/`):
- [correction/isolation_mu_tight.C](correction/isolation_mu_tight.C) — muon (current: pure TightID Δβ scan)
- [correction/isolation_ele.C](correction/isolation_ele.C) — electron
- [correction/isolation.C](correction/isolation.C) — muon legacy multi-cone/multi-def scan (kept untouched)
- plotted by [correction/plot_iso_summary.C](correction/plot_iso_summary.C) / [correction/PlotsIsoROC.C](correction/PlotsIsoROC.C) / [correction/PlotIsoROC_ele.C](correction/PlotIsoROC_ele.C). See the `correction/` section below.

**relIso definition (since 2026-07-06):** every isolation cut in the pipeline uses
the **Δβ PU-corrected** PF relIso (MuonPOG convention),
`(ch + max(0, neu + pho − 0.5·PU)) / pT`, computed from the
`{mu,ele}PF{Ch,Neu,Pho,PU}Iso` branches — `pOSkim::RelIsoPF` in `skim_common.h`
is the single source (5-arg; null PU vector ⇒ uncorrected fallback + `[WARN]`).
Verified to reproduce the ntuplizer's `elePFRelIsoWithDBeta` exactly. Cut values
unchanged (μ 0.15 / e 0.095 — both re-confirmed as Youden-J optima under Δβ).

**Selection (8-step cutflow in the W macros)**:
1. ≥1 PF lepton with pT > 25 GeV
2. pO event-quality filters (auto-detected from branches). **Since 2026-08-18
   this is `pprimaryVertexFilter` ONLY** — `pclusterCompatibilityFilter` is
   commented out in `PassEventSelection_pO` (`skim_common.h`; the filter is
   problematic in pO and will not be used). Branch wiring in `skim.C` and
   `correction/njet_WZ.C` is kept, so re-enabling is a one-line uncomment;
   `njet_WZ.C` shares the helper and follows automatically. **Re-skim + local
   downstream re-run DONE same day (2026-08-18)** — step 2 now retains ~99.8%
   (data) / ~99.7% (MC). **The filter was cutting ~10% of DATA ONLY** (MC passed
   it at ~100%: MC Z-peak totals moved 372.1→372.2 / 252.5→252.6, i.e. nothing),
   so beyond statistics it was a ~10% data-only depletion the absolute MC
   templates never modeled ⇒ **all pre-2026-08-18 fitted r's and σ's are biased
   ~10% LOW**; expect σ up ~10% on the next fit. New counts (supersede the
   2026-08-12 ones everywhere below): **data W candidates μ 6055 (was 5480,
   +10.5%) / e 7948 (was 7186, +10.6%)**; m_T>40 μ 5050 / e 5001 (was
   4382/4465); Z-window data μμ 388 / ee 284 (was 356/~248). ABCD re-run
   (MET-plane, ⁺/⁻): μ T 0.2215/0.2100, iso-pass QCD 629.1+604.5 (+28%);
   e T 0.3891/0.3409, iso-pass QCD 2224.2+1991.0 (+15%); m_T-plane T
   μ 0.1767/0.1707, e 0.3442/0.3282. **The +28% is NOT "the sideband growing
   faster than 10%"** (that earlier reading was wrong): region D is 99.8% pure
   QCD and grew exactly with the data, but region B is 58%/53% EWK (μ, MET
   plane), so a **data-only** +10.5% against a FIXED r=1 MC lands entirely in
   the difference and inflates QCD_B by ~2×. Same arithmetic widened the μ
   MET↔m_T-plane T transport ~10%→~20% — the m_T-plane B is only 21%/20% EWK
   (m_T<30 excludes the Jacobian better than MET<30), so its T moved just +3%.
   NB this is also why the trigger-ΔR fix did NOT do this: that cut moves data
   AND MC together, the cluster-compat filter was data-only (MC passed ~100%).
   **DIAGNOSIS 2026-08-19: the widened transport is dominated by the PREFIT
   (r=1) EWK subtraction, NOT by iso–recoil correlation.** Scanning the EWK
   scale, the two planes agree exactly at r≈1.18 (μ⁺) / 1.22 (μ⁻) — where the
   post-filter-change fit is heading; for e the spread is nearly r-independent
   (−11.6%→−7.5% over r=1.0–1.3, B only 24%/22% EWK) and is more likely
   genuine. **So do NOT raise κ_μ on the 20% number — it double-counts the
   signal strength the fit itself determines; κ_μ=1.15 stands pending a
   decision.** Options: (a) keep it, state the r-sensitivity; (b) normalize the
   mt40 pT template with T_mT instead of T_MET — nearly r-immune and the
   natural transport for an m_T>40 template (μ QCD 164.0/155.8 → ~131/127, e
   729.6/690.8 → ~645/665); (c) iterate the ABCD at the fitted r (two lxplus
   round-trips). ~~Decide BEFORE the next fit~~ **RESOLVED 2026-08-23 by the
   in-fit ABCD migration (`QCD_MODE=abcd`), which supersedes BOTH (b) and (c):
   the m_T-plane counts become fit channels, so the normalization is T_mT-based
   by construction and the EWK subtraction is continuously at the fitted r (no
   iteration); the transport row never enters the κ. See the `qcd_abcd.C`
   bullet and "Downstream fit".** All 6 `combine_input_W*.root` + both `combine_input_Z.root`
   rewritten; njet_WZ re-run confirms the selections stay event-identical
   (data counts exactly 6055/7948/388/284). Benign: mtandmet's MET-vs-m_T
   integral echo fired on two e bins (1-event MET>120 overflows). Fits NOT yet
   re-run (lxplus).
3. Trigger: `HLT_OxyL1SingleMuOpen_v1` (muon) or `HLT_OxyL1SingleEG10_v1` (electron) — L1-seeded paths (skim.C:362, 1048; the older "HLT_PAL3Mu12/PAL3Ele12" in this doc was WRONG, corrected 2026-08-12)
4. DY veto: reject events with an OS dilepton pair, both legs ID'd + isolated
   (μ: pT>15, Tight, relIso<0.15; e: pT>10, eleMVAIdWP95, relIso<0.095), mll in (80,110)
5. ≥1 Tight ID lepton. **Electrons since 2026-09-14: "ID'd" = `eleMVAIdWP95`
   AND NOT in the ECAL barrel–endcap crack, 1.4442 < |η_SC| < 1.566 evaluated
   on the SUPERCLUSTER η (`eleSCEta`; the track `eleEta` moves by
   ~z_vtx/(R_ECAL·cosh η) ≈ 0.02, comparable to the crack half-width).**
   **The |η| < 2.4 acceptance, by contrast, is cut on the electron's own
   `eleEta`** (user decision 2026-09-27: the same fiducial variable as the
   muon |η| < 2.4 and the gen fiducial; `skim_Wel`/`skim_Zee` and the four
   replicating macros already do this). Near the edge η_SC − η ≈ v_z/335 cm
   (RMS 0.015); in W⁺ MC the two definitions swap 0.24%/0.26% of electrons,
   net 0.01%. The separate ID+iso SF study (`correction/idiso_sf_*`) cuts 2.4
   on η_SC — kept deliberately (user), since the effect is negligible there.
   Single source `pOSkim::InEcalGap` / `kEcalGapLo,Hi` + `BuildEleIDNoGap`
   (`skim_common.h`): `skim_Wel` builds `eleIDNoGap` per event and passes it
   wherever `eleMVAIdWP95` went (DY-veto legs, this gate, the leading-electron
   pick), `skim_Zee` vetoes both legs, and the replicating macros
   (`njet_WZ.C`, `charge_flip.C`, `trig_eff_mb.C`, `isolation_ele.C`) apply
   the same cut; `eleSCEta` missing ⇒ FATAL (a file without it would silently
   get a different selection). RECO level only, data and MC alike — the gen
   fiducial stays |η| < 2.4, common with the muon channel (shared r's need one
   fiducial; the crack = 0.122 of each 0.4-wide |η| ∈ [1.2,1.6] bin is an
   acceptance hole filled by the MC's intra-bin rapidity shape; a gen-level
   veto would move the same extrapolation into the μ/e ratio and change
   nothing for the combined σ). Why: the crack efficiency (A·ε ≈ 0.62 there)
   is unmeasurable — tag-and-probe excludes the crack too. **Effect (electron
   re-skim 2026-09-14, backup `skim/rootfile_pre_gap/`): data W candidates
   7948 → 7726 (−2.8%), m_T>40 5001 → 4864; ONLY the two crack bins move —
   y2 [−1.6,−1.2] −22.8%, y9 [+1.2,+1.6] −18.3% (data), −18.5% both charges in
   W MC (charge-symmetric = geometric, as it must be), every other bin
   bit-identical; MC inclusive −2.5%; Z→ee data 284 → 273 (−3.9%), DY MC
   −4.7% (two legs); ABCD A40 e⁺ 652.1 → 627.5 ± 29.5, e⁻ 670.7 → 631.2 ± 29.8;
   a track-η estimate had predicted −3.8% inclusive (the SC η is what counts).
   Downstream re-run the same day: `lhe_updown` Wel/Zee, `qcd_abcd ele`
   (in-fit reduced κ_e now **1.12**, was 1.15 — the fork's
   `QCD_ABCD_LNN_ELE` default still says 1.15, decide before the refit),
   `mtandmet`/`dileptonpeak` ele (sidecars unchanged ⇒ no fork change),
   `syst_shapes all` (closure 1e-8, member 0 ≡ `signal`), `njet_WZ` Wel/Zee
   = 7726 / 4864 / 273 data events (event-identical with the skim; its MC
   counts 443216 / 414099 also equal `charge_flip`'s — three independent
   implementations agree), `charge_flip ele` f = 2.50×10⁻³, ΔA_e = −0.00061
   (unchanged), `trig_eff_mb ele` SF 0.9895 ± 0.0032 plain / 0.9975 ± 0.0037
   with m_T>40 (unchanged), `isolation ele` re-run: MVAIdWP95 AUC_MET 0.9105,
   J_MET optimum still exactly 0.095, ε_sig 0.917 / ε_QCD 0.187 at the cut
   (was 0.909 / 0.095 / 0.910 / 0.184 — the working point stands). Fit NOT
   yet re-run.**
6. Leading lepton pT > 25 GeV
7. Leading lepton PF relIso < 0.15
8. Trigger matching on leading lepton (ΔR < `trigMatchDR` = **0.4 since 2026-08-12**, was 0.1; skim.C:283/963, against the objects in `hltobject/HLT_Oxy*` = **L1-granularity candidates**). **Only the two W channels do this** — `skim_Zmm`/`skim_Zee` wire the trigger-object branches but never match on them.

**FIXED 2026-08-12: `trigMatchDR` 0.1 → 0.4 in both W channels + `correction/njet_WZ.C:194` (which replicates the selection), W re-skim and full downstream re-run DONE.** Post-fix cutflow step 7→8: μ⁺/μ⁻ **99.99%/99.99%** (was 93.88%/98.82%), e⁺/e⁻ **99.89%/99.88%** (was 97.79%/98.89%) — charge-symmetric in both flavours. **Data W candidates: μ 5203 → 5480 (+5.3%), e 7065 → 7186 (+1.7%)**; `njet_WZ` data counts still equal the skim's N[8] exactly (5480 / 7186), so the two independent selection implementations remain event-identical. ABCD moved with it: iso-pass QCD μ⁺/μ⁻ 413.7/453.8 → **491.6/473.2** (the artificial μ⁺ deficit is gone from the sideband too — ratio 0.91 → 1.04), e⁺/e⁻ 1887.2/1693.9 → **1937.9/1723.1**; MET-plane T μ 0.1895/0.1795, e 0.3704/0.3204. `ptmt_scan` baseline S/√(S+B) μ 53.1 → **54.40**, e 39.6 → **39.85** (same conclusions: optimum still at the pT floor with m_T > 40; e best (20, 50) = 50.27) — **both superseded 2026-08-20 by the post-filter-change re-run: 52.83 / 38.50, see the `ptmt_scan.C` bullet**. **The defect it removes:** The saved trigger object's φ is at the muon station, `muPhi` is the vertex-direction φ, and the difference decomposes **exactly** (measured in barrel pT slices) as

&nbsp;&nbsp;&nbsp;&nbsp;`Δφ(q, pT) = c + q·k/pT`, with **k = 2.55 rad·GeV** (solenoid bending out to R ≈ 4.5 m = the muon station; the fitted k is constant to 3% over pT 25–200) and a **charge-independent offset c = +0.0060 rad** (constant to 5% over the same range).

Both terms are understood. `k`: `muPhi` is the momentum direction **at the vertex**, the object φ is the muon's azimuth **at the muon station**, and in the 3.8 T solenoid the track's azimuth changes continuously between the two — the bend *direction* is set by F = qv×B, hence the ±q. `c`: the saved objects are **L1-granularity candidates** (pT quantized in 0.5 GeV steps; φ takes exactly **576 distinct values** = the μGMT grid 2π/576 = 10.9 mrad; η likewise on a ~11 mrad grid) — reporting the cell's lower edge instead of its centre biases φ by half a cell = **5.45 mrad**, matching the measured c = 6.0 mrad.

The bending term alone is charge-*antisymmetric* and would cost both charges equally — **`c` is what breaks the symmetry**: it adds to the μ⁺ bending and cancels part of the μ⁻ one, so |Δφ⁺| − |Δφ⁻| = 2c ≈ **+0.012 rad at every pT**. That 12 mrad is decisive only because ΔR < 0.1 slices straight through the bending scale at low pT (Δφ ≈ 0.09–0.10 rad at pT 25–30 in the barrel): matching efficiency by pT slice, μ⁺/μ⁻ = **0.478/0.905** (25–30), 0.956/0.988 (30–35), ≥0.985/0.993 above 35. The lepton pT spectra are NOT the cause (⟨1/pT⟩ differs by only 2% between charges). Net: ~10% charge-dependent barrel-only inefficiency, **entirely concentrated at pT 25–30** — i.e. right on the Jacobian turn-on, the worst place for the lepton-pT discriminant. Endcap is unaffected (Δφ ≈ 0.04, both charges ≥0.998). Cutflow confirms it inclusively (step 7→8: μ⁺ −5.72%, μ⁻ −1.11%; every other step agrees between charges to <0.2%), and a direct per-cut MC scan shows reco+TightID+iso is flat in η and charge-symmetric at 0.95, so step 8 is the sole source. It is what makes the μ⁺ (A·ε)(η) dip to 0.81 in the barrel while μ⁻ stays at 0.92 (endcap 0.96 both) in `xsec_fiducial_diag`. **Electrons are much less affected** (step 8: e⁺ −1.56%, e⁻ −0.79%) because L1 EG objects are ECAL-position-based, like the offline electron φ. Applied to data and MC alike, so it largely cancels in r — **but it rides directly on the charge-asymmetry observable and discards ~5.7% of μ⁺ events in a stats-limited measurement.** **The applied fix: `trigMatchDR` = 0.4** (chosen 2026-08-12; max |Δφ| is ~0.10 at threshold, so 0.4 clears it with margin). Measured in MC on the barrel pT 25–30 slice, matching efficiency μ⁺/μ⁻ goes 0.479/0.906 → **1.000/1.000** (ΔR<0.3 already saturates; |Δη|<0.1-only gives 0.996/0.996). NB correcting only the half-cell offset would fix the *charge* bias (diff −0.427 → −0.045) but still lose ~28% of BOTH charges at threshold — the window had to widen regardless. NB the +0.006 rad offset `c` is itself worth understanding — if it is a geometry/propagation artifact it need not be identical in data and MC, which would turn this into a direct charge-asymmetry bias rather than a cancelling one.

**Output**: per-sample ROOT files in `skim/rootfile/` (W: `{channel}_pO_PFMet_{sample}_hist.root`;
Z: `ZToMuMu_pO2025_*` / `ZToEE_pO2025_*`). W files contain 12 rapidity-binned MT
histograms (`h_mt_Wp_y0..y11`, `h_mt_Wm_y0..y11`), MET distributions, isolation
histos, **plus leading-lepton kinematics `h_lepPt/h_lepEta/h_lepPhi`**
(2026-07-29: leading lepton after the FULL selection incl. iso + trigger match,
charges combined; deliberately the same names as the Z-channel histos so
`correction/dataMC_kinematics.C` reads W and Z files uniformly). **Both `skim_Wmu` and `skim_Wel`** additionally store the ABCD inputs for
the QCD/low-MET background: per-charge 2D histos `h_iso_met_{mu,ele}{Plus,Minus}`
(relIso × PF MET) and `h_iso_mt_{mu,ele}{Plus,Minus}` (relIso × m_T), filled for
every event passing the full W selection **except** the isolation cut, so the
relIso axis spans iso-pass (signal) and the anti-iso sideband. Binning: relIso
0–1.0 in **0.005** steps (`kNIsoAB = 200`); MET 0–120/2 GeV; m_T 0–200/2.5 GeV.
Regions are chosen downstream by projecting these (no re-skim) in
`correction/qcd_abcd.C`. **The 0.005 width is load-bearing** (fixed 2026-07-30):
every boundary in use must be a bin EDGE — the electron cut 0.095 (= 19×0.005),
the muon 0.15, the anti-iso edges 0.20/0.30/0.60/0.65. With the previous 0.01
bins, 0.095 fell mid-bin and `qcd_abcd.C`'s `FindBin(isoCut−eps)` silently
integrated relIso < 0.10, leaking 170 events that fail `skim_Wel`'s own iso cut
into the ABCD iso-pass region and inflating T + every electron QCD template by
~4% (muon was always exact). Any new cut value applied by projection must land
on an edge.
**2026-07-29 additions (both W channels):** a third ABCD plane
`h_iso_pt_{mu,ele}{Plus,Minus}` (relIso × leading-lepton pT, 0–100/2 GeV) —
supplies the anti-iso lepton-pT *shape* for the QCD pT template (NOT used for
the 2×2 counting; relIso is pT-correlated) — and the rapidity-binned
per-charge lepton-pT histos `h_leppt_W{p,m}_y{0..11}` (+`_FB`), 2 GeV bins,
the pT twins of `h_met_*`/`h_mt_*` feeding the display-only lepton-pT stacks
in `mtandmet.C`. **2026-07-30:** joint per-charge (pT × m_T) 2Ds for the
cut-pair scan (`correction/ptmt_scan.C`): `h_pt_mt_*` (iso-pass) and
`h_pt_mt_antiiso_*` (the anti-iso sideband, relIso in the qcd_abcd windows,
with its actual MET/m_T — the scan's QCD model); pT 1 GeV bins, m_T 2.5 GeV.
**2026-08-04: these four scan planes are filled down to leading-lepton
pT > 20** (`kPtScanFloor` in `skim_common.h`, = the scan grid's lower edge) so
the scan can probe LOWERING the nominal pT cut. They are the ONLY histograms
with the relaxed floor — every other histogram AND the printed cutflow keep
the nominal pT > 25 selection (`passPtNominal` gates the fills; `hasPF25`
gates N[1..5] so the cutflow counts are unchanged; verified bin-identical on
the nominal set after the re-skim).
**The `_mt40` set (2026-07-30) = the LEPTON-pT-DISCRIMINANT selection**
`pT > 25 && m_T > 40` (local `mtCutForPtDisc` in both W functions; the m_T cut does
the QCD suppression that the MET *shape* does in the MET fit). Twins of the
existing histos with that cut applied: `h_leppt_mt40_W{p,m}_y{0..11}`(+`_FB`),
`h_lepPt_mt40`/`h_lepEta_mt40`/`h_lepPhi_mt40`, and the ABCD plane
`h_iso_pt_mt40_{mu,ele}{Plus,Minus}`. The nominal (no-m_T-cut) histograms, the
MET/m_T histos and the cutflow are untouched — the W skim still applies no
MET/m_T cut. Z files contain dilepton mass (`hMass`, `hMass_extended`, `hMass_vipul`)
plus kinematics of the dilepton system and its leptons (`h_Zpt/h_Zeta/h_Zy/h_Zphi`,
`h_lepPt/h_lepEta/h_lepPhi` — both legs; `h_Zy` is lab-frame rapidity), filled for
iso-selected OS pairs in the Z peak [60,120] GeV (both legs Tight/MVA-ID'd with
relIso < 0.15 μ / 0.095 e = the W cuts — the muon value was 0.20 until
2026-09-14; consumed by
`correction/dataMC_kinematics.C`). **`skim_Zmm` additionally** stores the hadronic
recoil for the MET correction (Sec 6 of AN2017_058): `h_uPar`/`h_uPerp` (1D
inclusive) and `h_uPar_qT`/`h_uPerp_qT` (2D, recoil component vs q_T), where
`u = −MET − q_T` and `q_T` is the dimuon pT. u-axis = AN's 2 GeV/c binning; q_T
axis fine so the fit binning can be chosen later by projecting (no re-skim). PF
tree is **soft-gated** (loud `[WARN]`, empty recoil histos if `pftree` absent).
**`skim_Zee` stores the same four recoil histos** (`u = −MET − q_T`, `q_T` = the
dielectron pT), wired identically (soft-gated PF tree); consumed by
`correction/recoil_raw_ele.C`.

**LHE member twins (2026-09-07, MC files only):** every fit-template histogram
— W: `h_met_W{p,m}_y*`(+`_FB`), `h_leppt_mt40_W{p,m}_y*`(+`_FB`) = 96 per
file; Z: `hMass` — gets three TH2D twins `<h>_epps21` (107 members) /
`<h>_scale` (9) / `<h>_alphas` (5), see `lhe_index.h`; plus, after
`run_lhe_updown.sh`, the TH1D pairs `<h>_nPDFUp/Down`, `<h>_qcdScaleUp/Down`
and `<h>_alphaSUp/Down`. W MC files grow 0.3 → 5 MB;
a W MC job takes ~30 s. Data files are unchanged (no `ttbar_w`, nothing
booked). The job-log line `[INFO] LHE member twins booked: 288 (107+9+5
members per template); events without usable ttbar_w: 0` is the record.
`compare_reskim.sh <backup>` is the bit-identity gate (2026-09-07 MC re-skim:
every pre-existing histogram IDENTICAL in all 24 MC files).

**Muon-SF twins (2026-09-14, muon MC files only):** the same fit-template set
additionally carries the per-source TH1D pairs `<h>_muIDUp/Down`,
`<h>_muIsoUp/Down`, `<h>_muTrigUp/Down` (one source at ±1σ, the others
nominal; diagnostics) and the COMBINED `<h>_muSFUp/Down` (the three added in
quadrature per bin = the one nuisance the fit sees), all written by the skim
itself (`skim/muon_sf.h`, NOT area-normalized) — 768 per W MC file (96 × 8),
8 per Z MC file; they survive `run_lhe_updown.sh` (UPDATE mode) and are wiped
by a re-skim like everything else. The job-log `[SF]` block is the record. The
2026-09-14 muon re-skim: data files IDENTICAL (`compare_reskim.sh
rootfile_pre_musf`), every MC histogram filled after step 8 DIFF by the SF.

Per-job logs in `skim/logs/`, cutflow text in `skim/output/`.

### Stage 2 — `plotting/`

Data/MC overlay, background estimation, plots.

- **[plotting/run_combine_inputs.sh](plotting/run_combine_inputs.sh) (2026-09-21) — the logging wrapper for the two macros that write the fit inputs.** `./run_combine_inputs.sh [mu|ele|both]` (default both) runs `mtandmet.C+(isElec)` and `dileptonpeak.C+(isElec)`, pre-building both once so two flavours cannot race on the ACLiC artifacts, and teeing to `plotting/logs/{mtandmet,dileptonpeak}_<chan>.log` = THE record. It echoes the headlines: the ABCD QCD totals, the in-fit `A0 = B0·C40/D0` with its rescale factor off the T_MET-normalized total, the per-y QCD split, and the per-file region + shape-systematic counts. Benign noise is filtered on purpose — the macros print `ERROR check plot N` as debug output (not a failure), and the `[WARN] … is empty -> flooring central bin` lines (empty `wtau` under the Z peak) are counted rather than listed, so any OTHER warning stands out. **Run `correction/run_qcd_abcd.sh` first** (mtandmet embeds its templates + `abcd_counts_*`), and re-run this after ANY re-skim. Verified idempotent 2026-09-21: re-running reproduced every integral exactly (μ `W_incl` data 4998 / signal 3678.9 / qcd 315.9 / qcd_abcd 254.7; e 4766 / 3035.5 / 1344.7 / 1259.9; Z 364 / 339.2 and 273 / 240.8) — the file checksums differ only through ROOT's embedded write timestamp.
- [plotting/mtandmet.C](plotting/mtandmet.C) — MT and MET data/MC plots. **MET stacks (muon AND electron) include the data-driven ABCD QCD** from `correction/rootfile/qcd_abcd_{mu,ele}.root` (run `correction/qcd_abcd.C[+(true)]` first): added to the inclusive MET stack (rigorous) and split across the per-y bins by the per-bin low-MET excess (data−EWK at MET<30) proxy, keeping the inclusive template shape (per-y sums back to inclusive). **Also writes the structured Combine input `plots[/Elec]/combine_input_W.root`** — one TDirectory per fit region (`Wp_lab_y3/`, `Wm_fb_y7/`, `Wp_incl/`, `W_incl/`, …) each with the 6 absolute templates `data_obs/signal/z/ztau/wtau/qcd` (consumed by the fork's `run_pO_fits.sh` — see "Downstream fit"). Channel auto-selected by `isElec`. **Info box of the MET + lepton-pT stacks (2026-09-24, user request): THREE counts** — `Passing Events` (data), `W signal MC`, and `Total MC+QCD` = the whole stack as drawn (the same MC total the pull pad uses; `Total MC` when no QCD template is stacked, the label derived from the legend names by `totalLine`), all over bins 1..N — so data vs total separates a normalization offset from a shape mismatch. The box band sits half a line lower (`psCounts`) so the first line stays clear of the header; the m_T stacks keep the data count only. Inclusive data/total: μ 1.10 (MET) / 1.17 (lepton pT), e 1.02 / 1.05. NB the lepton-pT stacks draw the DISPLAY QCD (`qcd_pt_mt40`, T_MET-normalized), not the fit's `qcd_abcd` prefit (A0): inclusive it is 61 (μ) / 85 (e) events higher, so the abcd-mode prefit total is 1.4% / 1.9% lower than the box says (data/total 1.19 μ / 1.07 e). Combine inputs regenerated bit-identical with the change (`compare_hists.C`, DIFF 0 in all six). **Lepton-pT stacks (2026-07-29):** the same absolute data-vs-stack comparisons in leading-lepton pT, driven by a `varStem[]`/`varTag[]` table that has held **ONE entry since 2026-09-21** (`kVarMt40`): **`plots[/Elec]/leppt_mt40/` = the pT-DISCRIMINANT selection pT>25 && m_T>40** (`h_leppt_mt40_*` + `qcd_pt_mt40_*`, 2026-07-30), per-y `leppt_mt40_W{p,m}_y*` + `_FB` + inclusive, per-y QCD split by the same low-MET-excess weights (the m_T cut doesn't move QCD in y). **RETIRED 2026-09-21: the second entry `kVarNom`** = the plain W selection (`h_leppt_*` + `qcd_pt_*` → `plots[/Elec]/leppt/` + `combine_input_W_leppt.root`). Rationale: the discriminant was dropped 2026-08-16, `pO_fit_out_leppt/` never had a simfit, and — the decisive part — `attachSysts` was applied only to `kVarMt40` because the skim stores LHE/SF twins for `h_met_*` and `h_leppt_mt40_*` only, so that input's sidecar listed **0 systematics** and the current card generator could not fit it. The writer now `Unlink`s any pre-retirement `combine_input_W_leppt.root[_systs.txt]` on every run (same guard the live inputs use — verified with planted decoys), `disc_variants.h` and `run_observables.sh` reject the tag by name, and the four surviving Combine inputs were re-generated **bit-identical** (`compare_hists.C`: 7178 histograms, DIFF 0). `qcdPtPlusBase`/`qcdPtMinusBase` are still loaded — they feed the `[QCD] m_T>40 lepton-pT QCD: … (vs no-m_T-cut …)` line, the measure of what the m_T cut removes. To restore: put the second entry back in the variant table + the `leppt` rows in the two consumers. NB `skim.C` still FILLS `h_leppt_W{p,m}_y*` (48 per W file); now unread, but dropping them needs a full re-skim so they were left. **`combine_input_W.root` itself stays PF-MET-only**, and since 2026-07-30 the same writer also emits the lepton-pT variant as a **separate file** — `combine_input_W_leppt_mt40.root` (pT>25 && m_T>40) — same 51 regions × 6 templates, pT binning 50/0–100 (**the fit sees exactly these 2 GeV bins, populated from 24 = the edge enclosing the 25 GeV cut; the overflow pT > 100 is DROPPED, data and MC alike — measured 2026-09-08 over the 24 lab regions: μ data 52/5050 = 1.0%, e 98/5001 = 2.0%, vs signal MC 0.7%/0.8% and QCD 4.8/34.3 events, i.e. the e tail excess continues past the window; NB Combine's saved `shapes_*` carry 60 unit bins because the common observable spans the widest channel, the 60-bin Z mass — bins 51–60 are outside the W templates, floored at 1e-9×integral**), consumed by the fork's `run_pO_fits.sh --disc leppt_mt40` (out-tree `pO_fit_out_leppt_mt40/`; `sync_lxplus.sh` uploads the variant when present and downloads all out-trees). (The 2026-07-30 "qcd_norm FREE for all three" decision was superseded 2026-08-17: QCD is now lnN-constrained at the ABCD prediction in the simfit — see "Downstream fit".) **In-fit ABCD inputs (2026-08-23, the `_leppt_mt40` file ONLY):** every SR dir gains a **7th template `qcd_abcd`** — same anti-iso m_T>40 pT shape and per-y split, total renormalized to A0 = qcdB0·qcdC0/qcdD0 read from `abcd_counts_*` (Σ_y = A0 exactly for lab AND fb since the weights sum to 1; staleness guard prints template-total/T_met vs C40, warns >2%; labels read by looping `GetBinLabel` — `TAxis::FindBin(label)` would APPEND on a miss) — and the file gains **6 CR dirs** `{Wp,Wm}_CR{B,C,D}` of 1-bin templates: `data_obs` (raw count), `qcd` (prefit EWK-subtracted count = the free CR scale's ×1 anchor), CRB `z`/`ztau` (swept into r_Z by the fork's map) plus BOTH W variants — per-y `w_lab_y0..11`/`w_fb_y0..11` (wB split by the SR per-y signal+wtau fractions of the matching binning → mapped to the r POIs) AND the frozen `wfix` — the card generator picks via `QCD_WCR`; CRC/CRD one frozen `ewk` (≤2.5%). The `met` file and all six existing template names are UNTOUCHED (verified: old `qcd` totals byte-identical, so `QCD_MODE=lnN|free` regress cleanly). **Composition after m_T>40** (data / W sig / DY+τ / QCD): μ 4382 / 3704 (84.5%) / 272 (6.2%) / 226 (5.2%), closure 96%; e 4465 / 3080 (69.0%) / 151 (3.4%) / 1202 (26.9%), closure 99.3%. (Before the cut: W purity 73.5% μ / 45.0% e.) (For electrons QCD is the *dominant* low-MET component — see `qcd_abcd.C`.)
- [plotting/mtandmet_overlay.C](plotting/mtandmet_overlay.C) — MC stack overlays
- [plotting/dileptonpeak.C](plotting/dileptonpeak.C) — Z peak plots
- (QCD sideband fit and isolation ROC curves moved to `correction/` — see below)
- [plotting/observables.C](plotting/observables.C) — final charge-asym + F/B plots of the fits, **built from the FIDUCIAL yields r × σ_gen, never from the raw fitted counts (user rule 2026-09-22 — see the `analysis/fiducial_yields.C` bullet)**, with **all four** nPDF theory bands (EPPS21/nCTEQ15HQ/nNNPDF3.0/TUJU21nlo). **Simfit-aware since 2026-08-04:** `observables_comb(disc)` = the PRIMARY μ+e-combined plots of the grand fit (reads `{charge_asym,FBratio}_fid_comb_<disc>.root` + `fidyields_comb_<disc>.root` for the W⁺/W⁻ weights of the summed-R_FB theory band — FIDUCIAL since 2026-09-22; until then the count-based `*_fit_comb_<disc>.root` made from `comb_fitted_yields.root` — label "W → l ν / μ + e simfit, fiducial") → `plots/comb/{charge_asym,FBratio}/<disc>/`; internals in `observables_run(chan, lepSym, outBase, …)`. **`observables_flav(disc, withComb=false)` (2026-09-22) = the μ-ONLY vs e-ONLY overlay** of the per-flavour fits (fork mode `flavfit`): A_ch and R_FB (sum, W⁺, W⁻), each fit drawn exactly like the comb one (point + stat bar + syst TBox) side by side in every bin (μ blue circle left, e red square right; `withComb` adds the grand fit in the middle, black diamond — conventions in `plotting/fit_variants.h`) via the new `plotting_helper.C::SaveNiceGraph_Overlay`, → `plots/flavfit/{charge_asym,FBratio}/<disc>/` + `mu_vs_e_chi2.csv` and a console χ² (Σ(x_μ−x_e)²/(σ²_μ+σ²_e) against the STAT errors — the two samples are independent — and, as a conservative bound, the total ones). **It reads the FIDUCIAL graphs `{charge_asym,FBratio}_fid_<fit>_<disc>.root`, NOT the count-based ones — see the `analysis/fiducial_yields.C` bullet: count-based R_FB of μ and e differ by up to 60% on identical physics.** `observables(disc)` (the file-named entry, since 2026-09-22) runs comb + flav. **REMOVED 2026-09-22 with the legacy per-flavour fits:** `observables(isElec, disc)` (per-channel plots) and `observables_overlay(disc)` (the old merged μ+e view, total bars only). **Discriminant-aware since 2026-08-03:** `disc` = `met`(default)|`leppt_mt40` (the third tag `leppt` was retired 2026-09-21 and is now rejected by name) drives BOTH the default inputs (the disc-tagged fiducial graph + yield files) and the per-disc output folders `plots/{comb,flavfit}/{charge_asym,FBratio}/<disc>/`, so variants never overwrite each other; the tag is stamped on every plot as a 4th header line ("PF MET fit" etc.). Mapping single-sourced in [plotting/disc_variants.h](plotting/disc_variants.h) (`pODisc::Spec`, unknown tag ⇒ hard error; `pODisc::GraphFile(stem, fit, disc)` = `../skim/rootfile/<stem>_fid_<fit>_<disc>.root`, no fallback — the count-based `*_fit_*` graph files, incl. the pre-2026-08-03 untagged ones, are never read). Reads `RpO_rootfile/RpO_FB_graphs.root` for theory (optional). The yields file reads `h_yield_*` first (`h_mt_*` legacy fallback). Whole Module-5 chain in one command: [analysis/run_observables.sh](analysis/run_observables.sh). The old "Projection with Electrons" pseudo-band is gone (real e overlaid instead). **STAT BAR + SYST BOX (2026-09-15):** picks up the `_stat` twins `g_chargeAsym_stat` / `g_RFB_{sum,Wp,Wm}_stat` by name and passes them as the trailing `gStat` of `SaveNiceGraph` / `SaveNiceGraph_ErrorBand`, which then draw statistical bars with the systematic as a `TBox` (see the `plotting_helper.C` bullet). Twins absent → one loud `[WARN]` naming the re-extraction command and a single total-error bar, unchanged.
- [plotting/xsec_fiducial.C](plotting/xsec_fiducial.C) — fiducial W cross sections of the simultaneous fits, σ_meas,i = r_i × σ_gen-fid,i: `xsec_fiducial_comb(disc)` (the grand fit), `xsec_fiducial_diag(disc)` (its A×ε diagnostic), **`xsec_fiducial_flav(disc, withComb=true, withCombDiff=false)` (2026-09-22) = the μ-ONLY vs e-ONLY overlay**, and `xsec_fiducial(disc)` = all three (the file-named entry). `xsec_fiducial_flav` loads each fit through the SAME `loadCombIngredients` + `sumWithCov` as the comb (per fit: its `<tag>_W_yields.csv` + `<tag>_fitted_yields.root` via `plotting/fit_variants.h`) with the SAME pooled σ_gen for every fit, so all results sit in ONE fiducial volume and **σ_e/σ_μ = r_e/r_μ exactly** (the ~1–2% bare-lepton FSR difference between flavours lives in each flavour's MC acceptance): `W_fiducial` (W⁺/W⁻/W, μ | comb | e side by side with stat bars + syst boxes, the r = 1 gen expectation as a dashed segment, the e/μ ratios with STAT errors printed under the legend — lumi cancels in them, the flavour-specific syst do not enter), `W_dsigma_deta_{Wp,Wm,W}` (per charge, μ vs e; the grand fit with `withCombDiff`), `xsec_flavfit.csv` (the `xsec_comb.csv` columns + a `fit` column) and `xsec_flavfit_ratio.csv` (e/μ per bin + inclusive, stat and total errors); console: per-fit σ with stat/syst, the e/μ ratios in σ from 1, and the per-bin μ-vs-e χ² of dσ/dη against stat and total errors. η_CM positions from the fiducial `charge_asym_fid_<fit>_<disc>.root`. → `plots/flavfit/xsec/<disc>/`. **REMOVED 2026-09-22 with the legacy per-flavour fits:** the N_fit/L `xsec_fiducial(disc, muCsv, eleCsv)` (it read the legacy `<chan>_summary.csv`) and `xsec_fiducial_diff(isElec, disc)`, with their CSV readers `readYield`/`readBinYields`. Run by `run_observables.sh`. **`xsec_fiducial_comb(disc)` — THE r×σ_gen MEASUREMENT since 2026-08-12** (2026-08-05..12 it was the N/L + prefit/gen-overlay view): **σ_meas,i = r_i × σ_gen-fid,i, PER LEPTON FLAVOUR** (μ/e-shared r's ⇒ σ(W→ℓν) under lepton universality; the old "×2 = μ+e summed" display convention is gone). Reads per-bin `r/rErr` from `comb_W_yields.csv` and σ_gen-fid from `gen_xsec.root`; inclusive W⁺/W⁻/W = Σᵢr_iσ_i with the full r-covariance recovered as cov_Y_ij/(S_i·S_j) from `h_cov_yield` × the `signal_prefit` column (diagonal-rErr fallback + WARN; gen MC-stat ~0.3%/bin neglected vs ~6% fit errors). Plots show measured points vs the gen-fiducial r=1 expectation (green diamonds / dashed lines — measured÷dashed per bin IS r_i; the prefit-open-marker overlay was dropped, that role is absorbed). Also writes the σ table `plots/comb/xsec/<disc>/xsec_comb.csv` (per-bin + inclusive rows; incl `r` col = gen-weighted mean r). **Uncertainty bars (2026-09-14):** the fitted `rErr` has been the TOTAL profiled error ever since the 2026-08-17 nuisances (the old "stat. unc. only" label was wrong from then on); the fork's extractor supplies the statistical part from the nominal fit's covariance matrix conditioned on the constrained nuisances (= the frozen-nuisance error, no refit; CSV column 19 `rErr_stat`, `h_cov_yield[_FB]_stat`, `simfit_<B>_stat` summary rows — see "STAT/SYST SPLIT" under "Downstream fit"), so `xsec_comb.csv` gains `sigma_meas_stat_nb` / `sigma_meas_syst_nb` (syst = √(total² − stat²)); without the stat column a single total bar is drawn and labelled as such. **DRAWING SETTLED 2026-09-15 (user decision): point + STATISTICAL error bar + the systematic as a `TBox` per point** (the 2026-09-14 inner-thick/outer-thin double bar is gone) — `plotting_helper.C::MakeSystBoxes` is the single source, see the `plotting_helper.C` bullet. `W_fiducial` sets `ps.systBoxHalfWidthAbs = 0.13` (its three points sit at discrete x with no x error); `W_dsigma_deta_comb` draws one box set per series in the series colour (W black / W⁺ azure / W⁻ red) — they share x but sit at different y, so they never overlap. **MARKERS (2026-09-15, user: "the data point looks a bit shifted" + "avoid triangles"):** the W⁺ square looked off its error bar in the legend — measured, and it is a ROOT RASTERIZATION offset, not a shift in the data: style 21 at size 1.3 renders a 10-px block whose centre lands **1.0 px left and 1.5 px high** of the width-2 error bar (an even-pixel glyph cannot straddle the bar). The legend makes it conspicuous because the `"lep"` option runs a long horizontal line through the glyph. Measured offsets (glyph centre − bar centre, px): circle 20 = (+0.5, 0.0) at EVERY size; square 21 = (−1.0, −1.5) at 1.3, (0.0, −0.5) at **1.5**; diamond 33 = (−0.5, 0.0) at **2.0**; triangle 22 = (0.0, +0.5…+1.0) and open triangle 26 = (+0.5, +1.0…+1.5). So **triangles are banned from these plots** — doubly so because a triangle's visual centre sits ~1/3 of the height above its base, which makes the eye read the point low whatever the rasterizer does. The convention is now single-sourced at the top of `xsec_fiducial.C` as `kMkW/kMkWp/kMkWm` (+ sizes, + the open `kMkWOpen/kMkWmOpen`): **W = filled circle 20 @ 1.3, W⁺ = filled square 21 @ 1.5, W⁻ = filled diamond 33 @ 2.0** (all point-symmetric, widths 11/12/11 px so they read at the same weight), and the `AxEps_vs_eta` table went `{{20,22},{24,26}}` → `{{20,33},{24,27}}`. Verified by re-measuring the rendered PNG: the blue square is now EXACTLY on the bar in x and y, the red diamond within 0.5 px. Every other square in the macro was moved off the uncentred 1.3/1.4 to 1.5. **Re-measure before changing any of these sizes** — the good values are not monotonic in size (21 is centred at 1.5 and 2.0 but off at 1.6–1.8).
**The lumi 3% stays INSIDE the syst box** (user decision 2026-09-15): on σ it is 3.19 nb of the 3.25 nb W-inclusive systematic, so the box is essentially the luminosity scale; quoting it separately would need the extractor to condition on `lumi` alone (a one-line variant of `ComputeStatCov`, not implemented). **The same treatment now runs on charge_asym / FBratio too** (2026-09-15, same decision) — see those bullets; there the box is a thin sliver because lumi cancels in a ratio. **The COUNT record is untouched** — `comb_W_yields.csv`/`comb_summary.csv` keep the fitted yields (the fit's record, echoed as a console cross-check; every observable takes only the r's from these files, the COUNT columns feed nothing but `xsec_fiducial_diag`). charge_asym/FBratio stayed count-based until 2026-09-22 (on the "A×ε cancels in a ratio" argument) and are built from r × σ_gen too since then — the argument holds for A_ch but not for R_FB, see the `fiducial_yields.C` bullet. First numbers, per flavour (post η-convention fix): met σ(W⁺)=50.36±0.91 / σ(W⁻)=39.68±0.80 / W=90.03±1.21 nb (eff. r≈1.01); leppt_mt40 54.20±1.37 / 43.89±1.08 / 98.09±1.75 (eff. r 1.09/1.11) — **the ~9% met↔leppt_mt40 spread = the e pT-tail mismodeling absorbed differently, now THE visible inter-discriminant systematic**. **2026-08-24 (post-filter data + in-fit-ABCD fit, leppt_mt40): σ(W⁺)=58.71±2.03 / σ(W⁻)=47.44±1.68 / W=106.15±3.44 nb, eff. r 1.180/1.198/1.188 — +8% vs the pre-filter era (the predicted filter-recovery shift) and eff. r right at the r-scan crossing 1.18–1.22; the met variant NOT rerun (its fit still pre-filter era — refit on lxplus before quoting).** → `plots/comb/xsec/<disc>/{W_fiducial,W_dsigma_deta_comb}`. **`xsec_fiducial_diag(disc)` — PER-FLAVOUR EXTRACTION DIAGNOSTIC (2026-08-12)**: two plots (μ, e), each with 4 series × {W⁺,W⁻,W} — σ_gen-fid, σ_reco=S_f/L, r×σ_gen (the measurement, identical in both since r is shared), r×σ_reco=N_fit,f/L (count-based) — so each gen/reco pair differs by exactly (A·ε)_MC,f, printed on the plot. Per-flavour S_f comes from THIS repo's `plots[/Elec]/combine_input_W<suffix>.root` (`W{p,m}_lab_y*/signal` integrals; μ+e reproduce the CSV `signal_prefit` exactly), r-errors via the same `sumWithCov`; Σ_flavour of r×σ_reco reproduces the fork's `*_sum_yield/L` exactly (printed cross-check) → `W_xsec_diag_{mu,ele}` + `xsec_diag.csv`. **(A·ε)_MC: μ 0.90/0.94/0.92 vs e 0.76/0.77/0.76 (W⁺/W⁻/W, met)** — the μ/e gap is the electron ID+ECAL-crack loss, i.e. the visual case for the ECAL-gap TODO; leppt_mt40 shifts them to μ 0.89 / e 0.74 (the m_T cut's efficiency). **The η-binned diagnostic (`W_dsigma_diag_<flav>_<charge>` + the `AxEps_vs_eta` summary) resolves the crack directly**: both electron curves dip to ≈0.62–0.64 in exactly the two |η_lab| ∈ [1.2,1.6] bins (containing 1.4442–1.566), with no such feature in muon — and, being detector response, the W⁺ and W⁻ curves now lie on top of each other (that charge-symmetry is what exposed the η-convention bug above; **use it as the standing sanity check on this plot**). Muon A·ε is ~0.96 at both η edges and dips to 0.81(W⁺)/0.92(W⁻) centrally — the edge rise is pT-threshold in-migration (the lepton pT spectrum is softer and falling at 25 GeV at large |η|, so smearing adds more than it removes), which is why it is charge-dependent while the crack is not. The bin-by-bin correction is not an unfolding, so per-bin σ carries this migration — and, since 2026-09-22, so do the A_ch and R_FB built from the same r × σ_gen (their count-based predecessors carried the raw detector A×ε instead, far larger in R_FB); the inclusive σ does not. Both `xsec_fiducial_comb` and `_diag` share `loadCombIngredients` + `sumWithCov` (single-source of the G/S/R/cov ingredients). Run by `run_observables.sh` (comb chain). The μ-vs-e comparison (meaningless under the grand fit's shared r's) comes from the per-flavour fits since 2026-09-22 — `xsec_fiducial_flav` above.
- [plotting/postfit_incl.C](plotting/postfit_incl.C) — **rapidity-INCLUSIVE simfit postfit stacks (2026-08-05)**: `postfit_incl(disc)` sums the 12 LAB bins per flavour for W⁺/W⁻/W (6 plots) → `plots/comb/postfit_incl/<disc>/postfit_{mu,ele}_{Wp,Wm,W}`. **`postfit_incl_fit(disc, fit)` (2026-09-22)** = the same for one per-flavour fit (`simfit_mu` / `simfit_ele`, fork mode `flavfit`), only its own flavour, from its own CSV, sidecar and `fitDiagnostics` (its work dir via `plotting/fit_variants.h`, which also makes the comb path honour `$FORK_TEST`) → `plots/flavfit/postfit_incl/<disc>/postfit_<flav>_{Wp,Wm,W}`; the header says "(#mu-only simfit postfit)". Built WITHOUT fitDiagnostics: since every simfit parameter is a pure normalization (no shape nuisances — the 2026-08-17 lnN nuisances are rate-only too), postfit shape ≡ prefit template × fitted scale, so the stacks are reconstructed exactly from `combine_input_W*.root` × the fitted params in `comb_W_yields.csv` (r per bin, r_Z, the qcd multipliers, and since 2026-08-17 the lumi multiplier on all MC — column absent ⇒ 1, so old CSVs still work; ndf auto-adapts by detecting the shared-QCD layout). **abcd-aware since 2026-08-23**: reads the CSV's 18th column `qcd_model` and, when `abcd`, multiplies the `qcd_abcd` template instead of `qcd` (the reconstruction stays exact — the CR scales are pure normalizations too; column absent = legacy behavior). **Shape nuisances (2026-09-07):** when the fitted cards' sidecar (`simfit/datacards/qcd_lnn_kappas.txt`, line `lheSysts`) lists LHE shape nuisances, the macro reads the per-channel POSTFIT shapes from `simfit/fits/simfit_lab/fitDiagnostics_simfit_lab.root` (`shapes_fit_s/<F>_<C>_lab_y<i>/<proc>`, unit-width bins remapped onto the input's axis as in the fork's `draw_postfit_pO.C`, summed over y; pulled by `sync_lxplus.sh download`) — and if that file is absent it falls back to prefit×scale with a loud `[WARN]` and an "APPROX: shape θ's ignored" label (the sidecar of a pre-09-07 fit has no `lheSysts` line ⇒ the exact shortcut, unchanged). Optional args 5/6 = `lheSysts` override ("none" forces the shortcut) and the fitDiagnostics path. Same cosmetics as the per-bin postfit plots (pull pad); info box: data/total, Baker-Cousins χ², r_Z. First results (leppt_mt40): μ W-incl χ²/ndf 1.46 (p 0.03); **e W-incl 3.43 (p≈0) — systematic data excess over the model across the pT tail 45–95 GeV**, i.e. the electron pT-shape mismodeling aggregated (links: e-scale/calibration shift, unembedded MC, anti-iso shape tilt; the per-bin fits partially absorb it via r/qcd_norm). Run by `run_observables.sh` (comb chain). **2026-08-24 (abcd fit, post-filter data): μ W-incl χ²/ndf 5.96, e 13.17** — NB ndf≈10 (38 bins − 28 params), per-bin rms residuals ~1.6σ (μ) / ~3.5σ (e); patterns = the diagnosed 25–32 turn-on deficit (BOTH flavors) + the e required-QCD tilt + tail excess; μ total overshoots by 95 evts while e undershoots by 67 (the shared-r compromise). Strengthens the TnP-SF + FF-shaped-e-template follow-ups. **2026-09-15, THE TURN-ON DEFICIT RE-MEASURED AFTER THE MUON SFs (user flagged that the μ plots now carry them — the 08-24 "flavor-universal ⇒ missing TnP SFs" reading is HALF-SUPERSEDED):** on the 09-14 night fit (μ has ID+ISO+trigger SFs since 09-14, ⟨SF⟩ ≈ 0.988; the electron still has NONE), data/postfit summed over the 12 lab y bins × 2 charges — μ [24,26) 0.873, [26,28) 0.949, [28,30) 0.907, [30,32) 0.923, then **1.061 already at [32,34)**; e 0.904 / 0.861 / 0.878 / 0.900 / 0.931, not back above 1 until [34,36). Over 24–36 GeV: **μ −3.4% (pulls −1.0…−1.9σ) vs e −7.6% (−1.4…−2.8σ)**; whole range μ −1.5%, e +1.0%. So **the muon deficit SURVIVED the SFs but is less than half the electron's and turns over one bin earlier** — it is no longer flavor-universal in magnitude, which weakens the common-cause reading. The remaining μ candidates are the SFs NOT applied: **lepton momentum scale/smearing** (nothing applied; migrates events straight across the 25 GeV cut, i.e. concentrated exactly here — the leading suspect), the **reco/tracking** efficiency (no SF supplied), and the pp-SF-in-pO caveat (cf. the Z→μμ iso observation: 0.20→0.15 cut 6.2% of data vs 2.9% of DY MC, i.e. the data efficiency sits BELOW the unembedded MC's while the pp iso SF pushes MC the other way). **NOT the trigger SF** — it is one flat inclusive number by design (0.9971), but `trig_eff_mb.C` measured the μ efficiency flat in pT from 10 GeV (L1 open muon, no threshold), so a flat SF loses nothing there. For the electron the 08-24 reading stands unchanged (no SFs of any kind yet, larger deficit).
- [plotting/syst_shapes.C](plotting/syst_shapes.C) + [plotting/run_syst_shapes.sh](plotting/run_syst_shapes.sh) — **diagnostics of the shape systematics in the Combine inputs (2026-09-07; since 2026-09-14 also the combined muon-SF nuisance muSF in the muon channel — the per-flavour runtime list `gSystNames` mirrors `mtandmet.C`, LHE names first)**, `./run_syst_shapes.sh [met|leppt_mt40|all]` → `plotting/logs/syst_shapes_<disc>.log` (the record). (1) `plots[/Elec]/syst_shapes/<disc>/perbin/<region>_<process>`: Up/nominal and Down/nominal of every systematic (nPDF/qcdScale/alphaS red/blue/green, muSF magenta; legend and info box grow with the count) on one canvas per MC template (416 per disc), info box = integral shifts, 68% range of the per-bin ratios, one-sided bins; (2) `summary_<syst>_{lab,fb}`: integral shifts of `signal` vs rapidity bin — **nPDF grows from ±2.5% centrally to +5.7/−8.2% in the most forward bin y11** (the oxygen x-range), qcdScale ≈ +4/−7% flat (μF×½ dominates), alphaS ±0.5%; (3) INCLUSIVE CONSISTENCY (`[INCL]` lines): (a) the true inclusive uncertainty = the member twins summed over the 24 lab templates then combined (`pOLhe::Hessian` with the exact LHAPDF CL factor √(χ²q(0.6827,1)/χ²q(0.9,1)) = 0.607957, i.e. "÷1.645" is a rounding) vs (c) the all-events reference of `lhe_weights.txt` — selected W⁺ nuclear +1.8/−3.2% vs +2.0/−3.8% all events, W⁻ 1.7/2.9 vs 1.7/3.1, DY 1.8/3.6 vs 1.7/3.3; scale +3.9/−6.7 vs +4.0/−6.6; alphaS ±0.46 both; the baseline products (no all-events reference) ≈ +1.9/−2.1% at 68% (= the 3.0–3.7% at 90%); and (b) the sum of the per-bin Up/Down templates over the 24 regions = what ONE collapsed nuisance implies in the fit: **(b)/(a) = 1.27 (up) / 1.16 (down) for the signal nPDF, 1.01 for qcdScale, 1.00 for alphaS** — the collapse over-estimates the inclusive nPDF uncertainty by ~20–30% (larger for the sparse z/ztau/wtau, 1.3–1.5); the per-eigen-direction treatment would remove that; the lepton-SF family prints a `(b)` line (for it (b) IS the inclusive shift — a coherent ±1σ of a normalization factor, nothing divided out: signal muSF +0.35/−0.38% = the three sources muID ±0.05%, muIso ±0.28%, muTrig +0.20/−0.24% in quadrature, 2026-09-14); (4) `[CLOSURE]`: the stored LHAPDF `_nPDFUp/Down` vs `pOLhe::Hessian` recomputed from the twins agree to 1e-8 (two implementations, same data); (5) **`members/<region>_signal` (2026-09-08, user request "12 plots per rapidity bin, nominal + all variations on the same plot")**: one canvas per W region (48 per flavor per disc) with the nominal `signal` (black) and ALL 106 EPPS21 member templates (grey, α = 0.25) on the same axes (log y, the fit's own 2 GeV bins; a ±3% band is 1–2 px at that scale, so a ZOOM INSET shows the peak bins ≥ 85% of the maximum on a linear y axis spanning nominal ± max(7%, 1.25 × the largest Up/Down deviation), window marked by a dashed blue box on the main plot), rebuilt from the skim `_epps21` twins with mtandmet's k_s and W⁺+W⁻ sample sum (member 0 ≡ the Combine-input `signal` bit for bit — the second `[CLOSURE]` line), ratio pad = member/nominal (grey) + the region's LHAPDF `nPDFUp/Down` in red; the red envelope lies OUTSIDE every single member (quadrature over 53 eigen-directions): central bins members within −1.6/+1.2% vs Up/Down +2.2/−2.6%, y11 members −6.3/+3.5% vs +5.7/−8.2%. One-sided bins (Up and Down on the same side; Combine warns, interpolates anyway): alphaS 61 (leppt_mt40) / 171 (met) in low-stat tails, nPDF/qcdScale a few in floored empty templates. **TAIL-BIN qcdScale SPIKES — DIAGNOSED 2026-09-15, benign (user asked about a +22/−25% spike at lepton pT 92–94 in `Wm_fb_y1/signal` while every neighbouring bin sat at ±2%):** NOT a Poisson count fluctuation but **a single POWHEG event with a freak scale weight**. Measured in the W⁻ mu ntuple (1.17M events, ⟨w₀⟩ = 5463, 0.69% negative): 98 events (0.008%) have |w₁| > 5⟨w₀⟩ for member 1 (μF×2), 19 have > 10⟨w₀⟩, the largest |w₁| = 1.26e6 = **231×⟨w₀⟩** (that event's w₁/w₀ = 228) — the nominal weights are ±constant (|w₀| = 5538.8 EXACTLY for every event, 0.69% negative) and the tail lives in the VARIED weights. **MECHANISM RESOLVED 2026-09-15 — and the constant |w₀| is what HIDES it, not evidence against a cancellation** (an intermediate reading here said "no near-zero-w₀ population ⇒ not an NLO cancellation"; that inference was wrong): POWHEG generates UNWEIGHTED events, so |w₀| is one number by construction and its sign is the sign of B̄ = B + V + ∫R. The varied weight is `ttbar_w[k] = w₀·B̄_k/B̄_0`, so the near-cancellation lives in B̄_0 — a denominator unweighting ERASED. Decisive evidence it is the shared denominator and not any one member: on the worst event **α_s ± 0.001 gives ρ = 11.5 / −10.2**, impossible for a 0.5% physical variation, and μF×2/μF×½/μR×2/μR×½ all blow up together (228 / −199 / −97 / 123) with four sign flips. Extreme events are 40× enriched in negative w₀ (28% vs 0.69%), i.e. they sit ON the B̄ = 0 boundary. Guarded at the source since the same day — see `pOLhe::kMaxMemberRatio` in the `lhe_index.h` bullet. In a well-populated bin such an event is diluted; in this bin (27 EFFECTIVE events, 0.100 events after k_s, 0.08% of the 131.9-event region) it is a large fraction: the whole variance excess is one event, √(Σw₁² − Σw₀²) = √2.22e9 ≈ 47000 ≈ 8.6⟨w₀⟩ with flipped sign (Σw₁ DOWN 121884 vs 156984 while Σw₁² is 3× up; member-1 N_eff collapses 27.2 → 4.7). **The rate closes quantitatively:** p = 98/1.17M per event ⇒ 0.23% per 27-event bin ⇒ **~4 affected bins predicted among the 1824 non-empty signal bins of the 48 W regions; 5 observed** with |Up−1| > 10% (qcdScale only; nPDF 0, alphaS 0, muSF 0), carrying **0.007% of the signal yield** — no effect on r or on the θ constraint (0.1 expected events = flat Poisson term). The max/min ENVELOPE is what exposes it (a max is one-sided: it seeks whichever member fluctuated highest); nPDF is milder (±8% worst bin) because PDF members only re-scale rather than reshape, and muSF is dead flat (a smooth per-event factor, no MC-stat content) = the control that proves the spike comes from the WEIGHTS, not the bin. **The 6-point choice already suppresses it**: the envelope here is Up = member 2 (1,½) 1.224 / Down = member 1 (1,2) 0.755, whereas `--scale-points 8` would give 1.281/0.616 set by the dropped member 7 (½,2), whose N_eff is 1.6. NB the twin STORES all 9 members; only [0,1,2,3,4,6,8] enter the envelope (`SCALE_MEMBERS` in `lhe_updown.py`, dropping the antagonistic corners 5 = (2,½) and 7 = (½,2)) — don't read a stored member as an envelope input. Also note the area normalization WIDENS this particular bin (raw 1.169/0.776 → 1.224/0.755, the member integrals moving the other way, 0.955/1.029). **FIXED AT THE SOURCE 2026-09-15 — `pOLhe::kMaxMemberRatio = 10` (see the `lhe_index.h` bullet), full MC re-skim + `run_lhe_updown.sh all` + all four Combine inputs regenerated the same day.** A per-bin N_eff floor in `lhe_updown.py` was considered and REJECTED: it repairs the symptom in the stored twins, would break this macro's `[CLOSURE]` line unless the logic were duplicated here, and did not even catch the worst case (that bin has 3687 effective entries — it is not a sparse bin). **The visible manifestation was not this tail bin but `summary_qcdScale_fb` / `_lab`, where W⁻ y5 (fb) and y6 (lab) stood at ±1.26% / ±1.09% while every other bin sat at ±0.3–0.5%** — 45% of that shift came from ONE event (`July_29_MC_Wm_mu` entry 70485, ρ = 228 scale / 261 nPDF, reco muon pT 36.5 in the 36–38 bin, η −0.23). After the guard: **y5_FB qcdScale +1.258/−1.262% → +0.37/−0.37%, nPDF +1.404/−1.421% → +0.52/−0.51%; y6 lab +1.086 → +0.36%, +1.237 → +0.51%** — level with the neighbours (y4_FB +0.46, y7 +0.37). At the Combine-input level the muon signal bins above 10% went nPDF 1 → **0**, qcdScale 6 → 4, the carried yield 0.77 → 0.46 of 6997 events; the electron is unchanged (its outliers are below the cap). The 4–6 survivors are pT 80–98 bins with ρ ≈ 4–5 in 25–43 effective entries — 0.007–0.010% of the signal, 0.04–0.25 expected events each, invisible in any integral plot and unable to move a flat Poisson term. A tighter cap would reach them at the price of removing a genuinely larger population; not done. NB `TString::operator()` returns a `TSubString` whose `Data()` points into the FULL string — copy to a `TString` first (bit this macro once).
- [plotting/xsec_contour.C](plotting/xsec_contour.C) — **the (σ_W, σ_Z) confidence ellipse from the simfit POIs (2026-09-15; the plane for comparing nPDF sets).** `xsec_contour_WZ(disc, binning="lab")`. **Per-flavour fits (2026-09-22):** the loading/propagation/scan/staleness-gate code is factored into `LoadContourFit` (a `ContourFit` struct per fit) and the plot + CSV into `DrawContourSingle`, so `xsec_contour_WZ` is unchanged in output (CSV byte-identical but the scan-path line, verified) and two entries were added: `xsec_contour_WZ_fit(disc, bn, fit)` = the same single-fit plot for `simfit_mu`/`simfit_ele` (their own `--contour` scans from `simfit_<flav>/contour/`) → `plots/flavfit/xsec/<disc>/xsec_contour_WZ_<bn>_{mu,ele}`, and **`xsec_contour_WZ_flav(disc, bn, withComb=true)` = the μ-only / e-only (/ grand) regions OVERLAID** → `plots/flavfit/xsec/<disc>/xsec_contour_WZ_<bn>.{png,pdf,csv}`: per fit in its colour the 68% region filled + solid and the 95% dashed, PROFILED when its scan exists and passes the gate, else Gaussian (the legend says which), over the common r = 1 point and EPPS21 cloud; a WIDE canvas (1150×800) with all text in a panel right of the frame, because the μ and e regions sit at different places and together cover the in-frame corners the single-fit layout relies on; the console/CSV give σ_e/σ_μ for σ_W and σ_Z with stat errors (lumi — the +ρ diagonal — moves both regions together). Both inclusive cross sections are LINEAR in the POIs — σ_W = Σᵢ r_i σ_gen-fid,i and σ_Z = r_Z · σ_gen-fid,Z — so their joint uncertainty is J V Jᵀ with V the **25×25** `h_cov_poi[_FB][_stat]` (POI space, order [r_Wp_y0..11, r_Wm_y0..11, **r_Z**], axis-labelled, read BY LABEL) that the fork's extractor gained the same day: `h_cov_yield` deliberately covers the 24 W POIs only, so the σ_W↔σ_Z cross term had no home. **No refit was needed** — `run_pO_fits.sh --extract-only` on the downloaded 09-14 tree produces it. **STAT-ONLY PROFILED CONTOUR (2026-09-15c, user request):** `--contour` now runs a SECOND grid scan per variant with every constrained nuisance frozen at its post-fit value (`--freezeParameters allConstrainedNuisances --setParameters <postfit>` — the identical recipe as the `--statonly` companion fit), with its OWN ±4σ window from a frozen `--algo singles` (the stat region is ~3× smaller per axis, so re-using the total window would spend most of the grid outside it) → `higgsCombine_contourstat_<B>.MultiDimFit.mH*.root`, gated for staleness exactly like the total one. **Colour convention on this panel: BLACK = total, DARK RED = stat only; DASHED = the Gaussian (covariance) region, SOLID = the profiled scan** — so the four objects read as a 2×2. The two scan stems cannot collide (`contour` vs `contourstat`, and neither matches the `…statfit`/`…fit` singles pre-fits — verified). Draws 68%/95% filled ellipses (Δχ²=2.296/5.991, eigen-decomposition), the stat-only 68% dashed, the r=1 theory point (POWHEG + the generation nPDF = EPPS21nlo_CT18Anlo_O16 central), and — when a `--contour` scan exists — the EXACT profiled contours on top (see below). Writes `xsec_contour_WZ_<binning>.csv` (the point, the 2×2, and the 3×3 over (W⁺,W⁻,Z), total + stat). **Numbers (leppt_mt40, lab, the 09-14 night fit): σ_W = 106.47 ± 3.51 nb (1.31 stat, 3.25 syst) — every digit identical to `xsec_comb.csv`, i.e. the 25×25 POI path and the 24×24 yield path agree, an independent check of both; σ_Z = 9.522 ± 0.472 nb (0.375 stat, 0.286 syst); ρ = +0.521.** ρ is **lumi**: 3% fully correlated on both gives ρ_lumi = 9/(3.30×4.96) = 0.55, slightly reduced by the DY-background anticorrelation — and the stat matrix shows exactly that decomposition, cov_stat(W,Z) = −0.020 nb² (NEGATIVE: DY is a background in the W channels) against cov_total = +0.862 (lumi), while cov_stat(W⁺,W⁻) = +0.003 ≈ 0 (disjoint event samples). **The lumi cancellation in the ratio is the validation**: σ_W/σ_Z = 11.18 ± 0.48, of which **0.083 is systematic (0.74%) vs 3.06%/3.00% on the individual σ** — a factor ~4. NB the ratio's TOTAL error is nonetheless *larger* than σ_W's, because σ_Z is stat-limited (~360 Z events, 3.9%); don't quote "the ratio is better measured" without that qualifier. `fb` runs too (a DIFFERENT lab window — Σ over the FB edge set: σ_W 91.82 ± 3.07 vs theory 76.60, r_eff 1.199 consistent with lab's 1.191). **PROVENANCE GUARDS on the profiled contour (2026-09-15b, user asked to make sure the contour code uses the CONTOUR output, not the nominal one):** (1) the scan file is GLOBBED as `higgsCombine_contour_<B>.MultiDimFit.mH*.root` in `simfit/contour/contour_<B>/` rather than hardcoding Combine's `mH120` — and the prefix cannot collide with the `--algo singles` pre-fit, which is named `higgsCombine_contour**fit**_<B>...` (verified). (2) **STALENESS GATE:** a reparametrization cannot move the minimum and a MultiDimFit grid's first row IS the exact best fit, so the scan's minimum must reproduce Σᵢ rᵢσ_gen,i and r_Zσ_gen,Z; if either disagrees by > 1e-3 relative the scan is from a DIFFERENT fit (or a different `gen_xsec.C` run — the σ_gen weights are baked into the contour workspace) and the macro **REFUSES to draw it**, printing both minima and what to do. The noise floor is the tree's `Float_t` storage, measured at 1.2e-07, so the tolerance sits ~4 orders of magnitude above it. (3) the plot CARRIES A STAMP — green "profiled contour: from the MultiDimFit scan", red "scan REJECTED (different fit)", grey "not available (Gaussian only)" — and the ellipse legend entries say "(Gaussian)", so a PNG can never be mistaken for having the profiled region when it does not; the CSV records `contour,status` = drawn|rejected_stale|absent plus the file and the scan minimum. (4) a wrong file that happens to have the branches (e.g. the singles pre-fit) still trips the `n < 10 usable points` check. All three states exercised 2026-09-15b with synthetic scans (matching / shifted by 1% / absent). Needs `h_gen_sig_Z` (see `gen_xsec.C`) and `h_cov_poi`; missing → a message naming the command to run. Run by `run_observables.sh` (comb chain) → `plots/comb/xsec/<disc>/xsec_contour_WZ_<binning>.{png,pdf,csv}`. **EPPS21 MEMBER SCATTER (2026-09-15c, user request "scattered point for all variations of EPPS21 member so we know the coverage"):** the 106 variations are drawn as small translucent green dots at (Σᵢσ_gen,i(m), σ_gen,Z(m)) from the `_epps21` twins above, UNDER the central diamond (= member 0); the axis range grows to contain them and the CSV records `n_epps21_members` + the min/max per axis. The cloud is elongated ALONG the data ellipse's major axis — nPDF variations move σ_W and σ_Z together, the same direction the lumi correlation points — which is exactly why the discriminating power of this plane sits in the SHORT axis. **Cosmetics (2026-09-15c, user request):** axis titles 0.036 / labels 0.030 (down from 0.045 / 0.040), legend 0.0245, numbers 0.026, stamp 0.023, one header line instead of two — the shrink is what buys room for the IN-FRAME CMS banner `CMS_lumi(c, 13, 10)`, as in every other plot of the repo. The numbers block hangs off the legend's computed bottom edge, so it follows when an entry appears or disappears (stat ellipse, profiled contour, member scatter). Before that (2026-09-15 only) this was the one plot using `CMS_lumi(c, 13, **0**)` (banner ABOVE the frame) instead of the repo's usual `10` — r_eff ≈ 1.19 pins the ellipses to the upper right and the whole left column is needed for header+legend+numbers; it also locally sets `relPosX = 0.13` around the call (CMS_lumi.C's own commented-out `if( iPosX == 0 ) relPosX = 0.12;`) so "Work in Progress" does not land on "CMS", restoring the global afterwards.
- [plotting/plotZcurve.C](plotting/plotZcurve.C), [plotting/plotRpOtheory.C](plotting/plotRpOtheory.C) — `plotRpOtheory.C` PRODUCES the theory graphs (`filelist_theory.txt` → `RpO_rootfile/RpO_FB_graphs.root`, channel/yield-independent); plotZcurve = Z curve plot
- [plotting/CMS_lumi.C](plotting/CMS_lumi.C), [plotting/plotting_helper.C](plotting/plotting_helper.C) — style / helpers. **Pull sub-pad (2026-07-19):** `SaveNicePlot1D_WithBkg` grew an opt-in bottom pull pad (`ps.pullPad`, default off): per-bin `(data − MC_total)/σ` with the Poisson-correct `σ² = MC_total + σ_MCstat²` (NOT the observed data error — that blows up empty-data bins against a small-error MC tail), drawn as zero-anchored bars (`"B"`), red dashed 0-line + dotted ±2σ guides, auto-symmetric y-range `max(3, 1.15·max|pull|)` (fix via `ps.pullYRange`). The canvas stays NEAR-SQUARE (height × `ps.pullCanvasScale`=1.125, so 800×900); the main pad compresses slightly — everything is pad-relative so the layout scales consistently. Enabled in `mtandmet.C` (all MT/MET stacks, both channels) and the fork's `draw_postfit_pO.C`; `dileptonpeak.C` left single-pad (flip `ps.pullPad` there to add it). **Opt-in x-zoom (2026-08-31):** `ps.xRangeLo/xRangeHi` (active when hi > lo) zooms the x axis in `SaveNicePlot1D_WithBkg` and `SaveDataMCRatio` — used to start every lepton-pT axis at the selection floor (24 = the 2-GeV bin edge enclosing the 25 GeV cut). `SaveNicePlot1D_WithBkg` RESTORES the full range on the input histogram after saving (it is often file-owned and written to Combine inputs later; persisted zoom bits would silently restrict no-arg `Integral()` calls downstream); `SaveDataMCRatio` zooms its clones only. Pull/ratio guide lines now span the VISIBLE range (bin edges of first/last shown bin), not the full axis. Re-synced to the fork same day. **Systematic boxes (2026-09-15):** `MakeSystBoxes(gTot, gStat, ps)` is the SINGLE SOURCE for "point + stat bar + syst box" — given a graph with the TOTAL error and its statistical twin it returns one `TBox` per point of half-height √(tot² − stat²), half-width `ps.systBoxWidthFrac ×` the point's own x error (falling back to `ps.systBoxHalfWidthAbs` for graphs at discrete x). `SaveNiceGraph` and `SaveNiceGraph_ErrorBand` take a trailing optional `gStat`: when given, **gStat becomes the drawn graph** (so the visible bars are statistical) and the boxes go under the markers — in `_ErrorBand` ON TOP of the translucent theory bands, which would otherwise tint them — with a legend entry added; when null both behave exactly as before, so every other caller (incl. the fork's `draw_postfit_pO.C`) is unaffected. Safe because every caller of this path sets its y-range explicitly in the graph tuner, so which graph establishes the frame does not change the axes. The stat > total guard warns only above 1e-3 relative: conditioning guarantees stat ≤ total exactly, but a point whose systematic is genuinely ~0 can come out a few 1e-6 negative through the numerical inversion (seen on `g_RFB_Wp` |y|bin 2), and a real graph mismatch is percent-level. Re-synced to the fork. **Overlays (2026-09-22):** `OverlaySeries` + `SaveNiceGraph_Overlay(series, …, g1..g4)` draw several measurements of one observable, each as point + stat bar + syst box in its own colour, shifted by `xShift` × the bin half-width with no horizontal bars (callers use ±0.35 for two series, −0.6/0/+0.6 for three, box half-width `systBoxWidthFrac` 0.20/0.17), frame = the unshifted bin range + TGraph's 10% margin, theory bands as in `_ErrorBand`; the tuner runs BEFORE the series styles (the shared tuners set a marker of their own); legend lower-left, or upper-left when theory bands take the lower-left. It REPLACED `SaveNiceGraph_ErrorBand_TwoData` (the legacy μ/e overlay, total bars only, removed with its only caller). Re-synced to the fork (compiles with `draw_postfit_pO.C`).
- [plotting/fit_variants.h](plotting/fit_variants.h) — **(2026-09-22) single source for WHICH fit a downstream macro reads**, the companion of `disc_variants.h`: `pOFit::Spec` per tag — `comb` (work dir `simfit/`, files `comb_*`), `simfit_mu`, `simfit_ele` (dirs and file prefixes = the tag) — with the flavour list, the plot label and the overlay style (μ blue circle 20 @1.3, e red square 21 @1.5, comb black diamond 33 @2.0: point-symmetric, centred sizes per the kMk* measurement, and no glyph+colour pair of the charge convention reused); `WorkDir/SummaryFile/Available` build the paths and honour **`$FORK_TEST`**, which `run_observables.sh` now exports ABSOLUTE so the macros read exactly the tree the driver checked.

### Stage 3 — `analysis/`

Final observable extraction from the rapidity-binned histograms:

- [analysis/analysis_helpers.h](analysis/analysis_helpers.h) — shared header. `pOAnalysis::YieldInRange`, `AsymErr`, `RatioErr`, and the `kPORapidityShift = 0.3466` constant (single source of truth for the pO→CM boost). **2026-08-04:** `AsymErr`/`RatioErr` take an optional trailing `cov` argument (default 0 = the old independent-yield formulas exactly) for the simfit covariance cross terms.
- [analysis/charge_asym.C](analysis/charge_asym.C) — `A = (N+ − N−) / (N+ + N−)` vs rapidity (12 bins, y ∈ [−2.4, 2.4]). **Production input since 2026-09-22 = a fit's FIDUCIAL yields `../skim/rootfile/fidyields_<fit>_<disc>.root` (`analysis/fiducial_yields.C`, same histogram names), so N± = r × σ_gen, never the fit's raw counts r × S** (the comb A_ch had been count-based until then; the switch moved it by ≤ 0.008 = ≤ 0.2 stat σ, because W⁺ and W⁻ of one bin share its A×ε). Any file with `h_yield_*` still works (e.g. raw skim histograms for a quick look). Reads `h_cov_yield` if present and includes cov(N+, N−) per bin; absent → legacy behavior. **STAT/SYST SPLIT (2026-09-15):** when the same file also carries `h_cov_yield_stat` (the fork extractor's conditioned covariance — see "STAT/SYST SPLIT" under "Downstream fit") the macro writes a SECOND graph **`g_chargeAsym_stat`**: identical points, statistical error only. **All three second moments must come from the stat matrix** — the `h_yield_*` bin errors are the TOTAL ones (= √ of the `h_cov_yield` diagonal, verified equal), so taking the diagonals from the histograms and only the cross term from the stat matrix would be wrong. `observables.C` then draws the bars from the stat graph and √(tot² − stat²) as a box. **Physics check the split must pass:** the lumi 3% is one nuisance fully correlated across every channel, so it rescales N⁺ and N⁻ coherently and CANCELS in A — measured syst 0.0026–0.0099 absolute vs stat 0.037–0.047 (2–10% of stat, largest in the edge bins y0/y1/y10/y11 where the per-(flavour,charge) QCD lnN bites; the same ranges on the fiducial graphs, 2026-09-22). A large syst here would mean the propagation is wrong, not that A is systematics-limited.
- [analysis/FBratio.C](analysis/FBratio.C) — forward/backward ratio in |y_CM|. Lab-frame `yEdges` chosen symmetric around Δy so the CM-frame bins are symmetric around 0 (required for F/B pairing). **Production input since 2026-09-22 = the FIDUCIAL yields `fidyields_<fit>_<disc>.root` (as for `charge_asym.C`) — for R_FB this is NOT cosmetic: F and B are different |η_lab| regions, so a count-based R_FB carries the detector A×ε ratio of two different places (the electron ECAL crack sits in F at |η_CM| ≈ 1.2 and in B at ≈ 1.9); the comb R_FB moved by up to 0.18 (−3.1 / +2.5 stat σ) in exactly those two bins, see the `fiducial_yields.C` bullet.** Reads `h_cov_yield_FB` if present: Var(F)/Var(B) get the within-sum cov terms (which `Yield::operator+` can't know) and the ratio error gets cov(F, B); absent → legacy behavior. **STAT/SYST SPLIT (2026-09-15):** `build_graph` gained a trailing covariance-matrix argument, so every graph is built TWICE — once from `h_cov_yield_FB` (the primary `g_RFB_*`, TOTAL) and, when present, once from `h_cov_yield_FB_stat` (**`g_RFB_{sum,Wp,Wm}_stat`**, statistical). The existing `setVar`/`covFB` machinery already takes every moment from the matrix, so the stat pass needs no other change; the `_stat` twins get NO legacy `g_RFB_mt_*` alias (nothing reads one). The console now prints the per-bin stat/syst breakdown. Same lumi cancellation as the asymmetry: syst 0.001–0.009 vs stat 0.058–0.105 (fiducial, 2026-09-22: syst 0.0015–0.0091 vs stat 0.057–0.108).
- [analysis/fiducial_yields.C](analysis/fiducial_yields.C) — **(2026-09-22) the ACCEPTANCE-CORRECTED yields of one fit**, Y_i = r_i × σ_gen-fid,i per (charge, y bin), lab AND fb, with the full r covariance (total + stat, from the fit's `h_cov_poi[_FB][_stat]`, read by axis label; an older extraction without it — e.g. the Aug-6 met tree — falls back to `h_cov_yield[_FB][_stat]`/(S_a S_b), then to the diagonal `rErr`/`rErr_stat`, and the console line `[fiducial_yields] <binning>: r covariance from …` says which), written under exactly the names `charge_asym.C` / `FBratio.C` read (`h_yield_W{p,m}_y*[_FB]`, `h_cov_yield[_FB][_stat]`) → those macros, unchanged, give the FIDUCIAL A_ch and R_FB. **WHY — a finding of the μ-vs-e work:** the extractor's yields r × S are RECO-level counts (S = L σ_gen (A×ε)_MC). In A_ch the W⁺/W⁻ A×ε of the SAME η bin nearly cancels, but R_FB divides two DIFFERENT |η_lab| regions (η_lab = η_CM + 0.3466), so each flavour's detector acceptance enters the ratio — the electron ECAL crack sits in the FORWARD bin at |η_CM| ≈ 1.2 and in the BACKWARD one at ≈ 1.9. Measured on the grand fit's own r's (identical physics for both flavours): the count-based R_FB of the μ and e yields differ by up to 60% (0.67 vs 1.11 at |η_CM| 1.2, 1.38 vs 0.84 at 1.9), A_ch by ≤ 0.015. With r × σ_gen (one pooled fiducial volume) both vanish. **POLICY (user, 2026-09-22): EVERY observable — σ, A_ch, R_FB, for the grand fit and the per-flavour fits alike — is built from r × σ_gen, NEVER from the raw fitted counts** ("we are not going to do a dedicated eff or acc correction" — r × σ_gen IS that correction, taken from MC); the count-based r × S survive only as the fit's record (CSVs, `<tag>_fitted_yields.root`), in the postfit stacks and in the `xsec_fiducial_diag` A×ε diagnostic. **So the PRIMARY comb A_ch/R_FB switched the same day** — they had been count-based since 2026-08-12 on the "A×ε cancels in a ratio" argument, while every σ has been r × σ_gen since that date. On the 09-21 grand fit R_FB_sum moved 1.026 → 0.848 at |η_CM| 1.88 (−3.1 stat σ) and 0.908 → 1.087 at 1.20 (+2.5σ), the other bins ≤ 0.04; W⁺/W⁻ alike (±0.16–0.20 in those bins, 1.7–2.3σ); A_ch ≤ 0.008 (≤ 0.2σ); every σ unchanged (`xsec_comb.csv` identical). Run by `run_observables.sh` for every fit (its `fid_observables` function); a partial output (< 48 yields) is deleted, never left for the next macro.
- [analysis/run_observables.sh](analysis/run_observables.sh) — **Module-5 driver (2026-08-03; simfit-aware 2026-08-04)**: `./run_observables.sh [met|leppt_mt40|all]` runs the whole fitted-yields→observables chain for ONE W-discriminant variant, carrying the disc tag through every filename/folder (`leppt` was retired 2026-09-21 and now exits 1 with a message naming the replacement; `all` loops the two survivors). Two conditional blocks per variant: **PRIMARY simfit chain** (when `pO_fit_out<suffix>/simfit/summary/comb_fitted_yields.root` exists): `fid_observables comb` = `fiducial_yields.C` on the grand fit → `../skim/rootfile/fidyields_comb_<disc>.root`, then `charge_asym.C` + `FBratio.C` on THAT (covariance-aware) → `{charge_asym,FBratio}_fid_comb_<disc>.root` (fiducial since 2026-09-22 — before, the two macros ran on the count-based `comb_fitted_yields.root` → `*_fit_comb_<disc>.root`, now orphaned) → `observables_comb(disc)` → `plots/comb/{charge_asym,FBratio}/<disc>/` → `xsec_fiducial_comb(disc)` + `xsec_fiducial_diag(disc)` + `xsec_contour_WZ(disc, bn)` **for BOTH `bn` ∈ {lab, fb}** → `plots/comb/xsec/<disc>/`, and `postfit_incl(disc)` → `plots/comb/postfit_incl/<disc>/`; **μ-vs-e chain (2026-09-22; replaces the legacy per-flavour chain, removed the same day with those fits)** (when `pO_fit_out<suffix>/simfit_{mu,ele}/summary/simfit_<flav>_fitted_yields.root` exist, either flavour alone works): `fid_observables` on each flavour fit (the grand fit's files come from the comb block and serve as the overlays' reference) → `../skim/rootfile/fidyields_<fit>_<disc>.root`, `charge_asym.C` + `FBratio.C` on those → `{charge_asym,FBratio}_fid_<fit>_<disc>.root`, then `observables_flav(disc)`, `xsec_fiducial_flav(disc)`, for BOTH `bn` the overlay `xsec_contour_WZ_flav(disc, bn)` + each fit's own `xsec_contour_WZ_fit(disc, bn, simfit_<flav>)`, and `postfit_incl_fit(disc, simfit_<flav>)` → `plots/flavfit/{charge_asym,FBratio,xsec,postfit_incl}/<disc>/`, every output `require_file`-checked. `FORK_TEST` is now resolved ABSOLUTE and EXPORTED, so the macros (`plotting/fit_variants.h`) read the tree the driver checked. Verified end to end 2026-09-22 on a scratch fork tree (the real grand fit symlinked + per-flavour extractions of it as fixtures, e ×0.95): every comb output bit-identical to before, the injected e/μ = 0.950 recovered, the contour staleness gate rejecting the perturbed fixture's scan. Real run the same day after the fiducial switch (leppt_mt40, the 09-21 grand fit): `[fiducial_yields] lab: r covariance from h_cov_poi`, `xsec_comb.csv` + `xsec_diag*.csv` identical to before, the contour CSVs identical but for the scan-path line (now absolute, `$FORK_TEST` being exported). Strict single-variant errors only when NEITHER block has inputs; `all` skips such variants. Post-checks outputs (root exits 0 even when a macro bails). bash-3.2-safe; `FORK_TEST` env overrides the fork location. **Latent bug fixed 2026-09-21:** `dsuf="$(disc_suffix "$disc")"` runs `disc_suffix`'s `exit 1` inside a command-substitution SUBSHELL, so an unknown tag did not stop the script — `dsuf` came back empty and the run silently proceeded with met's out-tree while writing into the other tag's folders. Now `|| exit 1`. **Second staleness bug fixed the same day:** the contour was called as `xsec_contour_WZ("$disc")` — no binning argument, so it only ever ran the default **lab**, and `xsec_contour_WZ_fb.{png,pdf,csv}` silently kept whatever the last MANUAL run had left. Caught on the 2026-09-21 fit cycle, where the fb files still carried the previous fit's σ_W = 92.01 nb against the correct 93.22. fb sums a DIFFERENT lab window (Σ over the FB edge set), so it is a genuinely different measurement, not a cosmetic twin of lab. Now a `for bn in lab fb` loop with a `require_file` on each PNG, so a missing or unwritten variant fails loudly instead of leaving a stale file that looks current.

### `correction/` — corrections & studies

Orthogonal to the main skim→plotting→analysis flow: code that derives/validates
corrections and background estimates. **All macros are run from `correction/`**
and write outputs there (`correction/plots/`, `correction/rootfile/`). The two
that need `plotting_helper.C` `#include "../plotting/plotting_helper.C"`; the
isolation study is self-contained (writes/reads `correction/rootfile/`); the QCD
fit reads the *main* W skim output at `../skim/rootfile/`.

- [correction/qcd_abcd.C](correction/qcd_abcd.C) — **current** data-driven QCD/low-MET background via the ABCD method, **both muon and electron** (`qcd_abcd.C+` = muon, `qcd_abcd.C+(true)` = electron). Plane = relIso × PF MET (also stores relIso × m_T). Regions: iso-pass (`relIso < isoCut`) / iso-fail (anti-iso window `[isoFailLo, isoFailHi)`) × MET high/low (`metCut`); all four counts are QCD-only (`data − Σ EWK·k_s`, absolute MC from `mc_norm.h`). Reports the transfer factor `T = B/D`, the signal-region QCD `A = B·C/D`, and writes the iso-pass QCD MET/m_T **template** (anti-iso shape × T) to `correction/rootfile/qcd_abcd_{mu,ele}.root` for Combine, plus closure overlays (data vs EWK ± ABCD QCD). Channel diffs: electron uses `isoCut=0.095`, anti-iso `[0.20,1.0)` (muon `0.15`, `[0.30,1.0)`) and the `Wp_ele/Wm_ele/DYee` MCScale labels. Boundaries are projections → retune `ABCDConfig` with no re-skim.
**Template error treatment (FIXED 2026-08-19, `ScaleToIsoPass` — shared by the
MET/m_T and the pT paths, which had duplicated the loop):** per-bin errors are
the sideband SHAPE stat only (`T·σ_v`); the transfer-factor error δT is a single
fully correlated NORMALIZATION uncertainty and is folded into the total alone,
`σ_tot² = (T·σ_V)² + (V·δT)²` with `V = Σv_i`. Before, δT was written into every
bin and `IntegralAndError` then summed it in quadrature — i.e. `Σv_i²` where the
total needs `(Σv_i)²` — diluting it by the effective √N_bins (μ⁺ m_T>40 template:
quoted 4.3% vs correct 7.8%, factor 1.8). Central values UNCHANGED (verified:
0 bin-content diffs across all 306 templates in every `combine_input_W*.root`,
errors only) so **no re-fit was needed**; the fit never read these errors anyway
(no `autoMCStats` in the cards) — but the `mtandmet` pull pads did, and QCD-heavy
bins previously got artificially small pulls. Keeping δT out of the bins also
avoids double-counting the `qcd_rate_<flav>_<chg>` lnN if `autoMCStats` is ever
enabled. Written up in [docs/AN_qcd_background.tex](docs/AN_qcd_background.tex) (the AN
"QCD background" subsection: method, region/purity tables, both discriminants,
validation, systematics).
**Full diagnostics in the log (2026-08-19)** — every number the AN quotes is now
printed by the macro, so nothing lives in throwaway scripts. Per plane+charge
(buffered into `ABCDResult::diag` so it files under the right header, since
`runPlaneCharge` runs before `printResult`): the region composition table
(data | EWK | QCD | EWK/data), the closure ratio data/(EWK+A_pred) in A, `T` in
thirds+halves of the low-y band (**factorisation test**), and the anti-iso window
scan (both edges; note A_pred moves ~8–15% while the TEMPLATE TOTAL moves <1.5%,
because the total is anchored by the measured N_B). Then one `CHANNEL REPORT` per
flavour: T on both planes + transport, the **r-scan** (T_MET and T_mT vs the
assumed W scale, with the crossing point auto-located — μ⁺ 1.18 / μ⁻ 1.22 /
e⁻ 1.20, e⁺ none), the multijet fraction of all three selections, and the
assembled systematic budget → **two κ's: with and without the transport row**
(μ 1.24 / **1.14**, e 1.20 / 1.17). The without-transport μ value reproduces the
shipped κ_μ=1.15, which is the quantitative argument for keeping it.
ROOT prints all of this to stdout only, so use the wrapper
[correction/run_qcd_abcd.sh](correction/run_qcd_abcd.sh) (`./run_qcd_abcd.sh
[mu|ele|both]`, default both) — it pre-builds once (so the two jobs can't race
on the ACLiC artifacts), tees each channel to `correction/logs/qcd_abcd_<chan>.log`
(~290 lines, gitignored like every output dir; ~2 s to regenerate) and echoes the
T's and both κ's to the terminal.
**IN-FIT ABCD (2026-08-23 — the leppt_mt40 QCD normalization moved into the
Combine likelihood, `QCD_MODE=abcd`; the Mattermost-agreed "option 2+3 mix":
(m_T × iso) regions per Andre's proposal so the fitted variable pT is never an
ABCD axis, executed with Combine's rateParam-formula ABCD so the subtraction
floats with the POIs per Sitian's circularity objection. The SHAPE is
unchanged — still the anti-iso m_T>40 pT distribution.)** The macro now also
computes the m_T-plane regions with **C re-cut at m_T>40** (`kSRMtCut` = the
SR's own cut; B/D stay at m_T<yCut=30, the 30–40 band a buffer):
**A40 = B·C40/D = μ⁺ 131.9±10.5 / μ⁻ 126.6±10.2 / e⁺ 652.1±29.9 /
e⁻ 670.7±30.8** (vs the T_MET-normalized `qcd_pt_mt40` totals
164.0/155.8/729.6/690.8 — the option-(b) shift, now the abcd-mode prefit).
Exports the machine-readable **`abcd_counts_<lep><chg>`** TH1D (15 labeled
bins: `dataB,qcdB0,wB,zB,ztauB,dataC40,qcdC0,wC40,zC40,dataD,qcdD0,wD,zD,
T_met,T_mt`; B-region DY split DY vs DYτ since CR-B carries them as separate
processes; self-check dataB−wB−zB−ztauB ≡ qcdB0) into `qcd_abcd_{mu,ele}.root`,
consumed by `mtandmet.C`. New in the log: the per-charge in-fit block, the
window scan of the in-fit observable A40 (μ 9.5/2.9%, e 4.3/7.0% — NOT the
MET-plane Apred numbers), and the **REDUCED κ budget** (the stat rows are
profiled in-fit, the plane-transport row does not apply — the m_T plane IS the
counting plane): window ⊕ tilt → μ 1.09 / e 1.11, and with the FF total shift
replacing the tilt → **μ 1.09 / e 1.15 = the card defaults** (the FF shift
measures the same iso-pT correlation directly); the frozen-W residual
(wB/dataB×0.2 ≈ μ 2.4%/2.0%, e 0.7%/0.6%) is printed for `QCD_WCR=frozen`.
**NEW DIAGNOSTIC `runFFCheck` — the per-pT-bin fake factor** (from the
existing `h_pt_mt[_antiiso]_*` scan planes, NO re-skim): F(pT) =
QCD(iso-pass)/QCD(anti-iso) at m_T<30 in 5 GeV bins, pT≥25 (bins with
sideband < 3 events fall back to flat-T, ~2–3% of the yield), applied to the
anti-iso m_T>40 spectrum and compared with flat T_mT → `ff_pt_*` +
`ff_closure_pt_*` plots, `ff_pt_*` histos in the rootfile, per-bin table in
the log. **Result: μ F(pT) is flat (FF/flat−1 = +0.8%/−5.5%); the ELECTRON
F(pT) rises genuinely (−11% at [25,30) → +34…+64% at 40–55), total FF/flat−1 =
+12.4%/+13.2%** — the direct measurement of the e iso–pT correlation, larger
than the ⟨pT⟩-tilt proxy (7–9%), hence κ_e 1.15. If a pT-dependent transfer
ever needs to enter the fit per bin, the documented upgrade is Combine's
RooParametricHist (docs/abcd_rooparametrichist_tutorial in the fork).
**F(pT) UNCERTAINTIES REVISED 2026-09-24** (a reviewer question on the error
bars that mentioned "binomial"; investigated with a 5000-toy study of the
actual per-bin counts, all four charges). The old bars = `TH1::Divide` of the
EWK-subtracted histograms. **That IS the correct first-order error, and it is
the binomial one:** iso-pass (relIso < 0.15 / 0.095) and anti-iso ([0.30,1.0)
/ [0.20,1.0)) are disjoint windows filled exclusively in the skim, so N and D
are independent, and conditioning on n_p + n_f (n_p binomial, F = ε/(1−ε))
gives the identical σ_F² = F²(1/n_p + 1/n_f). What would be WRONG is option
"B" (σ² = F(1−F)/D, numerator ⊂ denominator): 20–60% too small here and
meaningless for F > 1. **What WAS wrong:** (1) the SHAPE in the thin high-pT
bins — above 40 GeV 20–60% (μ) / 13–35% (e) of the iso-pass data is
subtracted EWK, F is skewed and √n_obs shrinks the bar when n fluctuates low:
total coverage of ±σ is right (0.67–0.73 in every bin) but at pT ≥ 50 the
truth sits ABOVE the bar 19–32% and BELOW it 2–12% of the time (16/16
nominal), and μ⁻ [55,60) (n_p = 1 < b_p = 2.4) had F = −0.11 ± 0.09 (0.59
coverage) and was silently clipped off the plot by its y = 0 floor;
(2) CORRELATIONS — the flat line T is the pooled ratio of the SAME events
([25,30) alone is ~60% of it; the e [25,30) deviation is 3.5–3.9σ with the
correlation vs 2.4–2.8σ read naively); the `ff_closure` ratio pad ran
`TH1::Divide` on two predictions sharing the anti-iso m_T>40 counts a_i
(FF/flat ≡ F_i/T, a_i cancels), so a_i's Poisson error entered twice (bars up
to 45% too large, plus bars on the fallback bins whose ratio is identically
1); the printed totals (e⁺: FF 710 ± 43 vs flat 627 ± 15) are correlated and
must not be combined in quadrature; (3) the prefit r = 1 EWK subtraction is a
COHERENT systematic the bars never contained. **Now:** the plot shows the 68%
PROFILE-LIKELIHOOD interval per bin (n_p ~ Pois(Fμ + b_p), n_f ~ Pois(μ + b_f),
EWK MC as a known background — its MC stat is ≤ 0.6% of the variance — μ
profiled in closed form; MLE = the subtraction estimate whenever ≥ 0, so the
points do not move; for b → 0 it is the Poisson-ratio "binomial" interval),
T ± δT as a band, the fitted slope, and the tests in a top margin; the
`ff_pt_*` TH1D keeps the Gaussian error, the graph is written as `ff_pt_pl_*`,
the log table carries both plus the raw data count and the subtracted EWK.
**Flatness = likelihood-ratio tests over the 9 drawn bins** (one common F vs a
free F per bin, and vs F0·exp(s(pT−40)/10)), calibration printed from 1000
H0 toys (asymptotic 5% cuts fire in 3.3–5.9% / 4.2–5.5%): **μ⁺/μ⁻ flat, p =
0.82/0.84, slope ×1.10/×0.90 per 10 GeV (0.8/0.9σ); e⁺/e⁻ flat REJECTED, p =
5.2×10⁻⁴/4.3×10⁻⁴, slope ×1.34/×1.30 per 10 GeV at 4.8σ/4.4σ.** The χ² of the
Gaussian bars vs the flat line — what a reader forms from the old plot — is
NOT a valid test at these counts (its 5% cut fires in 15–26% of H0 toys; for
e± toy p = 0.14/0.09 where χ²(8) says 0.02): the old plot could not
demonstrate the rise the likelihood sees at 4.4–4.8σ. **The FF shift now has
an uncertainty** (1000 bootstrap toys; one-sided p from the H0 toys): **μ⁺
+0.9 ± 5.9% (p 0.46), μ⁻ −6.0 ± 4.7% (p 0.08) — consistent with zero, so the
muon row of the κ budget is statistical noise; e⁺ +13.3 ± 4.9%, e⁻ +13.0 ±
4.9% (p < 0.001) — a genuine ~2.7σ effect each**; printed under the κ-budget
row (κ recipe unchanged). Coherent r_W = 1.2 line: T −3.4%/−2.7% (μ),
−0.8%/−0.6% (e); drawn F bins all move DOWN, −1…−20% (μ, most at high pT) /
−0.3…−6% (e); e shift +13.3 → +13.0%, +13.0 → +12.7%. Verified: all 12
histograms of each `qcd_abcd_{mu,ele}.root` IDENTICAL (`compare_hists.C`,
incl. `ff_pt_*`, the templates, `abcd_counts_*` ⇒ no downstream re-run), the
logs differ only in the FF blocks + that one budget line, and reproduce byte
for byte (seeded toys). **The AN is now behind:** `docs/AN_qcd_background_infit.tex`
check (v), Fig. `fig:qcd:ff` and Tab. `tab:qcd:syst` still quote the
pre-embedding shifts +0.8/−5.5/+12.4/+13.2% without uncertainties and judge
flatness by eye — update the numbers + caption and re-copy `ff_pt_*.pdf`.
Consistency prints: Σ(anti-iso m_T>40 pT spectrum) ≡ C40 exactly (same events,
two binnings — the h_iso_pt_mt40 ↔ h_iso_mt identity); no pT×MET skim plane
exists, so the MET-sideband FF is only available pT-integrated (= T_MET).
**SR PREDICTION-CHECK STACKS (2026-08-31, `runInfitCheck`, prefit r=1
throughout — the group's "check the ABCD relation vertically" request):** the
counting relation A = B·C40/D is algebraically SYMMETRIC — horizontal
(A = C40×T, T=B/D, migrate along relIso) and vertical (A = B×R, R=C40/D,
migrate along m_T) give the same number AND the same per-bin template — so the
two views differ only in which variable stays differential, which is where a
relIso×m_T correlation would show: **`srcheck_mt_horiz_*`** = iso-pass m_T
spectrum, data vs EWK(r=1)+QCD(anti-iso m_T shape × T), 30/40 boundaries +
shaded buffer; **`srcheck_iso_vert_*`** = relIso spectrum at m_T>40, data vs
EWK+QCD(= the m_T<30 relIso shape × R), coarse variable relIso bins all on the
0.005 grid (cut and window edges are bin edges), the [isoCut,isoFailLo) gap
predicted too (bonus closure region); **`srcheck_pt_*`** = the same SR in the
FITTED variable, QCD = the in-fit `qcd_abcd` normalization (template
renormalized to A_pred exactly — the anti-iso pT projection loses the pT>100
overflow that C40 keeps, ≤1.1%, the same renorm the fit template gets in
mtandmet). All three carry A_pred = B·C40/D ± and A_actual = data−EWK(r=1) ±
+ ratio in the info box, per flavor×charge → `leppt_fit/` + a log block per
charge (incl. the consistency prints: iso-view iso-pass part ≡ A_pred, pT-view
A_actual ≡ m_T-plane A_actual). **NB A_actual rides the prefit W scale** — at
r̂≠1 it absorbs (r̂−1)·S_SR, and the printed S_SR sizes that: μ A_act/A_pred ≈
3.5 is entirely the r̂≈1.18 excess on S_SR≈2000 per charge (0.18·S ≈ 330–390 of
the 330 excess), NOT QCD non-closure; e⁺ 1.47 / e⁻ 1.13 = the familiar
required-QCD charge tension. First physics read: μ R is flat across the
anti-iso window (vertical closure per-bin good); e shows the mild R(relIso)
tilt (+1–2σ sideband under-coverage at high relIso), the iso–m_T face of the
known e correlation. **Three more families (2026-08-31, user follow-ups):**
**`srcheck_mt_vert_*`** = the ×R twin of the horizontal m_T stack — the SR
band alone (x 40–200), QCD = anti-iso m_T>40 shape renormalized to B×R,
narrated vertically; within the SR the two factorizations are BIN-IDENTICAL
(B×[antiiso_j/D] ≡ T×antiiso_j), so this shows the same prediction zoomed to
the entire SR (the pair demonstrates the equivalence — do NOT present it as
an independent prediction). **`srcheck_iso_horiz_*`** = the mirror ×T twin of
the vertical relIso stack — the iso-pass region alone (x 0–isoCut, the SR's
own relIso bins), QCD = the low-m_T iso-pass relIso shape renormalized to
T×C₄₀ (the renorm scale ≡ R numerically, the identity again); horizontal
narration (C₄₀ measured in the sideband, T migrates it across). Same caveat. **The place R vs T differ TESTABLY is
`stab_T_mt_*` / `stab_R_relIso_*`** — each transfer factor measured
differentially along the axis it is assumed constant over, vs the flat value
± stat band (dashed + red band), χ²/ndf over the solid bins printed on the
plot and in the log; open markers = the approach regions (T: the [30,40)
buffer; R: the relIso gap), the SR bins excluded (signal-blinded at r=1).
**χ²/ndf vs flat: T(m_T) μ 0.63/1.28, e 2.38/1.91 (+/−); R(relIso) μ
1.48/1.62, e 1.42/1.84** — the e iso transfer RISES with m_T (the transport/
FF finding in yet another projection: T(m_T<5) 0.25 → 0.38 at 10–20 for e⁺,
buffer +21%), μ T is flat; both flavors show a mild +13–17% R uptick in the
last relIso bin [0.8,1.0).
**pT display floor (2026-08-31): all lepton-pT plots start at the selection
floor** — `kPtAxisLo` = 24 (the 2-GeV bin edge enclosing the 25 GeV cut; the
5-GeV-rebinned ff plots use the true edge 25) — no more empty [0,25) band.
Applied to: the qcd_abcd pT closures/tilt/ff_closure plots AND the pT-plane 2D
maps (y-axis), `mtandmet.C` leppt_mt40 stacks (per-y + inclusive),
`postfit_incl.C` pT discs, `dataMC_kinematics.C` W-channel lepPt (Z channels
keep the full axis — legs go to 10–15 GeV), and the fork's `draw_postfit_pO.C`
(gated on the pT x-title; takes effect on the next lxplus fit/--draw-only run —
no local fits/ tree). `ptmt_scan.C` (pT-cut-related) deliberately untouched.
Mechanism: opt-in `ps.xRangeLo/xRangeHi` in `plotting_helper.C` (see the
Stage-2 bullet).
**Plot layout (2026-08-24): two SELF-CONTAINED views per lepton, one per fit
choice** — `plots/qcd_abcd_<lep>/met_fit/` (MET plane 2D + closures, the
no-m_T-cut pT display templates, and COPIES of the m_T-plane plots as the
transport cross-check) and `plots/qcd_abcd_<lep>/leppt_fit/` (the in-fit
counting plane incl. **`qcd2D_mt_infit_*`** — the region map with BOTH m_T
boundaries drawn, the 30–40 buffer and the relIso gap SHADED, regions labeled
B/D/SR/C40 — plus the m_T closures, the `_mt40` pT template plots, the tilt
checks and the `ff_*` diagnostics). rootfile contents/names untouched by the
layout. **AN files (in `docs/` since 2026-09-27): `AN_qcd_background.tex`(+pdf) = the full trail (method +
lnN history + validation); `AN_qcd_background_infit.tex`(+pdf, 2026-08-24) =
the clean version presenting the in-fit ABCD as THE method** (region map with
the buffer, likelihood equations, 4-charge FF figure, results + postfit
diagnosis; both PDFs built via a standalone wrapper — external \ref's to
other AN sections render as ??).
**Lepton-pT template (2026-07-29):** also writes `qcd_pt_{mu,ele}{Plus,Minus}` =
(anti-iso pT shape from `h_iso_pt_*`) × (the **MET-plane** T) — relIso×pT is
correlated for QCD, so the pT plane is never counted with the 2×2; template
total = the MET-plane iso-pass QCD by construction. **`qcd_pt_mt40_*`
(2026-07-30)** is the same built from `h_iso_pt_mt40_*`, i.e. the QCD template
of the pT-discriminant selection (pT>25 && m_T>40), same MET-plane T (QCD
isolation efficiency assumed independent of the recoil variable) ⇒ its total
IS the ABCD prediction for QCD surviving the m_T cut. **m_T>40 keeps only
≈26%/25% of the μ⁺/μ⁻ QCD and ≈33%/34% of the e⁺/e⁻ QCD** (μ 413+454→111+116,
e 1881+1687→623+579). Ships an anti-iso
slice-stability shape check (`antiiso_shape_pt_*`): μ stable; e shows a mild
tilt (sideband pT slightly soft vs iso-pass) — the known systematic knob. **Result:** muon QCD is small (T≈0.17; ≈3% of the high-MET signal region). Electron QCD was originally huge with the loose `eleCutIdWP95` (T≈0.30; ≈24–26%), which **motivated the ID switch to `eleMVAIdWP95`** (see the isolation/ID study + FIXMEs); after the switch the electron QCD dropped to **T≈0.32–0.38, ≈10–11%** of the high-MET signal region (iso-pass QCD ≈10037→3634, W signal only −6.5%, S/B 0.34→0.87). Closure is good in both channels. (Electron T/QCD came down a further ~4% on 2026-07-30 with the relIso bin-edge fix — see the `kNIsoAB = 200` note above. **MET-plane** values, μ⁺/μ⁻ and e⁺/e⁻: μ T = 0.2001/0.1832 **unchanged** by construction, e T = 0.3967/0.3368 → **0.3802/0.3243**; iso-pass QCD μ 413.7+453.8, e 1969.5+1758.9 → **1887.2+1693.9**. NB the m_T-plane has its own T (μ 0.178/0.163, e 0.343/0.323) — don't quote those as the MET-plane numbers.)
- [correction/ptmt_scan.C](correction/ptmt_scan.C) — **(pT × m_T) cut-pair scan** (2026-07-30, iso NOT scanned — stays at nominal): which (lepton-pT, m_T) cut combination optimizes the W selection. Signal = absolute W MC; background = the **anti-iso sideband** (full selection except iso, relIso in the qcd_abcd windows [0.30,1.0) μ / [0.20,1.0) e, ACTUAL MET/m_T) EWK-subtracted and shape-normalized to the measured iso-pass QCD total, + non-W EWK counted directly. (An earlier MET<5 "low-MET proxy" convention was **removed same day**: conditioning on MET caps the proxy's m_T at ~2√(pT·MET)≈30 GeV — a scan of m_T must not use a MET-conditioned background. Residual anti-iso caveat: mild iso–pT shape correlation, bounded by the `antiiso_shape_pt_*` slice checks.) Reads the joint `h_pt_mt_{mu,ele}{Plus,Minus}` (iso-pass) / `h_pt_mt_antiiso_*` (sideband) 2Ds from the skim. Outputs S/√(S+B) & S/B maps (optimum ★, tentative (25,40) ✚), sweep ROCs + AUC, m_T profiles → `plots/ptmt_scan_{mu,ele}/` + `rootfile/ptmt_scan_{mu,ele}.root`. **Result (2026-08-20 re-run on the post-2026-08-18 skim — SUPERSEDES every earlier number in this bullet; the ~10% data-only recovery from dropping `pclusterCompatibilityFilter` landed almost entirely in the DATA-DRIVEN B, so every FOM fell ~3%):** normalization inputs (integrals over the full scan plane = pT>20, no m_T cut) — μ data 8293 − EWK 5119.4 → **QCD tot 3173.6**, sideband 14239 raw / 14192.2 EWK-subtracted; e 12946 − 4141.4 → **8804.6**, sideband 29710 / 29646.3. **μ**: baseline (25,0) S 3978.4 / B 1692.7 (QCD 1282.0 + EWK 410.7) / S/B 2.35 / **S/√(S+B) 52.83** (was 54.40); (25,40) 57.56 (B 622.6, ε_S 96.7%); (25,45) **57.67** = best at pT≥25; global optimum **(20,42.5) 59.71** (ε_S 105.9%); scan floor (20,0) 50.04. **e**: baseline (25,0) 38.50 (was 39.85; B 3812.5 of which QCD 3436.2, S/B 0.85); (25,40) 46.89 (ε_S 96.9%); (25,45) 48.25; (25,52.5) **49.01** = best at pT≥25; global optimum **(21,50) 49.58**; scan floor (20,0) 32.15. So **pT>25 costs 3.4% (μ) / 1.1% (e)** vs the unconstrained optimum, and **(25,45)** sits 0.0% (μ) / 1.5% (e) below each channel's own pT≥25 optimum — conclusions UNCHANGED: **pT>25 stands** (raising it only cuts the Jacobian peak), μ m_T plateau 40–45 (57.56/57.66/57.67 — flat to 0.2%), e rises slowly to ~52.5; **(25,40) is what the `leppt_mt40` discriminant actually uses**. AUC (rel. to the (20,0) floor): m_T sweep 0.954 μ / 0.948 e at pT>20, 0.890/0.895 at pT>25, vs pT sweep 0.821/0.834 — **the m_T cut does the discriminating, the pT cut does not**. m_T>40 at pT>25 keeps 26% of the μ QCD (335.7/1282.0) and 34% of the e QCD (1173.6/3436.2); the multijet fraction of the modelled S+B falls 23%→7.5% (μ) and 49%→26% (e) — consistent with the ABCD data-fraction figures in `qcd_abcd.C`. Sideband stats 14239 μ / 29710 e events. NB an m_T precut suits the lepton-pT-fit path; it is NOT compatible with the MET-shape fit (it would remove the QCD-constraining low-MET region). **pT grid extended down to 20 on 2026-08-04** (skim scan planes filled for pT > `kPtScanFloor` = 20 — see Stage 1; `kPtRef = 25` keeps ε_S/baseline quoted vs the CURRENT selection, so ε_S > 100% = signal gained). WITH an m_T cut the FOM still rises as pT drops (μ +3.5% from (25,45) to (20,42.5)); WITHOUT one, lowering pT only hurts — μ 52.83→50.04, e 38.50→32.15 — so pT→20 is a small win for the μ lepton-pT path, a wash for e, and NOT advisable for the MET-fit selection. Standing caveat: the sideband-shape normalization `qcdTot/sideTot` is anchored over the extended pT>20 plane, so the pT≥25 QCD inherits an extrapolation across the 20–25 slice where the transfer factor differs (the iso–pT-correlation systematic, bounded by `antiiso_shape_pt_*`); ROC/AUC are normalized to the loosest point (20, 0), NOT comparable to pre-2026-08-04 AUCs; per-pT-cut optimum table (pT 20–30) printed in the report. **Written up in [docs/AN_selection_optimization.tex](docs/AN_selection_optimization.tex)** (+ compiled `docs/AN_selection_optimization.pdf`) — the AN subsection "Optimisation of the isolation and lepton-pT requirements", which covers BOTH this scan and the isolation/ID study, with the S/√(S+B) background definition spelled out and the figure-file mapping in its header comment.
- [correction/qcd_sideband_fit_and_extrapolate.C](correction/qcd_sideband_fit_and_extrapolate.C) — **superseded** by `qcd_abcd.C`. QCD shape from anti-iso sideband, Rayleigh-like fit + linear shape-parameter extrapolation to signal iso; did not behave at pO statistics. Kept for reference.
- [correction/isolation_mu_tight.C](correction/isolation_mu_tight.C) (muon, **current**), [correction/isolation_ele.C](correction/isolation_ele.C) (electron) — isolation ROC study, both on the **Δβ-corrected** relIso since 2026-07-06. Signal = OS Z-window pairs; backgrounds = SS pairs (≈empty for μ) and single-lepton MET<5 (QCD-like). `isolation_mu_tight.C` is the fresh, pure muon study (2026-07-06): ONE ID (`muIDTight`, pT>15, |η|<2.4), ONE variable (branch-based Δβ relIso = the skim's `RelIsoPF`), continuous 200-point scan, dbeta + nodbeta(reference) tags → `rootfile/IsoStudyOutputs_muon_tight.root`. **Electron `isolation_ele.C`** scans the 8 ID WPs × 4 MVA-Iso WPs AND a continuous-relIso scan per ID. **Run both through the wrapper [correction/run_isolation.sh](correction/run_isolation.sh)** (`./run_isolation.sh [mu|ele|both]`, default both) — it echoes the input path in use, pre-builds once (ACLiC race), tees the scan tables to `correction/logs/isolation_{mu,ele}.log` (ROOT prints them to stdout ONLY, and they are what `AN_selection_optimization.tex` quotes) and regenerates the summary plots; ~10 min per channel (full ntuple loop). **Input path single-sourced since 2026-08-20** — all three iso macros now `#include ../skim/skim_common.h` and default to `pOSkim::kDefaultDataFile` (`DATA_FILE=… ./run_isolation.sh` overrides), so they follow the production automatically; before that they defaulted to May-26 (`isolation_mu_tight`) or the long-dead `pO_2025.root` on EOS (`isolation_ele`/`isolation`). NB `isolation_ele.C`/`isolation.C` had to rename their local `ComputePFMET` → `ComputePFMET_ele`/`_isoleg`: `skim_common.h` declares `pOSkim::ComputePFMET` and `njet_WZ.C` does `using namespace pOSkim`, so the unqualified call there became ambiguous once both dictionaries were loaded in one session (`isolation.C` + `isolation_ele.C` together still clash on their other shared statics — pre-existing, they are near-duplicates, and nothing loads both). The older multi-cone/multi-def muon scan [correction/isolation.C](correction/isolation.C) is kept untouched for reference (uncorrected iso). **Production-independence VERIFIED 2026-08-20:** re-running both studies on July-29 reproduces the May-26 results **bit-identically** (μ: 372 OS pairs / 2 SS / 10373 QCD-like muons, AUC_MET 0.9425, J 0.8054@0.156, ε 0.9341/0.1300; e: all 8 AUCs to 4 digits). Reason measured directly: the two productions cover the same runs [393952, 394007] and contain **exactly 99959 tight muons with pT>15** each — July-29 adds 1.46M events of which **914k have no lepton at all** (May-26 has zero such events, i.e. it was lepton-skimmed), and none of the additions supply a tight muon above 15 GeV. **So the isolation working points were never stale in content**, only in provenance — but the macros were pointing at a file that could vanish, and a future production change would have silently failed to propagate. Control-sample sizes: the multijet proxy is **10373 μ / 19027 e** (at `eleMVAIdWP95`) leptons; NB the per-ID `contPassMET_*` graphs end at relIso<1.0 (μ 7882, e 17161), which is NOT the sample size — divide by `contEffMET_*` for the total. **Conclusions (re-derived under Δβ, 2026-07-06; unchanged on July-29):** μ TightID AUC_MET=0.943, J(QCD) optimum 0.156 ⇒ relIso<0.15 stands (ε_sig=0.934, ε_QCD=0.130; Δβ vs uncorrected nearly identical — pO pileup tiny, PU-iso nonzero for only ~20% of leptons); e `MVAIdWP95` AUC_MET=0.909, J optimum exactly 0.095 ⇒ cut stands (ε_sig=0.910, ε_QCD=0.184; `CutWP95` remains the worst ID, AUC 0.868). **Control samples (the "background" of these ROCs — NOTE the FOM here is efficiency-only, NOT S/√(S+B), because the two samples have no common normalization):** signal = OS Z-window pairs 80<m<100, counted PER LEG (μ 372 pairs / 744 legs, e 310 / 620); background = **exactly one ID'd lepton + PF MET < 5 GeV** (μ 7882, e 17161 leptons at MVAIdWP95); SS pairs are a cross-check only (μ 1 pair, e 5). The anti-iso sideband CANNOT serve here — relIso is the scanned axis. **Why the ID switch mattered, quantitatively:** ε_QCD at fixed relIso is similar across IDs (0.165–0.209), but the number of fake candidates the ID admits BEFORE isolation differs 2.7× (CutIdWP95 45759 vs MVAIdWP95 17161) — the rejection happens upstream, in the ID. Full 8-ID table (AUC / J / ε at 0.095 / fake candidates) in [docs/AN_selection_optimization.tex](docs/AN_selection_optimization.tex), which writes up this study together with the (pT × m_T) scan.
- [correction/PlotsIsoROC.C](correction/PlotsIsoROC.C), [correction/PlotIsoROC_ele.C](correction/PlotIsoROC_ele.C) — isolation ROC curves → `correction/plotsROC[_ele]/`; plus `plotsROC_ele/roc_MET_continuous_allID.png` (one-off overlay of the 8 IDs' continuous-relIso QCD-ROC). **`PlotsIsoROC.C` is a LEGACY-ONLY plotter** (corrected 2026-09-21 — this bullet used to claim "`PlotsIsoROC.C+(false)` for the tight set", which is wrong): it reads `IsoStudyOutputs[_soft].root` / `ggbranchStudyOutputs*.root` / `PtcutStudyOutputs*.root`, i.e. the outputs of the superseded pre-Δβ multi-cone `isolation.C`, and never touches `IsoStudyOutputs_muon_tight.root`. `false` selects the non-soft LEGACY set, not the tight study; the DEFAULT argument (`isSoft = true`) points at `*_soft.root` files that **do not exist** in `correction/rootfile/`, so only `PlotsIsoROC.C+(false)` runs at all. The current muon study is plotted by `plot_iso_summary.C` → `plotsROC/summary_muon.*` — the only non-legacy files in that directory (the other 54, dated 2026-06-24, are the superseded multi-cone scan).
- [correction/plot_iso_summary.C](correction/plot_iso_summary.C) — **per-fixed-ID summary in one consistent cosmetic**: a 2-pad canvas (LEFT = continuous efficiency vs Δβ-relIso cut for signal/QCD/SS with the operating cut marked + ε_sig/ε_QCD; RIGHT = the corresponding ROC + AUC with the operating-point star). Electron → one per ID (`plotsROC_ele/summary_<ID>.png`, reads the continuous scan from `isolation_ele.C`); muon → single tight-ID summary (`plotsROC/summary_muon.png`, reads `IsoStudyOutputs_muon_tight.root` from `isolation_mu_tight.C` — since 2026-07-06; before that the old `ggbranch` scan). `root -l -q 'plot_iso_summary.C+'`.
- [correction/dataMC_kinematics.C](correction/dataMC_kinematics.C) — Data vs signal-MC overlay + ratio pad, shape-normalized. **Z channels** (`"Zmm"`/`"Zee"`): DY MC vs data for the Z kinematics (`h_Zpt/h_Zeta/h_Zy/h_Zphi`, `h_lepPt/h_lepEta/h_lepPhi`, `hMass`) — the Data/MC boson-pT check. **W channels** (`"Wmu"`/`"Wel"`, added 2026-07-29): leading-lepton `h_lepPt/h_lepEta/h_lepPhi` after the full W selection, data vs W⁺+W⁻ MC combined with the relative `k_s` from `mc_norm.h` (shape-norm cancels the absolute scale; backgrounds NOT subtracted — and since the W selection has NO MET cut, QCD here is ≈19% μ / ≈50% e of the plotted sample (data−EWK, 2026-07-29; the familiar 3%/10% are high-MET-region figures), so expect a low-pT data excess — it's a shape check). Probes whether lepton pT is well-enough modeled to serve as an alternative fit discriminant (Jacobian peak, MET-free). **`"Wmu_mt40"`/`"Wel_mt40"` (2026-07-30)** run the same check on the pT-discriminant selection (pT>25 && m_T>40) via the `_mt40` histos → `plots/dataMC_W{mu,el}_mt40/`; with QCD down to ~5% (μ) / ~28% (e) there, the shape comparison is far less background-dominated. `root -l -q 'dataMC_kinematics.C+("Wmu")'` etc.; W data file has no sample suffix (`<base>_hist.root`).
- [correction/recoil_raw.C](correction/recoil_raw.C) — raw look at the Z→μμ hadronic recoil (`u_par`/`u_perp`) for the MET recoil correction: data vs DY MC, inclusive and in q_T slices (projected from the 2D histos), **plus a printed entry-count table per q_T slice** to judge statistics/binning before fitting. Slice edges = editable `kQtEdges` array (projections → no re-skim to change). **Checked 2026-07-02: recoil looks stable for now** — no correction applied, the planned `recoil_fit.C` (double-Gaussian per q_T bin → μ(q_T)/σ(q_T), AN Sec 6) is deferred; the MET data/MC discrepancy is attributed to the unembedded MC samples (no pO underlying event), not recoil. `root -l -q 'recoil_raw.C+'`.
- [correction/recoil_raw_ele.C](correction/recoil_raw_ele.C) — electron-channel twin of `recoil_raw.C`: same raw-recoil look for Z→ee (reads `ZToEE_pO2025_*`, writes `plots/recoil_Zee/`). `root -l -q 'recoil_raw_ele.C+'`.
- [correction/njet_WZ.C](correction/njet_WZ.C) — **jet multiplicity in W-/Z-tagged events** (2026-08-05). Loops the NTUPLES directly (skim outputs carry no jet info), replicating the skim selections exactly: W = the full 8-step selection (verified event-identical — data counts equal the skim cutflow N[8]: 5203 μ / 7065 e), Z = first iso-selected OS pair in the peak [60,120] (356 μμ / 248 ee; event counted once; its `RunZ` iso cut follows the skim — μ 0.20 → 0.15 on 2026-09-14). Jets = `ak4PFJetAnalyzer/t` (entry-aligned with EventTree, equal-entry-counts enforced; verified run/evt-identical on July-29), calibrated `jtpt` > 30 GeV, |`jteta`| < 2.1, cleaned by ΔR > 0.4 against the selected lepton(s) (both Z legs — essential for electrons, which always double as a PF jet); no jet ID (PF-composition branches exist if wanted). Data vs signal MC (W: Wp+Wm; Z: DY), gen-weighted, shape-normalized `SaveDataMCRatio` overlays; W channels also fill the `_mt40` njet twin (pT>25 && m_T>40 — data 4382 μ / 4465 e, matching the known composition counts). Outputs `rootfile/njet_<chan>.root` (per-sample raw + `*_mc` combined at the ABSOLUTE `k_s` scale — the combined totals close on the known expectations: Zμμ 372.1 vs data 356, Zee 252.5 vs 248, W μ 3821.7 / e 3178.7) and `plots/njet_<chan>/{njet[,njet_mt40],jetpt,jeteta}` (header lower-right, `ps.yHeadroom` opt-in added to `plotting_helper.C` for the log-pad top room — legacy 1.4 default untouched elsewhere). Second arg = plots-only mode: `njet_WZ.C+("all", true)` redraws every plot from the stored rootfiles in seconds (cosmetics iterations, no ntuple loops). **Result (2026-08-05):** Z njet is well modeled (⟨njet⟩ data/MC: μμ 0.183/0.160, ee 0.177/0.167, frac(≥1 jet) ≈ 14%); W data sits far above signal MC (μ 0.291/0.138, e 0.372/0.146) — that's the jet-rich QCD in the no-MET-cut W sample, shrinking under m_T>40 (μ 0.222, e 0.300); the remaining excess = residual QCD (not subtracted) + real W+jets vs the unembedded MC. `root -l -q 'njet_WZ.C+'` (all four channels) or `njet_WZ.C+("Wmu")` etc.

- [correction/charge_flip.C](correction/charge_flip.C) — **lepton charge-misID (charge-flip) rate from the W signal MC (2026-09-01)**, μ and e. Replicates the skim's 8-step W selection (the `njet_WZ.C` replication, event-identical with the skim) on the Wp/Wm files, then matches the selected leading lepton to the gen lepton **CHARGE-BLIND**: gen |pdg| 13/11, ΔR<0.5, |ΔpT|/pT_gen<0.5 = the AN's gen-reco criteria = `skim_common.h::PassGenRecoMatchingWithAncestor` MINUS its same-charge requirement (that shared helper is a no-op in the skim and CANNOT serve a flip study as-is — with the charge inside the match a flipped lepton simply fails to match). Matched + opposite charge = FLIP, matched + same = correct, no candidate = unmatched (reported, not in the denominator); f = N_flip/N_matched with Clopper-Pearson 68% intervals on RAW counts (gen weight ≈ constant; the weighted f is printed alongside). Binned in y = −η_lab (`kYEdges`, so bin i = analysis bin i) and pT {25,30,35,40,45,50,60,80,120}. Also prints the **asymmetry bias per bin: A_reco − A_gen on the same matched leptons** (Wp/Wm k_s-combined) = the uncorrected flip bias (= −2fA to first order). **Two ntuple facts (2026-09-01):** (1) the W is NOT stored in the filtered gen tree `HiGenParticleAna/hi` (0 |pdg|=24 entries in 100k events; every W lepton has `nMothers`=1 with `motherIdx`=−999 = "mother not in the stored list"), so `HasAncestor(…,24,…)` is ALWAYS false there — the W-ancestor preference here and in `gen_xsec.C` is a structural no-op (ΔR-nearest / highest-pT candidate used; `gen_xsec.C`'s fallback WARN fires for the same reason, 100% of events). (2) The EventTree carries a SECOND gen block `nMC/mcPID/mcStatus/mcPt/mcEta/mcPhi/mcMomPID/mcGMomPID/…` which DOES store the W (status 62, one per event) plus direct parentage; the ntuplizer's `mu_/ele_genMatchedIndex` indexes THAT block (99.96% of pT>25 reco μ → a status-1 μ with `mcMomPID`=24), NOT `hi` (read against `hi` it lands on the neutrino). It is CHARGE-AWARE (0 charge disagreements in 45k matched e ⇒ flipped leptons return idx=−1), so like the skim helper it cannot see flips. If W ancestry is ever needed, read it from `mcMomPID`, not from `hi`'s `motherIdx`. **Verified equivalence (2026-09-01):** the `hi` tree is COMPLETE within its filter (the mc* W electron is found in `hi` in 224,106/224,106 events where it has pT>5, |η|<2.5), and a charge-blind match against `hi` vs against `mc*` gives bin-identical flip counts for |η_reco|<2.4 (42/117/146 flips in |η| <1.5 / 1.5–2.0 / 2.0–2.4 with both lists). They differ ONLY at |η_reco|>2.5, where the ntuple still stores reco electrons up to |η|=3.0 (no tracker; f≈11% there) and `hi` has no gen particle — the skim's |η|<2.4 cut removes those, so the choice of gen list is immaterial for the analysis objects. **Results (July-29 MC, post-2026-08-18 selection): μ 1 flip in 1,053,500 matched ⇒ f < 3×10⁻⁶ (CP68 upper), ΔA = 0 — negligible. e f = (2.53±0.05)×10⁻³ inclusive, charge-symmetric (e⁺ 2.57 / e⁻ 2.49 ×10⁻³), ≈1×10⁻⁴ at |y|<0.8 rising to 8.4–9.2×10⁻³ at |y|>2.0 (material/brem), mild pT rise 2.2→5.2×10⁻³ (25–30 → 80–120 GeV); ΔA_e = −0.00064 inclusive (−0.6% relative), −0.0034 (−1.5%) in the outermost bins; unmatched 0.7% e / 0.001% μ. Electron charge definitions: `eleTrkCharge` alone flips 2.4× more (6.1×10⁻³); requiring `eleCharge == eleTrkCharge` drops 0.5% of electrons and removes 29% of the flips (f → 1.8×10⁻³); the inconsistent 0.5% flip at 14%.** Run through [correction/run_charge_flip.sh](correction/run_charge_flip.sh) (`./run_charge_flip.sh [mu|ele|both]`: pre-build + tee → `logs/charge_flip_<chan>.log` = THE record, headlines echoed; ~1 min per MC file). Outputs `rootfile/charge_flip_{mu,ele}.root` (raw count histos per sample + combined, TEfficiency objects, labeled counters — `charge_flip.C+("both", true)` redraws tables + plots from it without ntuple loops) and `plots/charge_flip_{mu,ele}/` (`fliprate_vs_{y,pt}`, `fliprate_2D_pt_y`, `unmatched_vs_y`, `match_{dR,dpt}`; e: `fliprate_vs_{y,pt}_chargedef`, `incons_vs_y`). Caveats: unembedded MC (no UE track confusion); the data-driven cross-check is the SS/OS Z→ee ratio (Z→μμ has ~1 SS pair). Not yet propagated to the fit (a per-bin correction/systematic of ≤1.5% on A_e, if wanted).

- [correction/trig_eff_mb.C](correction/trig_eff_mb.C) + [correction/run_trig_eff_mb.sh](correction/run_trig_eff_mb.sh) — **single-lepton TRIGGER EFFICIENCY (turn-on) of the W analysis path from the MINIMUM-BIAS-triggered sample, data vs W signal MC on one plot (2026-09-09 — the first step of the SF-application phase, per the group's definition):** `den` = the W selection (skim steps 1,2,4,5,6,7 — i.e. WITHOUT step 3 trigger and step 8 match, the two things being measured; replicated as in `charge_flip.C`/`njet_WZ.C`, event-identical with the skim) **&& `HLT_MinimumBiasHF_OR_BptxAND_v1` fired**; `num` = den && the analysis path (`HLT_OxyL1SingleMuOpen_v1` / `HLT_OxyL1SingleEG10_v1`) fired && the leading lepton matched (ΔR < 0.4, the skim's step 8). **Since 2026-09-15 the MC `den` additionally requires the leading lepton to gen-match a PROMPT W LEPTON** (`MatchGenLeptonFromW`: |pdg| 13/11, `mcStatus` 1, `|mcMomPID| = 24` or FSR via `mcGMomPID`, **charge-blind**, ΔR < 0.5, |Δp_T|/p_T < 0.5 = `charge_flip.C`'s window) — user decision, and **the ONLY difference between the data and MC legs**. It reads the EventTree's OWN `mc*` block, never `HiGenParticleAna/hi`, whose `motherIdx` is −999 for every W lepton (verified on `July_29_MC_Wp_mu`: one |pdg| = 24 entry per event, 99.98% of gen muons with pT > 20, |η| < 2.4 carry `mcMomPID` = ±24 directly, all status 1). Charge-blind so a charge-misidentified lepton still counts as prompt — charge misID is not a trigger inefficiency and folding it in would double-count `charge_flip.C`. **Measured effect: 16/754039 (W⁺) and 27/640965 (W⁻) selected events rejected = 0.002–0.004%; no quoted number moved** (inclusive SF 0.9971 before and after), exactly as expected since the samples are W→μν by construction. It self-documents the definition and will matter for electrons (0.7% unmatched). **NB it does NOT symmetrize the two legs** — the DATA denominator is still a mixture (~5% QCD fakes with m_T > 40, ~19% without), which is the real asymmetry and the reason `mt40` is nominal; matching MC to the data's composition would mean adding fakes to MC, not gen-matching it. A missing `mc*` branch in an MC file is FATAL. The lepton-pT floor of steps 1/6 is lowered to 10 GeV so the turn-on is visible (25 is a bin edge; the analysis selection = the pT > 25 sub-range, same object). Two variants in one loop: `nom` (plain W selection) and `mt40` (+ m_T > 40 from PF MET, read only for selected events) = the leppt_mt40 selection. **MB availability (verified 2026-09-09):** `hltanalysis/HltTree` stores exactly 5 paths in data AND every MC file (the 4 lepton paths + MB; the `_PrescaleNumerator/Denominator` columns are always 1/1 = uninformative); there is no `hltobject/` tree for MB (nothing to match — "MB matched" ≡ the bit fired). In DATA the MB bit is a **run-dependent prescale**: fraction of lepton-triggered events with MB = 393952 100%, 393953 91%, **393974/393975 0%** (924k events = 24% of the data, no MB at all), 393976 67%, 394004 51%, 394005/6 68%, 394007 78%; 49–50% overall, identical for μ- and e-triggered events within each run and flat in lepton pT (the `mbfrac_pt_*` control: all-selected vs triggered MB fractions 0.5018/0.5006 μ, 0.5020/0.5023 e) ⇒ a pure prescale, uncorrelated with the lepton — it costs statistics (den = 3065 μ / 4147 e of 6108 / 8261 selected), nothing else. In MC MB fires on 99.93% of W events (no prescale). 390k data events fire none of the 5 stored paths (other IonPhysics0 paths). **Results (pT > 25, |η| < 2.4, CP 68% on raw counts; MC = Wp+Wm raw, the k_s-weighted and gen-weighted combinations agree to 1e-4):** μ data 0.9889 −0.0022/+0.0019 vs MC 0.9927 → **SF 0.9962 ± 0.0021** (mt40: 0.9898 vs 0.9927, SF 0.9971 ± 0.0022); flat in pT from 10 GeV (L1 open muon has no threshold), MC dips to 0.984 at |y| < 0.4 (η≈0) with data following (SF 0.97–0.98 there, ≈1.00 elsewhere), charge-symmetric (MC W⁺ 0.9923 / W⁻ 0.9931 — the 2026-08-12 ΔR fix confirmed; data 0.9865/0.9919). e data 0.9626 −0.0032/+0.0030 vs MC 0.9735 → **SF 0.9889 ± 0.0032 on the plain selection, but 0.9705 vs 0.9736 → SF 0.9969 ± 0.0036 with m_T > 40** (CURRENT values after the 2026-09-14 ECAL-gap veto and the 2026-09-15 gen match: data 0.9630 vs MC 0.9732 → **SF 0.9895 ± 0.0032** plain, 0.9709 vs 0.9733 → **SF 0.9975 ± 0.0037** mt40 — the gen match rejects 0.70–0.75% of MC events, 100× the muon rate, but leaves the SF unchanged to 4 decimals because those events have ~the same trigger efficiency as the kept ones; the electron |y| dependence is present but weaker, p = 0.020 vs the muon's 0.0006, and there 3 |y| bins would suffice, p = 0.22): the plain-selection deficit is the **QCD-fake composition of the data sample** (≈50% fakes without the m_T cut, ≈27% with — extrapolating to zero fakes gives ε_W,data ≈ 0.980, ε_fake ≈ 0.946), i.e. the MB method measures the MIXTURE's efficiency in data; the MC e efficiency rises slowly 0.956 (25–27.5) → 0.994 (80–120) and dips to 0.955 at |y| < 0.8, data agreeing within stats after m_T > 40; the true EG10 turn-on sits at 10–16 GeV (MC 0.67/0.85/0.93), well below the cut; the ΔR match costs 0.16% (μ) / 0.4% (e) in data on top of the bit; per-run efficiencies stable (μ 0.988–0.989, e 0.95–0.97; run 393952 e 0.916 on 119 events = 2.5σ, not significant). **vs tag-and-probe:** this is the per-EVENT trigger efficiency of the W selection itself (W kinematics, charge mix, the real fake content, no factorization assumption, no tag bias) — the ε in σ = N/(L·A·ε) for this selection — at the cost of the MB prescale and the data-sample composition; T&P on Z→ℓℓ gives the per-LEPTON efficiency on pure prompt leptons with Z kinematics (388 μμ / 284 ee events here), applied to MC per lepton assuming factorization. Roles: T&P → the SF on MC; MB method → closure of the SF-weighted MC event efficiency + the only handle on the fake-composition effect. No T&P code exists yet. Run through the wrapper (`./run_trig_eff_mb.sh [mu|ele|both]`, pre-build + tee → `logs/trig_eff_mb_<chan>.log` = THE record incl. the per-pT/per-y tables, the inclusive block, the MB-prescale checks and the data per-run table; ~70 s per flavour — the data pass reads only ~15 branches). Outputs `rootfile/trig_eff_mb_{mu,ele}.root` (raw count histos per sample + `_mc` sums, TEfficiency objects, `sf_pt/sf_y[_plus/_minus]` graphs **+ the |y|-folded twins `h_{all,den,bit,num,trg}_absy_<sel>[_plus/_minus]` and `sf_absy_<sel>[_plus/_minus]`, 6 bins of 0.4 — `skim/muon_sf.h` reads these (and has since 2026-09-15), but since 2026-09-21 only to PRINT them as the measured-but-not-applied record; it applies the inclusive SF**; `trig_eff_mb.C+("both", true)` redraws from it) and `plots/trig_eff_mb_{mu,ele}/`: `turnon_pt_<sel>[_zoom25]` (data vs MC + Data/MC pad), `eff_pt_<sel>_charge`, `eff_y_<sel>[_charge]` (pT > 25), **`eff_absy_<sel>[_charge]` — the rapidity-dependent SF: its ratio pad IS the per-|y| SF that is measured and deliberately not applied, i.e. THE plot that justifies the inclusive choice**, `bit_vs_match_pt_<sel>` (path fired vs fired+matched), `mbfrac_pt_<sel>` (the prescale control), `trig_eff_<sel>.csv` (per-bin den/num/eff/SF).
  **BINNING DECISION block (2026-09-15) — the record on whether and how the SF should be binned. Still printed on every run, and since the 2026-09-21 revert to an inclusive SF it is the JUSTIFICATION document rather than the recipe: it says the η dependence is real, and the analysis chooses not to correct for it.** Every run prints likelihood-ratio tests of "one flat SF" against each candidate binning (`FlatnessLRT`: k_i ~ Binomial(n_i, SF·ε_MC,i), ε_MC exact since MC denominators are ~200× the data ones; −2ΔlnL ~ χ²(nbins−1)). mt40 / nom: **pT 2 groups p = 0.71/0.80** (no pT dependence) · **|y| 6 bins p = 0.0006/0.0004** (NOT flat, ~3σ) · **|y| 3 → |y| 6 p = 0.0043/0.0132** (3 coarse bins are not enough — the dip is confined to |y| < 0.4) · **|y| 6 → signed y 12 p = 0.39/0.14** (folding justified, no left/right asymmetry) · **|y|3 → |y|3 × pT2 p = 0.68/0.68** (no pT×η interaction ⇒ 1D, not 2D) · per-charge |y| SFs 0.9969 vs 1.0003 (charge-inclusive; the interaction test gives p = 0.13). **This test REPLACED a Clopper-Pearson pull χ², which toys showed is ~16× under-powered** (⟨χ²⟩ = 3.15 for 5 dof, 0.3% rejection at a nominal 5%, because CP intervals over-cover as ε → 1) and had produced the 2026-09-14 "consistent with one flat SF, p = 0.41" conclusion. The LRT is calibrated on the same toys (⟨−2ΔlnL⟩ = 5.48 vs 5.00, 7.2% vs 5%); the toy-calibrated p for the observed 19.78 is 0.0020. The per-bin pT scan (11 bins) prints p = 0.087/0.0054 but answers a different question and is driven entirely by the [50,60) bin (SF 0.964/0.962 on ~150 events, the same ~2σ downward fluctuation in both selections, no trend) — quote the coarse test.

- [correction/idiso_sf_skim.C](correction/idiso_sf_skim.C) + [correction/idiso_sf_inputs.C](correction/idiso_sf_inputs.C) + [correction/idiso_sf_common.h](correction/idiso_sf_common.h) + [correction/run_idiso_sf.sh](correction/run_idiso_sf.sh) — **the ELECTRON ID+ISO SCALE-FACTOR CROSS-CHECK (2026-09-24/25), a SEPARATE stream: its own skim, inputs, fit (`test/run_pO_idiso_sf.sh` in the fork) and outputs; `skim/skim.C` and the nominal inputs are untouched.** Question (user): can the pp EGM `wp90iso` SF (one SF for ID + isolation, relative to a reco electron) be applied if the electron selection switches to the forest's `eleMVAIdWP90 && eleMVAIsoWP90`? The forest classifiers are not EGM VID (see the `project_electron_sf` memory), so only data/MC RATIOS are compared. Common electron: pT > 25, |η_SC| < 2.4, no crack, **relIso < 0.3** (frees 0.3–1.0 for the QCD templates; 99.3% of prompt W electrons pass it); pass = WP90 ID && WP90 iso. Everything (definition, bins, names) single-sourced in `idiso_sf_common.h`. `./run_idiso_sf.sh [skim|inputs|all]` → `logs/idiso_sf_ele_{skim_<sample>,inputs}.log`; the EGM table comes from `skim/sf/extract_electron_sf.py` → `skim/sf/electron_sf_2025Prompt_{wp90iso,wp90noiso}.csv` (Python stdlib, kept as Python — user 2026-09-25). **Two measurements:** (1) **W sample** (user's design, like `trig_eff_mb.C`): the leading trigger-matched electron of W candidates, ε = N_W(pass)/N_W(pass+fail) from POST-FIT W counts of pass/fail channels (per charge × coarse bin: incl, |η_SC| 4 bins, pT 3 bins — EGM's edges merged) + the Z_PP peak as the DY anchor. **QCD = the NOMINAL IN-FIT ABCD** (user 2026-09-25: "always the in-fit ABCD, so we don't assume the prefit normalization" — the version of that morning, a sideband template held by an lnN around a prefit normalization, was wrong; its MET/m_T inputs and plots are in `rootfile/idiso_sf_ele_pre_infit/` and `plots/idiso_sf_ele_pre_infit/`): the fitted variable is the **lepton pT at m_T > 40** (the nominal `leppt_mt40`; MET/m_T cannot be fitted, m_T is the ABCD axis), per W channel the 1-bin counting CRs CRB (the category, m_T < 30), CRC (the relIso sideband, m_T > 40 — also the SR QCD pT shape) and CRD (sideband, m_T < 30) with free scales sB/sC/sD, the SR QCD scaled by the formula sB·sC/sD, the CRB W riding the channel's r, and the residual lnN 1.15 on the SR QCD only; sideband window 0.3–1.0 / 0.3–0.6 / 0.6–1.0 ⇒ 9 variants; ε therefore refers to W electrons with m_T > 40. **Does not work**: in the fail SR the W is 0.4–2.7% of the ABCD QCD A0 while the prefit closure (data − EWK)/A0 runs 0.74–1.20 across bins (`[ABCD]` lines; data 6,059 pass / 133,692 fail overall). On data the fit pulls the fail QCD down ~5% (θ ≈ −0.4σ) and puts r_fail at its 0 boundary (SF 1.24 inclusive and in every pT bin); the 0.6–1.0 window gives SF 0.41 (multiplier 0.73); the |η| fits do not converge. The **expected** precision (Asimov) is SF ± 0.14 inclusive, ±0.13–0.22 per |η| bin, ±0.18–0.71 per pT bin — r_fail = 1 ± 1.0 (local RooFit stand-in of the Combine likelihood). (2) **INCLUSIVE Z TAG-AND-PROBE** (user 2026-09-25): tag = a passing leg matched to an EG10 object; the OS pair with ≥ 1 tag closest to m_Z; event categories Z_PP (both pass = 2 passing probes) / Z_PF (1 failing probe), ε = 2N_PP/(2N_PP + N_PF); **its own fit** (`combine_input_idiso_ztnp.root`, POIs r_ZPP / r_ZPF, flat background with free normalization) — in ONE likelihood with the W channels the shared DY scale let them pull it (r_Z 1.10–1.21, SF_Z 0.96–1.01 across W variants whose Z channels are identical). Data 201 PP / 109 PF (SS 1 / 5). **SF = 0.959 ± 0.025** (stand-in fit; cut-and-count SS-subtracted 0.956 ± 0.025; exact per-pair probes 0.960) **vs the EGM prediction for the same probes 0.956 ± 0.009** (`[EGMPRED]`: ⟨SF_EGM⟩ over the passing DY-MC probes; per W coarse bin the prediction is flat, 0.94–0.96) ⇒ the pp SF describes our WP90 ID+iso in pO at the inclusive 2.6% level. Outputs: `rootfile/idiso_sf_ele/{skim_<sample>,combine_input_idiso_<tag>}.root` (+ `_meta.txt`, read by the fork's card generator), `plots/idiso_sf_ele/prefit/`. **COMBINE FITS DONE 2026-09-26 (lxplus, all 10 fits status 0 / covQual 3, every Asimov closure PASS):** Z tag-and-probe **SF = 0.9588 ± 0.0251** (expected ± 0.026; r_ZPP 1.085 ± 0.077, r_ZPF 1.352 ± 0.141, Z_PF postfit χ²/ndf 0.88) vs EGM 0.956 ± 0.009 for the same probes — pull +0.11. W: SF 0.41–1.24 inclusive depending on the sideband window — **and BOTH ends are the limits of the POI range r ∈ [0, 10], not measurements** (corrected 2026-09-27: 0.3–1.0 and 0.3–0.6 put the fail r at 0 ⇒ SF = 1/ε_MC = 1.24; 0.6–1.0 puts it at its UPPER limit 10 ⇒ 0.41, which the first plots drew as a measured point with ±0.006); of the 24 per-bin W results (3 incl + 12 |η| + 9 pT) only 5 have an interior fail r, and those sit at 2–9× the MC-predicted fail W (SF 0.44–0.80); 0.75 vs 1.24 for the SAME data depending only on the binning scheme (|η| vs incl/pT combined), per bin 0.36–1.31; expected precision ±0.11–0.14 combined, ±0.12–0.67 per bin; the fail-channel postfit χ²/ndf 46 (inclusive W⁺) with ±20σ pulls — the sideband QCD pT spectrum is far too soft ((data − EWK)/QCD 0.6 at 25 GeV → 2.5–3 above 80 GeV, the electron fake-factor pT rise); the W channels drag the DY scale r_Z to 1.10–2.27 (2.27 in `abseta` 0.6–1.0). ⇒ the W route is a documented no-go; the Z tag-and-probe says the pp SF describes WP90 ID + iso in pO at 2.6%. **Plots:** [correction/idiso_sf_plots.C](correction/idiso_sf_plots.C) (`./run_idiso_sf.sh plots`, reads the downloaded `test/pO_idiso_sf_out/` via `$FORK_TEST` + the skim EGM maps; the `[expected]` Asimov precision is parsed from the extraction logs) → `plots/idiso_sf_ele/`: `sf_summary` (all 10 measurements on one axis with fitted errors, expected-precision boxes and the EGM bands), `sf_egm_{abseta,pt}` (the EGM overlay: top = our W SF for the 3 windows + the EGM prediction per coarse bin + the Z band; bottom = a zoom with the EGM fine-binned SF per slice), `abcd_multiplier_{incl,abseta,pt}` (κ^θ·sB·sC/sD per channel and window), `fail_qcdshape_{full,sblo,sbhi}` (the inclusive fail SR vs the in-fit ABCD QCD template); markers: window 0.3–1.0 / 0.3–0.6 / 0.6–1.0 = circle / square / diamond, OPEN = a fail r of the bin (either charge) at a limit of its range, 0 or 10 (no fit error drawn; the CSV `boundary` column = no|r0|r10|r0+r10 — until 2026-09-27 only the 0 side was flagged); + `idiso_sf_summary.csv` (every number drawn, with the pull against EGM using the expected error). The EGM table helpers (`LoadEgm`/`FindEgm`/`EgmPrediction`) moved into `idiso_sf_common.h` (the inputs record regenerated byte-identical); `plotting_helper.C` gained an opt-in `PlotStyle::ratioLo/ratioHi` for `SaveDataMCRatio` (default 0.5–1.5 unchanged; fork copy re-synced).

The shared ratio-pad helper `SaveDataMCRatio(...)` lives in [plotting/plotting_helper.C](plotting/plotting_helper.C).

### `merge_rootfile/` — input prep

Utilities to scan EOS and `hadd` raw ntuple files into the **per-sample** ROOT
files that `skim/run_all.sh` reads (now the Aug-20 embedded production — see
"Input data"; earlier it was a single consolidated `pO_2025.root`):
- [merge_rootfile/make_filelist.sh](merge_rootfile/make_filelist.sh) — discover files on EOS (ONE sample dir → one .txt; the original, still used ad hoc)
- [merge_rootfile/hadd_from_list.sh](merge_rootfile/hadd_from_list.sh) — merge one list into one ROOT file
- **[merge_rootfile/make_filelists_Aug_20.sh](merge_rootfile/make_filelists_Aug_20.sh) (2026-09-21)** — the per-sample, self-driving replacement for calling `make_filelist.sh` nine times: walks the production dir, maps each `HiForest_<token>_*` directory to the canonical output name via one `SAMPLES` table (`DYToMuMu_M_50|DY_mu_Z`, `WpToENu|Wp_ele`, …) and writes `Aug_20_MC_*.txt`. The DY tokens carry the explicit `_M_50` so they cannot match the `M_10_50` dirs (deliberately excluded — the analysis has never used low-mass DY). Prunes CRAB `failed/` and `log/` subtrees, requires the sample dir to be UNIQUE, and **WARNs on duplicated file basenames** across CRAB submission dirs — a resubmitted task would otherwise be `hadd`ed TWICE with no error anywhere downstream. `PREFIX=` overrides the naming.
- **[merge_rootfile/run_all_hadd_Aug_20.sh](merge_rootfile/run_all_hadd_Aug_20.sh) (2026-09-21)** — Aug-20 twin of `run_all_hadd_July_29.sh`: globs `Aug_20_*.txt`, skips already-merged outputs (`FORCE=1` to redo), per-sample logs under `logs_hadd/`, one failure does not abort the rest.
- [merge_rootfile/run_all_hadd_July_29.sh](merge_rootfile/run_all_hadd_July_29.sh), `run_all_hadd.sh` — the previous productions' drivers
- `MC_*.txt`, `DATA_pass_*.txt`, `version_*.txt` — per-sample filelists

### `docs/` — write-ups (since 2026-09-27)

Every generated document lives here (user request):
- the AN-style LaTeX notes and their compiled PDFs: `AN_qcd_background` (the full trail), `AN_qcd_background_infit`, `AN_selection_optimization`, and `AN_lhe_systematics` (+ `AN_lhe_systematics_wrapper.tex`), moved from the repo root;
- the HTML write-ups, one folder per page with its images: `idiso_sf_crosscheck/` = the electron ID+iso SF cross-check (published at https://claude.ai/artifact/S3MkFm8dUnUZNKFf1SK5mG).

The `.tex` files are behind the analysis and do NOT compile as they stand: they reference a `figures/{qcd,opt}/` folder that no longer exists, so the PDFs are the record. `*.tex`, `*.pdf` and `*.png` are gitignored, so most of the folder is local only. The reference documents (AN2013_136, AN2017_058, HIN-17-007) stay in the repo root, as do the tracked `AN_outline.md` and `abstract_QM2027.md`. New write-ups go here; to update a published page, publish this copy to its URL.

## Input data

**Per-sample** ROOT files, not a single consolidated file. Single source of
truth: `kDefaultDataFile` (data) and `ResolveMCSample` (MC) in
`skim/skim_common.h` (`run_all.sh` greps the data path; override per-run via the
`DATA_FILE` env var).

**MC = the Aug-20 EMBEDDED production since 2026-09-21** (POWHEG + **Angantyr**
underlying event = the first production with a pO UE embedded in the
simulation; HiForest `2026_08_20`, MC ONLY, no low-mass DY). **DATA is still
the July-29 reconstruction** — that production ships no data. Everything sits
under `~/pO_2026_Aug_20/` (absolute paths; ROOT's `TFile::Open` doesn't
reliably expand `~`): 9 MC files `Aug_20_MC_DY_{mu,ele,tau}_Z.root`,
`Aug_20_MC_W{p,m}_{mu,ele,tau}.root`, plus `July_29_DATA.root` as a **HARD
LINK** to the July-29 copy (`ln`, no `-s`: same inode, link count 2, no second
45 GB — the volume had only 41 GB free, so a real copy would have failed; both
production directories stay self-contained and deleting either name keeps the
data). To read off EOS/lxplus again, repoint `inputBase` at
`root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2026_Aug_20/`
— but the DATA path needs the July-29 prefix there (no hard link on EOS).
Merged from EOS filelists by `merge_rootfile/make_filelists_Aug_20.sh` (walks
the production dir, writes one `Aug_20_MC_*.txt` per sample, prunes CRAB
`failed/`+`log/`, WARNs on duplicate basenames = a resubmitted task that would
be merged twice) + `merge_rootfile/run_all_hadd_Aug_20.sh` (glob-autodetects,
skips already-merged, `FORCE=1` to redo) — the Aug-20 twins of the July-29
scripts. Previous productions: July-29 (2026-07-30..2026-09-21) at
`~/pO_2026_July_29/`, May-26 at `~/pO_2026_May_26/`.

**VERIFIED BEFORE THE SWITCH (2026-09-21), and worth re-checking on any new
production:** all 7 trees present and entry-aligned (1.16–1.20M per file —
`njet_WZ.C` requires equal entry counts), `hltobject/` carries the 4 Oxy paths,
`ttbar_w` is **217** long with the July-29 layout, and **⟨weight⟩ = 6376.1 /
5463 / 1174.5 reproduces σ = 6.376 / 5.464 / 1.175 nb**, so `mc_norm.h` needs
no change. `|ttbar_w[0]|` is the same constant as July-29 (5538.8 for W⁻),
i.e. **the SAME generated events re-reconstructed with embedding** — confirmed
at gen level: `gen_xsec.C` gives σ_fid 49.80 / 39.63 / 8.769 nb vs July-29's
49.8 / 39.6 / 8.757, and `max|member0 − nominal| = 0` still holds. So every
downstream change is RECONSTRUCTION. N_gen moved 0–6% per sample (different
event counts) — irrelevant to yields, which are N_gen-invariant.

**N_gen labels stay canonical** (`Wp_mu`, `DYee`, `Wp_tau`, …):
`count_ngen.C::LabelFromFname` strips the dated prefixes (`Aug_20_MC_`,
`July_29_MC_`, `MC_`) so every `pONorm::MCScale(label)` call site is
production-independent.

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

### Current pipeline — `test/run_pO_fits.sh` (grand simultaneous fit since 2026-08-04)

Run under `cmsenv`:
```bash
./run_pO_fits.sh [mu|ele|both] [simfit|flavfit|all] [--dry-run] [--no-postfit] [--draw-only] [--asimov] [--no-statonly] [--no-contour] [--extract-only]
```

**`flavfit` (2026-09-22, user request: "a dedicated muon-only / electron-only
fit, still a sim fit but not μ and e combined", with the SAME treatment as the
grand fit) = the PER-FLAVOUR SIMULTANEOUS FITS.** The simfit model below with
ONE lepton flavour per likelihood: that flavour's 24 W channels + its own Z
peak (+ its 6 ABCD CR channels) per binning variant, 25 POIs (`r_<C>_y<i>` +
`r_Z`, measured by that flavour alone — the μ-only and e-only r's are
independent; the grand fit forces one r on both), and everything else
identical: the nuisances that act on the flavour (lumi, its two
`qcd_rate_<F>_*` rows, nPDF/qcdScale/alphaS, `muSF` in the muon fit only), the
three passes (nominal, `--statonly` companion = the stat error, `--contour`),
`--asimov`, the same extraction. `both flavfit` = one μ fit then one e fit;
`mu flavfit` one; mode `all` = grand + both flavfits. Work dirs
`pO_fit_out<suffix>/simfit_mu/`, `simfit_ele/` (the grand fit's layout), summary
files `simfit_<flav>_{W_yields.csv,summary.csv,fitted_yields.root}`, extraction
log `summary/extract_simfit_<flav>.log`. Mechanics: `run_simfit_set "<flavours>"`
parametrizes the unchanged simfit functions (work dir, flavour list, file tag);
`make_pO_simfit_cards.sh` takes `SIMFIT_FLAVS` and writes the flavour half of the
grand card (only that flavour's channels, QCD rows and the shape rows ITS
sidecars list — rows that would be all '-' are not written; map regexes
`mu_…` instead of `(mu|ele)_…`; sidecar line `flavours mu`), and
`extract_pO_simfit.C` takes trailing `flavours` / `outTag` args (only that
flavour's input is opened, S sums over the fit's flavours, the absent flavour's
two qcd CSV columns are 0,0). **Validated 2026-09-22 without cmsenv:** the grand
cards + maps + σ-POI maps regenerated BYTE-IDENTICAL (both discs; the sidecar
only gains `flavours mu,ele`); every row of the per-flavour cards
column-aligned (162 = 324/2 columns in abcd); the fork's own `MultiSignalModel`
(real RooWorkspace behind a minimal modelBuilder) scales every per-flavour
column by exactly the grand card's parameter, the σ-POI maps give 25 POIs,
r_i = 1 at init and Σ r_i G_i = sigmaW to 1e-16; the new extractor reproduces
the lxplus `comb_W_yields.csv` / `comb_summary.csv` of the 09-21 fit BYTE for
BYTE (and all 104 histograms) with the current inputs, and its per-flavour path
passes the prefit-S check. **FITTED on lxplus 2026-09-22** (`both flavfit
--asimov` + `run_pO_impacts.sh --fit all`; results in the flavfit
Current-state bullet). Consumers: the
analysis repo's μ-vs-e chain (`run_observables.sh`, plots/flavfit/).
**The legacy per-flavour per-bin pipeline (`perbin|incl|combined`) was REMOVED
the same day (user: "we can remove the old mu and ele option")** — those modes
exit with a pointer to `flavfit`; `make_pO_datacards.sh`, `extract_pO_yields.C`
and `make_yields_from_csv.C` were deleted (git history).

**`simfit` is the DEFAULT mode (2026-08-04) — the GRAND SIMULTANEOUS FIT.**
One likelihood per binning variant (lab, fb — same events rebinned, so fitted
separately): all 48 W channels ({mu,ele} × {Wp,Wm} × y0..11, channel names
`<F>_<C>_<B>_y<i>`) + BOTH Z peaks (`mu_Z_incl`, `ele_Z_incl`). 2N+1 = 25 POIs:
**`r_<C>_y<i>`** (24) scales that bin's W-related MC (`signal`+`wtau`) in BOTH
flavours' channels (μ/e SHARED — lepton universality; relative μ/e acc×eff from
MC, lepton SFs still not applied), and **`r_Z`** (the "+1") is ONE global scale
on all DY-related MC (`z`/`ztau` in every W channel + `zsig`/`ztau` under both
peaks — DY rapidity dependence trusted from MC, normalization pinned by the
peaks). **QCD: lnN-constrained since 2026-08-17** — one log-normal nuisance
per (flavour, charge), `qcd_rate_{mu,ele}_{Wp,Wm}`, correlated across that
flavour+charge's 12 y bins and uncorrelated across flavour/charge (T + template
are measured once inclusively in y; no unearned cancellation in the asymmetry).
κ **derived from the 2026-08-17 qcd_abcd printout** (stat ⊕ m_T-vs-MET-plane-T
transport ~10% ⊕ anti-iso tilt): **μ 1.15, e 1.20**; env-overridable
(`QCD_LNN_MU`/`QCD_LNN_ELE`, for robustness scans) and `QCD_MODE=free` restores
the pre-2026-08-17 48 free `qcd_norm_<channel>` rateParams.
**DEFAULT FLIPPED 2026-09-15b (user decision): `QCD_MODE` now defaults to
`abcd`** — but DISC-DEPENDENTLY, because abcd exits 2 on anything but
leppt_mt40 (the `qcd_abcd` template and the CR dirs live only in
`combine_input_W_leppt_mt40.root`) and a flat default would have taken the
whole simfit down on `--disc met`. So: **abcd for leppt_mt40, lnN for
met**, decided in `make_pO_simfit_cards.sh` (which knows `DISC`) and
echoed as `QCD_MODE not set -> default '<mode>' for disc '<disc>'`; an explicit
`QCD_MODE=` always wins, so `QCD_MODE=lnN --disc leppt_mt40` is still the
like-for-like comparison. This closes the "lnN→abcd default flip" decision that
had been open since 2026-08-24. **STILL OPEN and relevant to the same refit:
`QCD_ABCD_LNN_ELE` is 1.15, but the 2026-09-14 ECAL-gap re-run gives a reduced
κ_e of 1.12** — left as it was, deliberately, rather than changed unasked.
**`QCD_MODE=abcd` (2026-08-23, leppt_mt40 ONLY — hard error otherwise): the
IN-FIT ABCD.** 12 counting CR channels `<F>_<C>_CR{B,C,D}` join the likelihood
(imax 62; 1-bin shapes from the input file's CR dirs — CRB iso-pass m_T<30,
CRC anti-iso m_T>40 = the SR template's source region, CRD anti-iso m_T<30);
free scales `qcd_s{B,C,D}_<F>_<C>` (init 1, [0,10]) float the three prefit QCD
counts, and the SR qcd — the **`qcd_abcd`** template, total B0·C40/D0 — is
scaled by the FORMULA rateParam `qcd_abcd_<F>_<C> rateParam <F>_<C>_<B>_y* qcd
(@0*@1/@2) sB,sC,sD` (Combine's documented ABCD pattern,
docs/part2/settinguptheanalysis.md; a glob wildcard, matching only that card's
12 SR channels). No baked constants: prefit ≡ the ABCD prediction and Asimov
closure demands all scales = 1. The EWK subtraction rides the POIs: CRB
`z`/`ztau` are swept into r_Z by the existing catch-all map; the CRB W content
is 12 per-y processes `w_y0..11` (← histograms `w_<B>_y*` — the lab/fb wiring
is load-bearing) mapped to `r_<C>_y<i>` via 24 new map lines (49 total;
`QCD_WCR=frozen` switches to the single frozen `wfix` with no extra maps —
residual ~2%/1% then belongs in κ); CRC/CRD carry one frozen `ewk` (≤2.5% of
those regions). The 4 `qcd_rate_*` lnN rows STAY but with the REDUCED κ
(`QCD_ABCD_LNN_MU`=1.09 / `QCD_ABCD_LNN_ELE`=1.15 from the 2026-08-23
qcd_abcd log: window ⊕ FF-shift; on the SR qcd columns ONLY — never on CR qcd,
whose yields are measurements). Sidecar gains a `qcdMode <abcd|lnN|free>` line
(3-line legacy sidecars still inferred from the κ's); the driver passes it to
the extractor as an 8th argument. Extraction: the formula is a RooFormulaVar —
ABSENT from `floatParsFinal` — so `extract_pO_simfit.C::AbcdScale` evaluates
M = κ^θ·sB·sC/sD from the floating constituents with the full-correlation
error propagation (g = M·(lnκ, 1/sB, 1/sC, −1/sD)); M lands in the same CSV
qcd columns ("factor on the input-file qcd template", which is now qcd_abcd)
plus a new 18th column `qcd_model`, the scales + `qcd_abcd_mult_*` go to
comb_summary.csv, and the Asimov closure additionally requires the 12 scales
= 1. **Local verification 2026-08-23 all-green**: lnN/free cards byte-identical
to the pre-change baseline; abcd card imax 62, 386 shapes lines, all rows
column-aligned (324 process columns), correct lab/fb `w_lab_y*`/`w_fb_y*`
wiring; postfit_incl + run_observables run on old CSVs unchanged. **FIT RUN
on lxplus 2026-08-24: Asimov closure EXACT (<1e-7), status 0 / covQual 3;
multipliers μ 0.894/0.919 (sB≈0.94 = the floating subtraction working),
e 0.833/0.654 (θ_e⁻ −2.6σ — see the postfit diagnosis in the active-work
bullet); r_Z 1.095±0.053. `sync_lxplus.sh download` now also pulls
`datacards/` (cards+maps+sidecar as ACTUALLY fitted — local dry-runs
overwrite them with local env defaults). The lnN→abcd DEFAULT FLIP for
leppt_mt40 is still an open decision, as is the B/D boundary question
(m_T<30+buffer as implemented vs tiled m_T<40 — user to revisit; tiling
moves A0 by only +2.0–2.6%).** **`lumi` lnN 1.03
(±3% on L = 46.5 nb⁻¹) on EVERY MC template in all 50 channels, never on the
data-driven `qcd`** — degenerate direction absorbed by the constraint, so each
r gains a fully-correlated ~3% (propagates to σ via r×σ_gen, cancels in
charge-asym/F/B through the r-covariance automatically). The κ values the cards
were built with are recorded in the `datacards/qcd_lnn_kappas.txt` sidecar,
read back by the driver at extraction time (kQcd=0 ⇒ legacy free mode).
**LHE SHAPE SYSTEMATICS (2026-09-07, `LHE_SYST=auto|off|<list>`, default auto):**
the inputs' `<process>_{nPDF,qcdScale,alphaS}Up/Down` templates (see
"Structured inputs") become three `shape` rows — entries `1` on
`signal/z/ztau/wtau` (W channels) and `zsig/w/wtau/ztau` (Z peaks), `-` on
the data-driven `qcd` and on EVERY CR column — plus the fifth `shapes` token
`<dir>/<proc>_$SYSTEMATIC` on those MC processes (explicit paths: `$PROCESS`
would not match `zsig`/`w_y*`) and a `lhe group = nPDF qcdScale alphaS` line
(`--freezeNuisanceGroups lhe` = the stat-only comparison). The generator reads
the FOUR inputs' `_systs.txt` sidecars (copied next to the input copies by the
driver, uploaded by `sync_lxplus.sh`) and uses their common list; `off` gives
cards identical to the pre-09-07 ones up to one header comment; a subset not
present in all four is a hard error (exit 2). One row name = one θ for the
whole card, so the variation is fully correlated across channels, flavours,
charges and processes (what makes the charge asymmetry nearly nPDF-free). The
`qcd_lnn_kappas.txt` sidecar gains `lheSysts a,b,c|none`; the driver passes it
as the extractor's 9th argument → `<name>_theta` rows in `comb_summary.csv`
(value = pull, error = post-fit constraint; < 1 means the data constrain that
shape — for `qcdScale`, which re-shapes the lepton pT, that is expected and
means the data are tuning μR/μF), Asimov closure requires every θ = 0, and a
final sweep prints ANY floating parameter of `fit_s` not otherwise reported.
Combine builds one morphing pdf per flagged column: nominal at θ=0, Up at
θ=+1, Down at θ=−1, quadratic interpolation inside, linear beyond, the
normalization change carried as an asymmetric log-normal (`shape`, not
`shapeN`). Consequence for THIS fit: the raw variation's normalization part is
degenerate with r_i bin by bin, so θ_nPDF/θ_alphaS should come back ≈ 0 ± 1
and the r_i errors grow by roughly the normalization shifts, correlated across
bins — and since σ_meas = r_i × σ_gen(nominal), that counts the theory
normalization uncertainty as σ uncertainty (the deferred gen-level twins would
remove it). `run_pO_impacts.sh` discovers the three nuisances automatically
(+6 fits per variant). Local checks 2026-09-07: `LHE_SYST=off` cards ≡
baseline (+1 comment line) in lnN/abcd/met, all 8 systematic rows aligned
(248 / 324 columns), dry-run copies the sidecars, extractor + `postfit_incl.C`
compile. **Not yet fitted — needs the next lxplus run (`--asimov`); then
`sync_lxplus.sh download` also pulls the two `fitDiagnostics_simfit_<B>.root`
so `postfit_incl.C` can use `shapes_fit_s`.**
**LEPTON-SF SHAPE SYSTEMATIC (2026-09-14):** the muon inputs' sidecars
additionally list **`muSF`** (`<process>_muSFUp/Down` from `skim/muon_sf.h`:
the ID, ISO and trigger SF shifts added in quadrature per bin — ONE combined
nuisance by user decision, see Stage 1) — it enters the cards EXACTLY like
nPDF/qcdScale/alphaS: one `shape` row = one Gaussian-constrained nuisance θ,
templates at θ = ±1, pull and constraint in `comb_summary.csv`
(`muSF_theta`), Asimov closure at θ = 0. The generator now takes the UNION
of the four sidecars and keeps per-INPUT process lists
(`LHEP_{W,Z}_{MU,ELE}`; `lhe_append <W|Z|CR> <mu|ele> …`): a systematic gets
`1` only on the columns of the inputs that list it and `-` everywhere else —
the `muSF` row has entries on the 96 muon W SR MC columns + the 4 muon Z
columns only (verified: 0 on electron, qcd and CR columns; the theory rows
are unchanged on both flavours), and a theory systematic missing from an
input now WARNs instead of being dropped (`LHE_SYST=a,b` requires each name
in SOME sidecar). Groups: `lhe` = the theory names (`LHE_THEORY`), `lepsf` =
every other listed name, today `muSF` (`--freezeNuisanceGroups lhe,lepsf` =
stat-only). The sidecar's `lheSysts` line lists every shape nuisance
(`nPDF,qcdScale,alphaS,muSF`) plus a `sfTrigCorr` line; the extractor,
`postfit_incl.C` and the impacts consume it unchanged. **Dormant machinery
for the per-source alternative** (should the analysis repo list `muID`,
`muIso`, `muTrig` separately again): `SF_TRIG_CORR=auto|perbin|coherent`
(default auto = the muon W sidecar's directive `#! muTrig corr`, `coherent`
when absent) — `perbin` appends 12 lines `nuisance edit rename *
mu_W[pm]_<B>_y<i> muTrig muTrig_y<i>` to the card; Combine's
`DatacardParser.systematicsShapeMap` maps each new name back to the SAME
`<proc>_muTrigUp/Down` histograms (no duplicated templates;
`NuisanceModifier.fullmatch` on the channel regex, so y1 ≠ y10), the Z-peak
piece stays under `muTrig`, and `ModelTools.doNuisancesGroups` validates the
`lepsf` members after the edits, so post-edit names are listed. Local
dry-runs 2026-09-14 on the regenerated inputs (`QCD_MODE=abcd` leppt_mt40 and
lnN met): 324 / 248 aligned columns in every row, `lepsf group = muSF`, 0
edit lines (the `perbin` path was exercised earlier the same day with the
three separate sources: 12 edits, 96+4 columns each). NOT yet fitted (next
lxplus run).
**STAT/SYST SPLIT OF THE POI ERRORS (2026-09-15; supersedes the 2026-09-14
companion-fit design after the user asked why a second fit was needed):** the
error Minuit quotes for a POI is the TOTAL (profiled) one — every nuisance
may move while the POI is displaced — and its statistical part is NOT a
separate number in `fit_s`, but it IS contained in the full post-fit
covariance matrix `fit_s->covarianceMatrix()` (46 × 46 in abcd mode: 25 POIs
+ 12 ABCD scales + 9 constrained nuisances). Splitting V into the
statistical block k (POIs + unconstrained rateParams `r_*`, `qcd_s*`,
`qcd_norm*` — data-driven quantities) and the nuisance block n, the
covariance of k with n held FIXED is the Schur complement
V_stat = V_kk − V_kn V_nn⁻¹ V_nk (= (H_kk)⁻¹ with H = V⁻¹), i.e. exactly what
a refit with `--freezeParameters allConstrainedNuisances` at the post-fit
values returns in the Gaussian (Hesse) approximation — the Combine
"breakdown" recipe without the refit. `extract_pO_simfit.C::ComputeStatCov`
does this per variant and writes the 19th CSV column `rErr_stat`, the
`simfit_<B>_stat` rows (POIs with their stat error, r_Z, W⁺/W⁻/W sums) in
`comb_summary.csv` and the matrices `h_cov_yield_stat` /
`h_cov_yield_FB_stat`; it prints the kept/conditioned parameter lists (a new
unconstrained parameter under another name would show up as "conditioned
out" instead of silently shrinking the stat errors) and the mean stat/total
ratio. Checked on the 09-14 night fit: Schur complement vs full inversion
agree to 1e-16, √V_ii reproduces every quoted error to 2e-4 (the errors ARE
the Hesse diagonal), and stat ≤ total holds by construction (a PSD term is
subtracted), so syst = √(total² − stat²) is always real. **Numbers (lab):
stat/total = 0.85–0.91 per r_i, 0.80 for r_Z; per-bin syst ≈ 0.031–0.047
in r (lumi 3% dominates); the 12 ABCD scales lose < 4% of their error when
the nuisances are fixed (they are statistics, and stay in k).** Consumer:
`xsec_fiducial_comb` (inner/outer bars, `sigma_meas_stat/syst_nb`). The
**SOURCE FLIPPED 2026-09-15b (user decision): the `--statonly` COMPANION FIT is
now the PRIMARY stat source and runs BY DEFAULT**, with syst = √(total² − stat²)
downstream; the conditioned covariance stays as the FALLBACK when the companion
is absent and is always computed as the CROSS-CHECK (the two agreeing = the
likelihood is parabolic in the POI directions). Everything reads through one
pair of lambdas (`statErr`/`statCov`), so `rErr_stat`, `h_cov_yield[_FB]_stat`,
`h_cov_poi[_FB]_stat`, the `simfit_<B>_stat` rows and the inclusive sums all
switch together, and one printed line says which source was used. **NB the
conditioning guarantees stat ≤ total (a PSD term is subtracted) but a separate
refit does NOT** — POIs with stat > total are counted and WARNed ("syst floored
at 0") instead of being silently hidden by a downstream `max(0, ·)`. Verified
2026-09-15b: with no companion present the extraction is BIT-IDENTICAL to the
pre-change one (both CSVs), and standing in a copy of the nominal result as the
companion drives stat = total / syst = 0 all the way through to
`xsec_contour.C` — i.e. the branch really is what feeds every consumer.
Between 2026-09-15 and 09-15b the roles were reversed: the companion fit was
OPT-IN (`--statonly`, default off) and a
CROSS-CHECK only: when `fitDiagnostics_simfit_<B>_statonly.root` exists the
extractor writes `simfit_<B>_statonlyfit` rows and prints max |r_fit − r|
(WARN > 0.01) and max |rErr_fit/rErr_stat − 1| (WARN > 2%: non-Gaussian
likelihood or a mis-classified parameter). **`--extract-only` (same day)**
re-runs only the extraction on an existing `fits/` tree (only `root`; prefit
integrals from the work-dir input copies, else the analysis plots dir), and
the extractor now compares those integrals with the fit's own
`shapes_prefit/<F>_<C>_<B>_y<i>/signal` sums (agree to 1e-8 for the fitted
inputs; WARN > 1e-5 = regenerated inputs). `sync_lxplus.sh download` pulls
the nominal, `_statonly` and `_asimov` fitDiagnostics so a local
re-extraction reproduces the full summary (the 09-14 night download predates
the `_asimov` addition — its local re-extraction lacks the closure rows until
the next download). syst here = everything profiled in the fit INCLUDING
lumi 3%; a lumi-separate presentation would condition on `lumi` only (a
one-line variant of the same conditioning, not implemented).
`w`/`wtau` under the Z peaks FROZEN at absolute MC (measured 0.03–0.06 events
vs 372/252 signal — decision 2026-08-04; the lumi lnN does ride on them).
Implemented WITHOUT combineCards.py:
`my_script/make_pO_simfit_cards.sh` writes one 50-channel card + the
multiSignalModel map file per variant (`t2w_maps_simfit_<B>.txt` = THE model
definition; regex maps like `(mu|ele)_Wp_lab_y3/(signal|wtau)$:r_Wp_y3[1,0,10]`
— the trailing `/...` keeps y1 from matching y10, the `$` keeps `z` from
matching `ztau`); `text2workspace.py -P ...:multiSignalModel` +
`combine -M FitDiagnostics --skipBOnlyFit` (no `--rMin/--rMax` — there is no
POI named `r`). `--asimov` adds a prefit-Asimov (`-t -1`) closure fit per
variant: every POI must return 1, checked/PASS-FAIL'd by the extraction.
`my_script/extract_pO_simfit.C` reads the one `fit_s` per variant → all POIs +
the r-correlation matrix → `simfit/summary/comb_W_yields.csv` (first 9 columns
= legacy layout; qcd columns carry the MULTIPLIER κ^θ in lnN mode, same
semantics as the old rateParam values; `,lumi,lumiErr` appended 2026-08-17),
`comb_summary.csv` (POIs, status/covQual, covariance-propagated Wp/Wm/W
inclusive sums, Asimov rows, and the **nuisance pulls `<name>_theta`** — the
ABCD-consistency check: |θ̂|≳1 ⇒ κ too small or template biased; Asimov closure
also requires every θ = 0), and `comb_fitted_yields.root`:
`h_yield_W{p,m}_y{0..11}(_FB)` with yields = r×(S_mu+S_ele) (+`h_mt_*` aliases)
PLUS the 24×24 yield-covariance TH2Ds **`h_cov_yield` (lab) / `h_cov_yield_FB`
(fb)**, fixed order [Wp_y0..11, Wm_y0..11] (axis labels set). (These counts
are the fit's RECORD — since 2026-09-22 no observable is built on them;
`analysis/fiducial_yields.C` first turns the r's into r × σ_gen.) Postfit plots per
grand-fit channel (draw_postfit_pO.C grew optional trailing args poi/dy/qcd
name + ndfParams; defaults = legacy behavior). Statistical motivation: the
legacy per-bin scheme re-fits the SAME Z data 48× per flavour ignoring the
induced correlations, and never combines μ/e; the grand fit fixes both and
supplies the covariance the downstream error propagation now uses. μ/e-shared
r's ⇒ per-flavour observables are 100% correlated — the comb result is THE
result; independent per-flavour results come from `flavfit` (above; until
2026-09-22 from the removed legacy per-bin fits).

Mode `all` = simfit (channel `both`) + flavfit. The `simfit` mode ignores the
channel argument (always μ+e); its out-tree is `pO_fit_out<suffix>/simfit/`
{datacards, fits/simfit_{lab,fb}, contour, postfit, summary}.
(`--draw-only` redraws postfit plots from an existing run — needs only `root` +
the `fits/` tree, e.g. after cosmetic `draw_postfit_pO.C` changes. **`--disc
met|leppt_mt40`** (2026-07-30; default leppt_mt40 since 2026-08-16) selects
the W discriminant: reads `combine_input_W[_leppt_mt40].root`, writes to
`pO_fit_out[_leppt_mt40]/`, sets the postfit x-title; datacards/fit model
identical across discs (QCD lnN + lumi lnN everywhere since 2026-08-17 — the
lnN matters MOST for the pT variants, which lack the low-MET in-fit QCD
anchor; in the met fit the same nuisance is data-constrained, a cross-check of
κ). **Downstream continuation of the same
tag (2026-08-03):** `analysis/run_observables.sh <disc>` runs the whole
Module-5 observables chain per variant into disc-tagged files/folders — see
Stage 3 and the `observables.C`/`xsec_fiducial.C` bullets.)
It (1) locates this repo's **structured** Combine inputs, (2) generates the
datacards + multiSignalModel maps, (3) runs `text2workspace`+`combine -M
FitDiagnostics` per binning variant into a clean output tree, (4) extracts
POIs, fitted yields and covariances, (5) draws postfit plots in the same
cosmetics as `plotting/mtandmet.C`.

**Full step-by-step runbook (both repos):** this repo's `README.md` is the
single end-to-end procedure (skim → ngen → ABCD QCD → structured inputs → fit →
charge_asym/FBratio → observables). The fork's `test/README_pO_fits.md` is now
just the fit-stage quick reference pointing back to it.

Scripts on `zheng/po-analysis`:
- `test/run_pO_fits.sh` — master driver (bash-3.2 safe; `--dry-run` builds
  datacards without `cmsenv`). Modes `simfit | flavfit | all`; one work dir per
  fit: `test/pO_fit_out<suffix>/simfit/` (grand) and `simfit_{mu,ele}/`
  (flavfit), all set up by `run_simfit_set "<flavours>"` (2026-09-22).
- `test/my_script/make_pO_simfit_cards.sh` — **2026-09-22: `SIMFIT_FLAVS`
  (default `"mu ele"` = the grand card, byte-identical to before; `mu`/`ele` =
  the per-flavour card of `flavfit`: that flavour's channels, its two QCD rows,
  only the shape rows its OWN sidecars list, map regexes `<flav>_…`, sidecar
  line `flavours <list>`; the other flavour's inputs are not read).** —
  **2026-09-07: LHE shape rows
  from the inputs' `_systs.txt` sidecars (`LHE_SYST`, `lhe_append` helper
  extends the ALIGNMENT RULE to the shape rows; `lheSysts` sidecar line)** —
  simfit 50-channel datacards +
  multiSignalModel maps (the 25-POI model definition), lab + fb. Since
  2026-08-17 also the systematics: 4 QCD lnN rows (κ μ 1.15 / e 1.20, env
  `QCD_LNN_MU/QCD_LNN_ELE/QCD_MODE`) + the `lumi` lnN 1.03 row (`LUMI_LNN`),
  plus the `qcd_lnn_kappas.txt` sidecar. **2026-08-23: `QCD_MODE=abcd`**
  (62-channel card: + the 12 CR channels, the `qcd_s*` scales + formula
  rateParams, the `w_y*` map lines; env `QCD_ABCD_LNN_MU/ELE`, `QCD_WCR`;
  hard-errors unless disc = leppt_mt40; sidecar gains `qcdMode`). ALIGNMENT
  RULE unchanged: every new (channel, process) column appends one entry to
  MB/MP/MI/MR AND all five positional systematics rows in the same block.
- `test/my_script/extract_pO_simfit.C` — **2026-09-22: trailing args
  `flavours` ("mu,ele" default | "mu" | "ele") and `outTag` ("comb" default |
  "simfit_<flav>")** for the per-flavour fits: only the fit's flavours' W
  inputs are opened (the other may be "none"), S sums over them, the prefit-S
  check compares with that many fitted channels, the per-flavour QCD
  reporting covers them alone, the absent flavour's CSV qcd columns are 0,0,
  and the three output files carry the tag; the 19-column CSV layout is
  unchanged. Regression 2026-09-22: default args reproduce the old extractor's
  CSVs and all 104 histograms exactly. — **2026-09-15: also writes the 25×25
  `h_cov_poi[_FB][_stat]`** — cov(p_a, p_b) in PARAMETER space, order
  [r_Wp_y0..11, r_Wm_y0..11, **r_Z**], axis-labelled with the POI names. The
  existing `h_cov_yield` is in yield space and covers the 24 W POIs only, so
  the σ_W↔σ_Z cross term that `plotting/xsec_contour.C` needs had no home;
  this one is written unconditionally (no prefit-template dependence) and is
  the cleaner input for any r propagation. Cheap: `--extract-only` on an
  already-downloaded tree produces it, no refit. — **2026-09-07: 9th arg `lheSysts`
  (from the sidecar via the driver) → LHE nuisance pull + constraint rows,
  Asimov closure on them, safety-net sweep of any unreported floating
  parameter** — simfit extraction: POIs + covariance →
  `comb_*` CSVs + `comb_fitted_yields.root` (+`h_cov_yield[_FB]`), Asimov
  closure check (incl. nuisance θ=0 since 2026-08-17); trailing κ args (from
  the sidecar via the driver) switch qcd readout between κ^θ and legacy
  rateParams; nuisance pulls dumped to `comb_summary.csv` + console.
  **2026-08-23: 8th arg `qcdMode`** — "abcd" makes `AbcdScale` evaluate the
  formula multiplier κ^θ·sB·sC/sD from `floatParsFinal` with full-correlation
  errors (the RooFormulaVar itself is not there), stamps the `qcd_model` CSV
  column, writes the scales + `qcd_abcd_mult_*` summary rows, and extends the
  Asimov closure to the 12 scales (= 1).
- `test/run_pO_impacts.sh` (2026-08-17) — nuisance impacts + covariance plots
  for the simfit, run AFTER `run_pO_fits.sh`: `./run_pO_impacts.sh [--disc …]
  [--variant lab|fb|both] [--fit comb|mu|ele|all] [--pois "…"|all]
  [--impacts-only|--cov-only] [--dry-run]` (`--fit`, 2026-09-22: the
  per-flavour fits' work dirs `simfit_{mu,ele}/`, outputs inside them;
  **`--fit all` = grand + μ + e in one call**, or a list like `mu,ele`,
  validated before any work; a fit whose work dir is absent is skipped with a
  note, so `all` also works where flavfit never ran — exit 1 only if NO fit
  was found. Dry-run + a real `--cov-only` on the flavfit fixture tree
  verified 2026-09-22.) Impacts via `combineTool.py -M Impacts` + `plotImpacts.py
  --POI` — both SHIPPED with Combine v10 in `<fork>/scripts/` (no
  CombineHarvester) — one PDF per POI (default all 25) →
  `simfit/impacts/impacts_simfit_<B>[_<poi>].{json,pdf}` (intermediates in
  `impacts/wd_<B>/`, not synced). Covariance via `my_script/plot_pO_cov.C`
  (root-only, also runs locally): full parameter correlation from `fit_s`
  (`corr_params_simfit_<B>`; skipped+WARN when fits/ absent) + the 24×24
  fitted-yield correlation from `h_cov_yield[_FB]` (`corr_yield_<B>`) →
  `simfit/cov/`. **EXTRACTION LOG (2026-09-15c):** the `[extract-simfit]` console output — which stat source was used, the statonly-vs-conditioned cross-check, the Asimov closure, every WARN — is printed, NOT stored in the CSVs, and the fit runs on lxplus, so it used to exist only in that terminal. It is now teed to `simfit/summary/extract_simfit.log` (and the legacy path to `<chan>/summary/extract_<chan>.log`); `summary/` is downloaded wholesale, so it travels with the results. The download also pulls the per-variant `fits/simfit_<B>/*.log` (t2w / fit / fit_statonly / fit_asimov) — small text, and the only place a badly converged fit explains itself. `sync_lxplus.sh download` pulls both dirs (excluding `wd_*`);
  `upload-scripts` ships the two new files.
- **`--contour` (2026-09-15; ON BY DEFAULT since 2026-09-15b) — the PROFILED
  (σ_W, σ_Z) region.** A simfit run now does THREE fits per binning variant:
  the nominal one, the `--statonly` companion, and this scan (`--no-statonly`
  / `--no-contour` opt out; `CONTOUR_POINTS` default 2500). **It is purely
  ADDITIVE** — the SAME datacard (only an extra map file), its own workspace,
  its own output dir `simfit/contour/contour_<B>/`; nothing in `simfit/fits/`
  is touched and the extraction never reads it. It runs AFTER the nominal fits
  in the same invocation. A missing gen sidecar is a WARN that skips this pass
  ONLY — it must never take the nominal fit down with it. It reparametrizes the SAME 62-channel datacard through a
  SECOND map file, `t2w_maps_simfit_<B>_sigma.txt`, written by
  `make_pO_simfit_cards.sh::write_sigma_maps` when `SIGMA_POI=1` (the driver
  sets it, and `GEN_XSEC_FID`, defaulting to `$PO_PLOTS/../../skim/output/
  gen_xsec_fid.txt`):

      POIs   sigmaW, a_<C>_y<i>   (23 free; ONE bin anchors with a == 1)
      funcs  Wnorm = G_anchor + Sum_j a_j G_j
             r_i   = sigmaW * a_i / Wnorm      (anchor: r = sigmaW / Wnorm)

  so Σᵢ rᵢ Gᵢ ≡ sigmaW identically (Gᵢ = the gen fiducial σ of
  `gen_xsec_fid.txt`). Still 25 POIs, still 24 dof in the W sector. The
  **datacard is untouched** — POIs live only in the map — so nothing about the
  nominal fit changes and one card serves both parametrizations. The profile
  likelihood is invariant under a bijective reparametrization, so the best fit
  is identical and what this buys is the exact profiled 2D region (of which
  `plotting/xsec_contour.C`'s ellipse is the Gaussian approximation).
  Mechanics that matter: `MultiSignalModel` supports factory expressions via
  `map=<patterns>:<name>=expr;;<name>("…",args)` (`;` → `:`; declarations are
  `doVar`'d before factories, and factories run in MAP-FILE ORDER, so `Wnorm`
  must come first); each POI gets **ONE line carrying BOTH the SR and the CRB
  pattern** (comma-separated) because re-listing a factory-defined name on a
  second plain `map=` line would push it into the POI set as a bare variable
  and collide with the RooFormulaVar. Declarations hang off never-matching
  patterns (`__poi_decl__/__poi_decl__`). Two combine calls per variant:
  `--algo singles` locates the minimum and sets the grid window (±4σ), then
  `--algo grid` scans (sigmaW, r_Z) with everything else profiled — σ_Z is a
  one-to-one rescaling of r_Z, so no second POI is needed. Output
  `simfit/contour/contour_<B>/` (`sync_lxplus.sh download` pulls it, excluding
  the ~100 MB `workspace_sigma.root`). **The standard extractor is NOT run on
  it** — under this parametrization `r_<C>_y<i>` are RooFormulaVars, absent
  from `floatParsFinal`; every per-bin number keeps coming from the nominal
  fit, which by invariance has the same minimum. **Chosen over the user's
  "redefine one bin as σ_total − Σ(others)"**: the subtraction drives that
  bin's yield negative once σ_W is scanned ~1.4σ down (the largest bin carries
  only ~5% of the total), well inside a 95% contour, whereas the ratio form
  keeps every rᵢ > 0 by construction. **Validated locally 2026-09-15 without
  cmsenv** (three independent checks, all green): the RooFit algebra closes to
  6e-16 over 200 random points; the generated map files run through the fork's
  OWN `MultiSignalModel` parser with a real `RooWorkspace` behind
  `modelBuilder` and give 25 POIs, all `RooRealVar`, with **all 324 (bin,
  process) columns of the datacard assigned the IDENTICAL scale name** as the
  nominal map (the reparametrization changes nothing about what scales what)
  and rᵢ = 1.00000000 at the init point (Asimov closure unmoved); and a
  synthetic exactly-Gaussian scan file exercises `xsec_contour.C`'s reader
  end-to-end (it caught a real bug: `quantileExpected` is −1 on EVERY
  MultiDimFit row, so an early attempt to use it to drop the duplicate
  best-fit row would have discarded the entire scan). **Not yet run on
  lxplus.**
- `test/run_pO_idiso_sf.sh` + `my_script/make_pO_idiso_sf_cards.sh` +
  `my_script/extract_pO_idiso_sf.C` (2026-09-24/25) — **the electron ID+iso
  SF cross-check, a SEPARATE stream** (see the `correction/idiso_sf_*`
  bullet): reads `correction/rootfile/idiso_sf_ele/`, writes
  `test/pO_idiso_sf_out/<tag>/` — never `pO_fit_out*`, and `run_pO_fits.sh`
  knows nothing about it. 9 W fits (lepton pT at m_T > 40, W pass/fail SR +
  the in-fit ABCD CRs `ele_<R>_CR{B,C,D}` per channel with `qcd_s{B,C,D}_ele_<R>`
  and the formula `qcd_abcd_ele_<R>` = (@0*@1/@2) — the nominal QCD_MODE=abcd
  syntax — + the Z_PP DY anchor; POIs `r_<C>_<cat>_k<k>` (SR and CRB W) +
  `r_Z`; residual lnN `qcd_rate_ele_<C>_<cat>`, κ `IDISO_QCD_ABCD_LNN` default
  1.15) + ONE Z tag-and-probe fit (tag `ztnp`: Z_PP/Z_PF, `r_ZPP`/`r_ZPF`);
  extraction → `summary/idiso_sf_<tag>.csv` (ε_data, ε_MC, SF per bin +
  charge; the Z row is `ztnp`) + `extract_idiso_sf_<tag>.log` (+ the SR ABCD
  multiplier κ^θ·sB·sC/sD per channel); `--asimov` closure (every r and scale
  = 1, θ = 0). Synced by `sync_lxplus.sh upload-idiso|download-idiso` (never
  by plain upload/download). Validated locally 2026-09-25 (cards aligned, all
  1399 (bin, process) columns mapped as intended by the fork's own
  MultiSignalModel, stand-in Asimov closure exact); run on lxplus 2026-09-26
  (all 10 fits clean, Asimov PASS — results in the `correction/idiso_sf_*`
  bullet). The extractor also prints `[expected] … SF 1 +- X` from each Asimov
  fit: on data a fail r at its 0 boundary has a meaningless error.
- (REMOVED 2026-09-22 with the legacy per-bin pipeline, in git history:
  `test/my_script/make_pO_datacards.sh` — the 53 legacy cards per flavour;
  `test/my_script/extract_pO_yields.C` — their `<chan>_{W_yields.csv,
  summary.csv,fitted_yields.root}`; `test/my_script/make_yields_from_csv.C`.)
- `test/my_script/draw_postfit_pO.C` — generic postfit data/MC plotter routed
  through the synced `plotting_helper.C::SaveNicePlot1D_WithBkg` (α=0.65 fills,
  width-3 outlines, `CMS_lumi`) → identical look to `mtandmet.C`.
- `test/my_script/plotting_helper.C` — **synced** from this repo's
  `plotting/plotting_helper.C` (was an older solid-fill copy; `CMS_lumi.C` was
  already identical). Re-synced 2026-07-19 with the pull-pad helper;
  `draw_postfit_pO.C` sets `ps.pullPad=true`, so postfit plots gain the pull
  sub-pad on the next fit run (the local `pO_fit_out/` tree has no `fits/`
  dirs left, so `--draw-only` can't regenerate them without re-fitting).

**Structured inputs** (produced HERE, consumed by the fork):
- `plotting/mtandmet.C` → `plots[/Elec]/combine_input_W.root`: one **TDirectory
  per fit region** (`Wp_lab_y3/`, `Wm_fb_y7/`, `Wp_incl/`, `W_incl/`, …), each
  holding the 6 **absolute** templates `data_obs/signal/z/ztau/wtau/qcd` (MET
  discriminant; per-y ABCD `qcd` split from the inclusive template). `signal` =
  both W MC samples in that reco-charge region; `wtau` = Wptau+Wmtau.
- `plotting/dileptonpeak.C` → `plots[/Elec]/combine_input_Z.root`: a `Z_incl/`
  dir with `data_obs/signal/w/wtau/ztau` (mass peak, absolute).
- **LHE shape systematics (2026-09-07, both writers):** every MC process of
  every region additionally carries `<process>_<syst>Up` / `<process>_<syst>Down`
  for `syst` ∈ {`nPDF`, `qcdScale`, `alphaS`} (`pOLhe::kLheSystNames` in
  `skim/lhe_index.h`, = the twins `skim/lhe_updown.py` wrote into the MC skim
  files). W inputs: `signal/z/ztau/wtau` in all 51 dirs (48 SR + the 3
  inclusive; the CR dirs and the data-driven `qcd`/`qcd_abcd` carry none) =
  1224 extra TH1Ds per file; Z input: `signal/ztau/w/wtau` in `Z_incl/` = 24.
  Built by `attachSysts()` in `mtandmet.C` (twin block in `dileptonpeak.C`)
  from the per-sample `<hist>_<syst>Up/Down` with the SAME k_s scaling and
  W⁺+W⁻ sample summing as the nominal — the nominal templates are untouched
  (verified bit-identical with `compare_hists.C`, which recurses into the
  dirs). Present in the `met` and `leppt_mt40` files; the plain `leppt` file
  had none (no twins in the skim) — which is why it was retired 2026-09-21
  rather than merely left unfitted. A syst is written only when
  it is complete in EVERY SR region, and the **sidecar `<input>_systs.txt`**
  next to each file (`<syst> <processes…>`, `#` comments) records exactly
  what was written — the fork's card generator must read it to emit the
  `shape` rows, never assume. Combine convention: `shapes <proc> <ch> <file>
  <dir>/<proc> <dir>/<proc>_$SYSTEMATIC` (Combine appends `Up`/`Down`).
- **Muon-SF shape systematic (2026-09-14, muon inputs only):** the same
  writers also carry the ONE combined `<process>_muSFUp/Down` for every MC
  process of every SR + inclusive region (W) and of `Z_incl` (Z), from the
  skim's `<hist>_muSFUp/Down` twins (`skim/muon_sf.h`: ID, ISO and trigger
  shifts in quadrature per bin; NOT area-normalized — the normalization is
  the effect). The per-flavour list is built in `mtandmet.C`/`dileptonpeak.C`
  (`systNames` = `pOLhe::kLheSystNames` + `pOSF::kMuonSFSystNames` when
  `!isElec`), so the muon sidecars list 4 systematics and the electron ones 3
  — the card generator puts a flavour-specific row on that flavour's columns
  only. (A `#!` line in a sidecar is a directive for the card generator, not a
  systematic — `sidecar_systs()` skips every `#` line; the only one defined,
  `#! muTrig corr <coherent|perbin>` from `pOSF::kMuTrigCorr`, is written
  only when a separate `muTrig` nuisance is shipped, i.e. not today.)

**Legacy fit model — REMOVED 2026-09-22** (it ran as modes
`perbin|incl|combined`): per flavour, 48 separate two-channel cards (one W
region + `Z_incl` each) plus incl/WZ cards, with the TWO-PARAMETER model — POI
`r` on all W-related MC, `dy_norm` on all DY-related MC (shared with the Z
peak), a free `qcd_norm` per region — and no systematics; its statistical flaw
was re-fitting the same Z data 48× per flavour with the induced correlations
ignored. The per-flavour fit is `flavfit` now (the simfit model and
treatment, one flavour per likelihood).

All templates are **absolutely** normalized (`k_s = A·σ·L/N_gen`) — NO area
normalization (the old `make_combine_input*.C` `data_int/mc_sum` rescale is gone).

**Superseded** (kept for reference, not in the new pipeline): `test/run_fit.sh`,
`test/my_script/make_combine_input{,_Z}.C` (area-normalized + Rayleigh `pdfbkg`),
`testdatacard_{inclusive,Zmumu,Zee}.txt`, `draw_postfit_{inclusive,Zmumu,Zee}.C`.

Fit method: `combine -M FitDiagnostics --saveShapes --saveWithUncertainties`.

## Current state of the work

Active (in flight, recent commits):
- **ELECTRON ID+ISO SF CROSS-CHECK (2026-09-24..26, separate stream — the
  `correction/idiso_sf_*` bullet), FITTED WITH COMBINE 2026-09-26:** the
  W-sample route cannot measure the fail category even with the nominal in-fit
  ABCD (W SF 0.41–1.24 across windows — both ends being the r_fail = 10 and
  r_fail = 0 limits of the POI range, not measurements — 0.75 vs 1.24 across
  binning schemes of the same data, fail postfit χ²/ndf 46; expected
  ±0.11–0.14); the inclusive Z
  tag-and-probe, its own fit, gives **SF 0.9588 ± 0.0251 vs EGM `wp90iso` 0.956 ±
  0.009** for the same probes. Plots in `correction/plots/idiso_sf_ele/`.
  **Next (user's call):** whether the electron selection moves to WP90 ID+iso
  with the pp SF applied.
- **μ-ONLY / e-ONLY SIMULTANEOUS FITS (`flavfit`) + μ-vs-e OVERLAYS — BUILT
  2026-09-22, FITTED on lxplus the same evening, observables regenerated
  2026-09-23.** **RESULTS (leppt_mt40, abcd, all shape nuisances):** fit
  quality clean for both flavours and both binnings — status 0 / covQual 3,
  Asimov closure PASS, `--statonly` vs conditioned covariance ≤ 0.06%, all six
  profiled contours (comb/μ/e × lab/fb) pass the staleness gate (scan minimum
  vs covariance best fit ≤ 0.009 nb). **μ-only: every pull |θ| < 0.25**
  (qcdScale +0.12, nPDF −0.13, muSF −0.10, qcd_rate_mu ±0.2, QCD multipliers
  0.93/0.93), r_Z 1.073 ± 0.065. **e-only carries ALL the tension:**
  qcd_rate_ele_Wp −1.57 ± 0.76 / _Wm −2.78 ± 0.72 (fb −1.40 / −3.55),
  multipliers 0.76/0.63, nPDF +0.59, qcdScale +0.55, r_Z 1.190 ± 0.080 — so
  the grand fit's theory pulls (qcdScale +0.63, nPDF +0.41) come from the
  electron channel. **σ (nb, stat ± syst):** μ W⁺ 59.66 ± 1.28 ± 1.82 / W⁻
  48.62 ± 1.16 ± 1.49 / W 108.28 ± 1.73 ± 3.30; e W⁺ 61.39 ± 1.61 ± 2.20 / W⁻
  47.69 ± 1.39 ± 1.68 / W 109.08 ± 2.13 ± 3.59; comb unchanged (108.23 ± 1.34 ±
  3.30 — it sits 0.05 nb below both flavours, not an error: per charge it lies
  between them, W⁺ nearer μ, W⁻ midway). **e/μ (stat): W⁺ 1.029 ± 0.035, W⁻
  0.981 ± 0.037, W 1.007 ± 0.025** — lepton universality at 2.5% with no
  electron SF applied; **σ_Z μ 9.41 ± 0.57 vs e 10.43 ± 0.70 nb, e/μ 1.108 ±
  0.088 (+1.2σ)**, the largest μ/e gap (the Z→ee peak's higher data/MC).
  Shape tests μ vs e (stat χ²): A_ch 11.6/12, R_FB 6.0/6 (W⁺ 10.4/6, W⁻
  2.9/6), dσ/dη 8.0/12, 11.3/12, 7.6/12 — every p ≥ 0.11. **Inclusive postfit
  stacks:** μ-only 5003 vs 4998 data, χ²/ndf 3.07, pulls scattered, NO
  threshold trend left; e-only χ²/ndf 6.36 with a deficit at 26–34 GeV and an
  excess over 45–95 GeV — the e pT spectrum is harder than the model, which
  is what the e QCD pulls absorb (e scale / missing e SFs / soft flat-T QCD
  template, the known suspects). **Impacts (`--fit all`, 2026-09-22):** lumi
  dominates every r (31–46% of σ_r); in the e-only fit the QCD rates reach
  44% (r_Wp_y11) / 37% in the edge bins; in the μ-only fit every other
  nuisance is ≤ 6%. The earlier build notes follow. User request:
  since the electron SFs are not applied, check μ and e separately with a
  dedicated per-flavour fit (still a simultaneous fit, not μ+e combined) that
  gets the SAME treatment as the grand fit (stat/syst split, contour, closure),
  and overlay μ vs e on every downstream plot (σ, A_ch, R_FB, contours). Fork:
  mode `flavfit` (see "Downstream fit"); analysis repo: `observables_flav`,
  `xsec_fiducial_flav`, `xsec_contour_WZ_flav` + `_fit`, `postfit_incl_fit`,
  `analysis/fiducial_yields.C`, `plotting/fit_variants.h`, the flavfit block of
  `run_observables.sh` → `plots/flavfit/`. **The legacy per-flavour per-bin
  pipeline was removed the same day** (user: "we can remove the old mu and ele
  option"): fork modes `perbin|incl|combined` + 3 scripts, the downstream
  `observables(isElec)`, `observables_overlay`, legacy `xsec_fiducial`/`_diff`,
  `SaveNiceGraph_ErrorBand_TwoData`. Their last OUTPUTS are orphaned and were
  left in place for the user to decide: `plotting/plots/{charge_asym,FBratio}/`,
  `plots/Elec/{charge_asym,FBratio}/`, `plots/merged/`, `plots/xsec/`,
  `skim/rootfile/{charge_asym,FBratio}_fit_{mu,ele}*.root`, and the fork's
  `pO_fit_out*/{mu,ele}/` trees (partly git-tracked in the fork).
  **Physics finding on the way — count-based R_FB is acceptance-distorted:**
  the count-based R_FB of the two flavours differs by up to 60% on identical
  physics (`fiducial_yields.C` bullet), and the same effect moved the PRIMARY
  comb R_FB by up to ±0.18 (2.5–3 stat σ) in the |η_CM| 1.2 and 1.9 bins.
  **DECIDED the same day (user): "all the observables should be using r × gen
  xsec, we should never use raw counts … since we are not going to do a
  dedicated eff or acc correction"** — so every A_ch / R_FB, comb included, is
  now built from `fidyields_<fit>_<disc>.root` (the σ's always were); the comb
  plots in `plots/comb/{charge_asym,FBratio}/leppt_mt40/` were regenerated
  from the 09-21 fit (R_FB_sum 1.026 → 0.848 at |η_CM| 1.88, 0.908 → 1.087 at
  1.20; A_ch ≤ 0.008; σ identical). Pre-switch comb outputs backed up in the
  session scratchpad only; the count-based `skim/rootfile/{charge_asym,FBratio}_fit_comb_*.root`
  are orphaned (nothing reads them). The met variant was NOT regenerated (its
  out-tree is the stale Aug-6 fit). (The fits were then run — results at
  the top of this bullet.) **PDF note (2026-09-23, checked because the wide
  1150×800 contour canvas looked broken):** ROOT writes every canvas that is
  wider than tall — incl. every 800×800 one, since batch mode shrinks it to
  796×772 — as a landscape page with `/Rotate 90`. Preview, Quick Look and
  pdflatex honour that flag and show it upright and complete; `sips` ignores
  it and renders a rotated, cropped fragment. So check PDFs with
  `qlmanage -t`, never `sips`. The one real cosmetic difference: dashed
  contour lines render dashed in the PDF but look solid in the PNG (ROOT's
  PNG backend restarts the dash pattern on every short segment of a dense
  polyline).
- **EMBEDDED MC (Aug-20, POWHEG + Angantyr) — SWITCHED AND FULLY RE-RUN
  2026-09-21, up to but NOT including the fit.** Repoint (`skim_common.h`
  `inputBase` + the `Aug_20_MC_` prefix in `count_ngen.C::LabelFromFname`) →
  N_gen → `gen_xsec.C` → `run_trig_eff_mb.sh both` → full 28-job re-skim
  (failures=0) → `run_lhe_updown.sh all` → `run_qcd_abcd.sh` → `mtandmet` ×2 +
  `dileptonpeak` ×2 → `run_syst_shapes.sh all`. **Sidecars unchanged (4
  systematics μ incl. `muSF`, 3 e) ⇒ the fork needs NO change**; κ moved only
  μ 1.09 → **1.10** (e stays 1.15 under the FF-shift recipe; the window⊕tilt
  variant gives 1.12), ABCD A40 μ 132.4 → 129.5, e 627.5 → 627.0. Old skims
  backed up in `skim/rootfile_pre_embed/`, comparison log
  `skim/logs/compare_embed.log`, chain log `skim/logs/rerun_embed_chain.log`.
  **REGRESSION PASSES: all four DATA files bit-identical** (`compare_reskim.sh`
  `DIFF 0`), and the Combine-input integrals reproduce the documented counts
  once the dropped overflow is added back (μ 5050 − 52 = 4998, e 4864 − 98 =
  4766). **WHAT THE EMBEDDING DID — measured as ε = Σw_sel/N_gen, which is
  N_gen-invariant, so it isolates reconstruction from the N_gen shift:**
  | | isolation ε | everything else | net |
  |---|---|---|---|
  | W→μν | **−2.1%** | +0.1% | **−2.2%** |
  | W→eν | **−2.8%** | **+3.5%** | **+0.6%** |
  The **isolation loss is the embedding working as intended** (Angantyr UE
  lands in the iso cone; both flavours lose 2–3% per lepton). Z→μμ confirms the
  size independently: ε fell **−4.6%** = (1 − 0.023)², the same per-muon loss
  counted twice for two legs, and `Z_incl/signal` 355.4 → **339.2** against
  unchanged data 364, so **Z→μμ data/MC worsened 1.024 → 1.072**. The
  **+3.5% electron gain in the pre-isolation chain (reco/ID/trigger/gap) is NOT
  an embedding effect** — embedding can only reduce efficiency — and it is NOT
  an energy-scale change (Z→ee MC mass peak 90.717 → 90.734, RMS identical).
  Measured directly as the `h_iso_met_*` denominator / N_gen: μ +0.06/+0.20%,
  **e +3.61/+3.42%**. **OPEN: ask the producer what changed in the electron
  reco/ID between July-29 and Aug-20 besides the embedding** — this moves the
  relative μ/e acceptance by ~3.5%, which is exactly what the simfit's
  μ/e-SHARED `r` takes from MC. Ruled out as the cause: the muon SFs
  (⟨SF⟩ 0.98792 vs 0.98794 documented). **The W normalization gap did NOT
  close**: inclusive data/MC on `leppt_mt40` is μ **1.169** / e **1.050**,
  i.e. still the r ≈ 1.19 of the 09-14 fit, and since the muon MC shrank
  expect r_μ to come out slightly HIGHER. So embedding alone does not explain
  the overshoot; the MET *shape* question is for the fit to answer.
  **FITTED AND FULLY PROCESSED THE SAME DAY (2026-09-21, abcd/leppt_mt40, all
  four shape nuisances):** status 0 / covQual 3 both variants, Asimov closure
  exact, and the `--statonly` companion agrees with the conditioned covariance
  to **1e-4** (max |rErr_fit/rErr_cond − 1| = 0.0001, max |r_fit − r| = 0.0000)
  — the likelihood is parabolic in the POI directions, so either stat source is
  defensible. **σ(W⁺) = 60.09 ± 1.00 (stat) ± 1.85 (syst) / σ(W⁻) = 48.14 ±
  0.89 ± 1.49 / σ(W) = 108.23 ± 1.34 ± 3.30 nb** (eff. r 1.207/1.215/1.210);
  **σ(Z) = 9.817 ± 0.387 ± 0.294 nb** (r_Z 1.1195 ± 0.0554), ρ(σ_W, σ_Z) =
  +0.52; fb σ_W = 93.22 ± 3.11 nb (r_eff 1.216). Both profiled contours drawn,
  scan minimum reproducing the covariance best fit to **7e-5** (lab), and the
  ellipse misstating the 1σ radius by −3.8%…+3.2%. **vs the July-29 fit: σ_W
  106.47 → 108.23 (+1.7%, r 1.191 → 1.210), σ_Z 9.530 → 9.817 (+3.0%, r_Z
  1.087 → 1.120) — exactly what the measured efficiency drops require (muon
  −2.2% in W, −4.6% in Z), i.e. the skim-level efficiency numbers and the fit
  close on each other.** Pulls (lab): nPDF +0.41 ± 0.985, **qcdScale +0.63 ±
  0.945 (was +1.16 — the embedded MC needs LESS help from the scale
  nuisance)**, alphaS −0.08, muSF −0.27, lumi +0.02, qcd_rate_mu_Wp/Wm
  −0.26/−0.10, **qcd_rate_ele_Wm −2.81 ± 0.63 (was −2.67 — still THE
  outlier)**. **THE RESULT THAT MATTERS — the turn-on deficit, re-measured as
  data/postfit summed over the 12 lab y bins × 2 charges: μ over 24–36 GeV
  −3.4% → −1.3% (more than halved), e −7.6% → −7.6% (unchanged to the
  digit).** That is exactly the efficiency decomposition above playing out:
  embedding cut muon isolation efficiency by 2.1% with nothing offsetting it,
  so the MC stopped over-predicting at threshold, while in the electron
  channel the 2.8% isolation loss was cancelled by the unexplained +3.5%
  pre-isolation gain, leaving the net efficiency — and hence the threshold
  shape — where it was. **So the embedding helps precisely where it changed
  the efficiency and does nothing where it did not, which sharpens the
  electron question: whatever produced the +3.5% may be cancelling a real
  improvement in that channel.** **BUG FOUND AND FIXED while processing
  (2026-09-21):** `run_observables.sh` called `xsec_contour_WZ("$disc")` with
  no binning argument, so only **lab** was ever regenerated and
  `xsec_contour_WZ_fb.{png,pdf,csv}` silently kept whatever the last manual
  run left — caught because the fb files still showed the PREVIOUS fit's
  σ_W = 92.01 nb against the correct 93.22. fb sums a DIFFERENT lab window,
  so it is a genuinely different number, not a cosmetic twin. The driver now
  loops both binnings with a `require_file` check on each. **Still open /
  next:** the electron pre-isolation +3.5% (ask the producer what changed in
  the electron reco/ID between July-29 and Aug-20 — it is not the energy
  scale, the Z→ee MC peak moved 90.717 → 90.734 with identical RMS), the
  `met` variant (its out-tree is still the Aug-6 pre-filter fit — refit before
  quoting anything from it), and `QCD_ABCD_LNN_MU` 1.09 → 1.10.
- **INCLUSIVE (σ_W, σ_Z) PLANE — ellipse DONE 2026-09-15, profiled contour
  BUILT AND VALIDATED LOCALLY, not yet run on lxplus.** Motivation (user): get
  the covariance of the total inclusive W and Z cross sections so the two can
  be shown as a contour with nPDF predictions overlaid. **First question
  answered: the current inclusive σ_W is NOT an approximation of a fitted
  inclusive value — it IS one.** σ_W = Σᵢ rᵢ σ_gen,i is LINEAR in the POIs, so
  (a) the MLE is invariant under the reparametrization (the best fit cannot
  move, by construction) and (b) `sumWithCov`'s √(gᵀVg) is exactly the Hesse
  error a reparametrized fit would report. The only difference is profile
  non-parabolicity, bounded analytically: the systematic is 3.25 of 3.51 nb
  and **3.19 nb of that is the lumi log-normal**, whose exact profiled interval
  is σ̂κ^±1 = +3.00/−2.91% instead of the linearized ±2.96% ⇒
  **+3.51/−3.42 nb vs the quoted ±3.51, i.e. ~0.1 nb of asymmetry, ~2% of the
  error.** Everything else is comfortably parabolic (r̂ᵢ ≈ 1.0–1.4 with ~0.07
  errors = 14σ from the r ≥ 0 bound; the ABCD ratio's ŝ_D = 1.003 ± 0.025;
  all four shape nuisances have post-fit widths ≈ 1). **The honest caveat is
  not statistical**: the number inherits the model's known tension
  (`qcd_rate_ele_Wm_theta` = −2.66 ± 0.62, inclusive postfit χ²/ndf 3.3 μ /
  7.2 e) and THAT is not in the error bar. Delivered: `h_gen_sig_Z` +
  `gen_xsec_fid.txt` (`skim/gen_xsec.C`), the 25×25 `h_cov_poi` (fork
  extractor, obtained by `--extract-only` — no refit), and
  `plotting/xsec_contour.C` (ellipse + CSV, wired into `run_observables.sh`).
  σ_W = 106.47 ± 3.51 / σ_Z = 9.522 ± 0.472 nb, ρ = +0.521 (≈ all lumi);
  the σ_W path reproduces `xsec_comb.csv` to every digit via a completely
  different matrix. **Also built: `run_pO_fits.sh --contour`**, the
  reparametrized workspace that makes σ_W a POI for the exact profiled
  contour (see "Downstream fit"; validated three ways offline). **NEXT:**
  run it on lxplus with the pending refit; then decide what goes on the
  theory axis — EPPS21's own uncertainty ellipse is obtainable in-house from
  the 107 member weights at gen level (= the deferred `gen_xsec.C` gen-twin
  TODO), but nCTEQ15HQ / nNNPDF3.0 / TUJU21 would need external ABSOLUTE
  fiducial predictions (`RpO_FB_graphs.root` holds ratios only).
  **Side effect to be aware of:** the local dry-runs replaced the
  `simfit/combine_input_*.root` work-dir copies with the CURRENT (2026-09-15,
  post-η-SF) inputs, which are NOT the ones the 09-14 night fit used — the
  extractor's own guard reports it ("prefit signal integrals … differ by up
  to 8.23e-03"). Confined to the `signal_prefit` / `fitted_yield` diagnostic
  columns of `comb_W_yields.csv`; r, rErr and every covariance are untouched
  (σ_W/σ_Z verified identical to the last digit), and the pending refit
  erases it.
- **FIT NOW RUNS THREE PASSES BY DEFAULT + STAT SOURCE FLIPPED — 2026-09-15b
  (user decision, BEFORE the refit).** `run_pO_fits.sh … simfit` per binning
  variant: (1) the nominal fit, (2) the `--statonly` companion (constrained
  nuisances frozen at their post-fit values), (3) the `--contour` scan
  (σ_W as a POI). (2) and (3) are ADDITIVE — separate files, separate dirs,
  nothing of (1) overwritten; opt out with `--no-statonly` / `--no-contour`.
  **The quoted statistical error now comes from (2), not from the covariance
  conditioning**, with syst = √(total² − stat²); the conditioning becomes the
  fallback and the always-computed cross-check. `sync_lxplus.sh download`
  already pulled the `_statonly` fitDiagnostics and now also pulls
  `simfit/contour/` (excluding the ~100 MB `workspace_sigma.root`).
  **FITTED AND VALIDATED 2026-09-15c (the first three-pass run):** status 0 /
  covQual 3 both variants, Asimov closure PASS (max |POI − 1| = 0.0000).
  Stat source is the companion fit, and **it agrees with the conditioned
  covariance to 0.0001** (max |rErr_fit/rErr_conditioned − 1|, both variants)
  with max |r_fit − r| = 0.0000 — so the refit-free conditioning was right all
  along and either source is defensible. No stat > total anywhere; the prefit
  signal-integral WARN is gone (the uploaded inputs ARE the fitted ones).
  **σ(W⁺) = 59.09 ± 0.98 ± 1.84, σ(W⁻) = 47.66 ± 0.88 ± 1.49, σ(W) = 106.75 ±
  1.32 (stat) ± 3.29 (syst) nb** (eff. r 1.187/1.204/1.195); σ_Z = 9.530 ±
  0.376 ± 0.286 nb, ρ = +0.521. Pulls: qcdScale +1.06 ± 0.96, nPDF +0.60 ±
  0.96, muSF −0.76 ± 0.95, alphaS −0.03, lumi +0.01, qcd_rate_ele_Wm −2.67 ±
  0.62 (still THE outlier). **THE PROFILED CONTOUR RAN AND MATCHES:** the
  reparametrized scan's minimum reproduces the covariance best fit to
  d(σ_W) = +0.0036 nb (3e-5 relative) and d(σ_Z) = −0.0007 nb — the invariance
  argument confirmed empirically, and the staleness gate passes. **Measured
  non-Gaussianity: the ellipse misstates the 1σ radius by −3.7% … +3.2%**
  (Mahalanobis d² over the profiled 68% contour runs 2.128–2.444 vs 2.296),
  ≈ 8× the 0.4% interpolation floor, so it is real but small — the analytic
  lumi-lognormal estimate was ~2% and nPDF's non-parabolic direction (post-fit
  width 1.08 > 1, see below) supplies the rest. **BUG FOUND ON THE FIRST REAL
  FILE:** `xsec_contour.C` bound a `double` to `deltaNLL`, which combine writes
  as **Float_t** (as it does sigmaW and r_Z), so `SetBranchAddress` returned a
  mismatch code and the macro reported "no scan". The synthetic Gaussian
  fixture could not catch it — it was written with the types the reader
  assumed. Now read through `TLeaf::GetValue()`, which converts from any
  numeric storage type. **Lesson: a synthetic fixture validates logic, never
  an interface.**
  Regression-checked locally: with no companion file present the extraction is
  BIT-IDENTICAL to the pre-change one. **`sync_lxplus.sh` audited the same
  day** — every file the remote side reads is covered (the 4+4 inputs, their
  `_systs.txt` sidecars incl. the `#! muTrig corr` directive, and the new
  `skim/output/gen_xsec_fid.txt`, whose upload destination matches the path
  `--contour` derives from `$PO_PLOTS/../../`), and the download now pulls
  `simfit/contour/` too. Four defects fixed: the post-upload hint labelled
  `./run_pO_fits.sh both simfit --asimov` as "backup: PF MET" although
  `--disc` has defaulted to leppt_mt40 since 2026-08-16 (so it ran the
  primary variant twice and never MET — it needs `--disc met`); the
  post-download hint still pointed at the LEGACY per-flavour
  `charge_asym.C`/`FBratio.C` calls instead of `run_observables.sh`; the
  fitDiagnostics loop re-used `sfx`, the OUTER discriminant loop's variable
  (harmless only because bash pre-expands a for-loop's word list — renamed to
  `kind`); and `-h` printed a FIXED line range `2,29p` that had already been
  outgrown, now `3,/^# ===/p`. **`QCD_MODE` also defaults to `abcd`
  now** (disc-dependently — see "Downstream fit"), so the refit no longer
  needs the env var: `./run_pO_fits.sh both simfit --disc leppt_mt40
  --asimov`. **NOT yet run.**
- **LHE FREAK-WEIGHT GUARD — DONE 2026-09-15** (`pOLhe::kMaxMemberRatio = 10`
  in `skim/lhe_index.h`; full details in that bullet and in the
  `syst_shapes.C` one). The user asked why `summary_qcdScale_fb` showed W⁻ y5
  at ±1.26% while every other bin sat at ±0.3–0.5%, in both qcdScale and nPDF
  and in exactly one lab bin (y6) and one fb bin (y5). **Cause: ONE event** —
  `July_29_MC_Wm_mu` entry 70485, reco muon pT 36.5 / η −0.23 (so one pT bin,
  one rapidity, two binnings), carrying the file's largest weight ratios
  (μF×2 ρ = 228, nPDF baseline product 261) and supplying **45% of that
  template's entire qcdScale shift**. Mechanism: POWHEG unweights, so |w₀| is
  a single number and the near-cancellation of B̄ that makes the ratio explode
  is invisible in the nominal weight — it only shows in ρ = ttbar_w[k]/w₀,
  and it shows in EVERY member at once (α_s ± 0.001 giving ρ = 11.5 is the
  proof). Guard neutralizes such an event's variations (ρ → 1) while keeping
  its full nominal weight; fires 41–69×/file (0.0035–0.0058%), costs ≤0.082%
  on any inclusive member sum. **Chain re-run the same day:** MC re-skim
  24/24 OK → `compare_reskim.sh rootfile_pre_rhoguard` (every DIFF is a
  member twin, **0 nominals**, data files IDENTICAL) → `run_lhe_updown.sh
  all` (maxdev0 = 0 everywhere) → all four Combine inputs regenerated
  (**sidecars unchanged ⇒ NO fork change**) → `run_syst_shapes.sh all`.
  Result: y5_FB qcdScale +1.26% → **+0.37%**, nPDF +1.40% → **+0.52%**, level
  with the neighbours; muon signal bins >10% went nPDF 1 → **0**, qcdScale
  6 → 4. **NOT yet refitted (lxplus)** — folds into the pending refit; the
  templates moved by ≤0.9% on any real template, so nothing else changes.
- **MUON TRIGGER SF BACK TO INCLUSIVE — DONE 2026-09-21 (user decision after
  discussion; reverts the application, NOT the measurement).** `kTrigBinning`
  = `kTrigInclusive`, so every muon gets the one number **0.9971
  (+0.0020 −0.0024)**. **The η diagnostics were deliberately preserved so the
  choice can be justified** — that was the explicit requirement — and they are
  regenerated on every run, not frozen copies: `muon_sf.h` still READS the
  6-bin |y| table and prints it in the `[SF]` block under
  `per-bin table (NOT applied)` with its pulls, plus the signed-y folding
  cross-check; `correction/trig_eff_mb.C` still writes `eff_absy_<sel>`,
  `eff_absy_<sel>_charge`, `eff_y_<sel>[_charge]`, the `sf_absy_/sf_y_`
  graphs, `trig_eff_<sel>.csv` and the full `BINNING DECISION`
  likelihood-ratio block (verified by re-running it in plots-only mode). Its
  labels were reworded from "THE APPLIED binning" to "measured; an inclusive
  SF is applied" so a plot cannot be misread. **The physics finding stands
  and is stated as such in both files:** the η dependence is real at ~3σ
  (LRT p = 0.0006 mt40 / 0.0004 nom), pT is flat (p = 0.71), there is no
  pT×η interaction (p = 0.68), folding is justified (p = 0.39) and 3 bins
  would not be enough (p = 0.0043) — MC itself dips at |η| < 0.4 (ε_MC 0.985
  vs 0.994, the η ≈ 0 barrel wheel gap) and data dips further (0.967).
  Being inclusive therefore means **leaving a measured ~2.4% rapidity effect
  uncorrected**, and the docs say exactly that rather than claiming flatness.
  In particular the Clopper-Pearson pull χ² (11.5/11 p = 0.41 signed, 11.5/5
  p = 0.04 folded) must NOT be quoted as the justification — it is ~16×
  under-powered near ε = 1 — and both the header and the macro carry that
  warning in place. **What it costs / buys:** the inclusive normalization is
  unchanged either way (⟨SF⟩ 0.98794 → 0.98802, ⟨TRIG⟩ 0.99699 → 0.99708), so
  the effect is pure rapidity shape — W templates ×1.0156 (|y| < 0.4) to
  ×0.9902 (|y| > 2.0), |y|- and charge-symmetric to all printed digits (σ_incl
  and the charge asymmetry untouched, dσ/dη and R_FB move) — and `muSF`
  shrinks ±0.762% → **+0.36/−0.39%**, now flat in rapidity (+0.32…+0.40%, the
  residual being ID×ISO), because muTrig goes back to its single CP error
  ±0.2% instead of six independent errors added coherently. Full muon re-skim
  (Wmu + Zmm, 14/14 jobs) + downstream re-run done: backup in
  `skim/rootfile_pre_incltrig/`, data files bit-identical, every
  still-IDENTICAL MC histogram empty (checked: 0 non-empty), per-bin ratios
  matching SF_incl/SF_bin to 5 decimals, `run_lhe_updown.sh` Wmu+Zmm
  (`max|member0 − nominal| = 0`), `run_qcd_abcd.sh mu` (A40 132.4/127.0,
  κ 1.09 — unchanged), Combine inputs regenerated with **sidecars unchanged ⇒
  no fork change needed**, `run_syst_shapes.sh` both discs (closures 1.6e-8 /
  7.6e-8, member 0 ≡ `signal`). Switching back = that one constant + the same
  sequence. **NOT yet refitted (lxplus)** — same commands as the 09-14 refit.
  Carried over from 2026-09-15 and unchanged: the MC leg of `trig_eff_mb.C`
  requires a prompt-W gen match on the leading lepton (charge-blind, via the
  EventTree `mc*` block's `mcMomPID`) — 0.002–0.004% of events, no number
  moved; and the Z uses the **standard OR formula**
  `[1−Π(1−ε_data,i)]/[1−Π(1−ε_MC,i)]` on the path-fired (not matched) per-leg
  efficiencies, whose pT>25-derived value applied to 10–25 GeV Z legs is
  harmless for muons (flat from 10 GeV) but will not be for electrons.
- **SF APPLICATION PHASE — step 2 DONE 2026-09-14: MUON SFs FOLDED INTO THE
  MC SKIM WEIGHT** (`skim/muon_sf.h`, see the Stage-1 bullet): the pp POG
  TightID + TightPFIso SFs (η × pT) and our MB-derived trigger SF (ONE
  inclusive value 0.9971 +0.0020/−0.0024 — user decision: the per-y values
  are consistent with a flat SF, χ²/ndf 11.5/11) on the leading muon (W) /
  both legs + the two-leg trigger factor (Z); **the Z→μμ iso cut was
  harmonized to the W's 0.15 (was 0.20) so one iso SF serves both; Z→ee
  already at 0.095.** ⟨SF⟩ W⁺ 0.988 / W⁻ 0.988 / DY-in-W 0.987 / Z→μμ 0.983;
  **ONE combined nuisance `muSF`** (user decision: the three independent
  sources muID ±0.05% / muIso ±0.28% / muTrig +0.20−0.24% added in quadrature
  per bin → +0.35/−0.38% on the W templates, ±0.52% on the Z peak; the
  per-source twins stay in the skim files as diagnostics), treated in the
  cards exactly like the theory shape nuisances; the LHE member weights
  rescaled so the theory variations sit on the SF-weighted nominal.
  Downstream re-run the same day: `run_lhe_updown.sh Wmu Zmm`,
  `run_qcd_abcd.sh mu` (A40 μ⁺ 132.4 ± 10.6 / μ⁻ 127.0 ± 10.2, was
  131.9/126.6; κ 1.09 unchanged), `mtandmet` / `dileptonpeak` μ (inputs +
  4-line sidecars; nominal templates bit-identical across the three-vs-one
  nuisance switch; Z peak data 388 → 364, `Z_incl/signal` 372.1 → 355.4),
  `run_syst_shapes.sh all`, `njet_WZ Zmm` (364 = event-identical), fork
  dry-runs OK (per-flavour shape rows, `lepsf group = muSF`, 0 edit lines).
  NB the iso-cut change
  removed 6.2% of the Z data but 2.9% of the DY MC — the data iso efficiency
  sits below the unembedded MC's, opposite to the pp iso SF's +1.6%/pair: a
  pO tag-and-probe iso measurement is the natural follow-up. Step 1
  (2026-09-09) was the trigger turn-on
  itself (`correction/trig_eff_mb.C`: MB path verified in data — run-dependent
  prescale, none in runs 393974/393975 — and MC 99.9%; pT > 25 SF μ 0.996 ±
  0.002, e 0.989 ± 0.003 plain / 0.997 ± 0.004 with m_T > 40, the difference
  = the data sample's fake composition). **REFIT DONE on lxplus 2026-09-14
  night** (abcd/leppt_mt40, all four shape nuisances; status 0 / covQual 3,
  Asimov closure exact; pulls lab: nPDF +0.51 ± 0.98, qcdScale +1.16 ± 0.93,
  alphaS −0.03, **muSF −0.54 ± 0.99**, lumi 0.01; r_Z 1.087 ± 0.054; QCD
  multipliers μ 0.897/0.925, e 0.820/0.651) — **σ = r×σ_gen: W⁺ 58.96 ± 0.98
  (stat) ± 1.82 (syst) / W⁻ 47.51 ± 0.87 ± 1.48 / W 106.47 ± 1.31 ± 3.25 nb**
  (eff. r 1.185/1.200/1.191; back at the 08-24 level, so the 09-07 −6% was
  the raw-variation artifact; the stat/syst split is the 2026-09-15
  covariance conditioning, see "STAT/SYST SPLIT" under "Downstream fit" —
  lumi 3% alone is 3.19 of the 3.25 nb on W, every other profiled nuisance
  together ≈ 0.6 nb; per bin stat dominates, stat/total 0.85–0.91).
  Inclusive postfit χ²/ndf μ 3.33 / e 7.17.
  **FINAL-PLOT UNCERTAINTY PRESENTATION DONE 2026-09-15** (user: "data point +
  stat error bar + syst as boxes using TBox"): the chain was audited end to
  end and **no refit is needed** — the 09-14 lxplus fit was already
  re-extracted locally that night with the `ComputeStatCov` extractor, so
  `comb_W_yields.csv` carries `rErr_stat` (19 columns) and
  `comb_fitted_yields.root` carries `h_cov_yield[_FB]_stat`; the only gap was
  downstream, in `charge_asym.C` / `FBratio.C` / `observables.C`, which read
  the TOTAL covariance only. Now every observable draws point + stat bar +
  `TBox` syst (`plotting_helper.C::MakeSystBoxes`, single source; re-synced to
  the fork): σ W⁺ 58.96 ± 0.98 (stat) ± 1.82 (syst) / W⁻ 47.51 ± 0.87 ± 1.48
  / W 106.47 ± 1.31 ± 3.25 nb; A_ch syst 0.0026–0.0099 vs stat 0.037–0.047;
  R_FB syst 0.001–0.009 vs stat 0.058–0.105. **The near-vanishing syst on
  A_ch and R_FB is the validation, not a disappointment** — lumi 3% is the
  dominant nuisance and cancels exactly in a ratio of coherently-scaled
  yields, while on σ it is 3.19 nb of the 3.25 nb W-inclusive systematic.
  User decisions: lumi stays INSIDE the syst box (separating it needs a
  `lumi`-only conditioning variant of `ComputeStatCov`, one line, not
  implemented), and boxes are drawn on every observable, not just σ. NB the
  `met` variant (`pO_fit_out/`, Aug-6 extraction) has NO `h_cov_yield_stat`
  and its fit is pre-filter-era stale anyway, and the legacy per-flavour
  `{mu,ele}_fitted_yields.root` are Aug-3 with no covariance at all — both
  fall back to a single total bar with a `[WARN]`. Only the `leppt_mt40` comb
  plots are live.
  **NEXT:** the electron SFs (same pattern;
  `trig_eff_mb.C` ele exists, ID/ISO need a pp file or pO tag-and-probe) and,
  still open, the `skim_Zmm` iso cut 0.20 vs the Tight-SF WP and the
  pp-SF-in-pO caveat (no pile-up in pO).
- **LHE WEIGHTS → nPDF SYSTEMATIC, SKIM STAGE DONE 2026-09-07** — every MC
  skim file now carries the member twins `<h>_epps21/_scale/_alphas` of the
  fit templates (`skim/lhe_index.h`, filled in `skim.C` with
  `w·ttbar_w[i]/ttbar_w[0]`; member 0 ≡ nominal; MC-only re-skim verified
  bit-identical on every pre-existing histogram) and the three Up/Down pairs
  per template written by `skim/lhe_updown.py`: `<h>_nPDFUp/Down` (LHAPDF
  `PDFSet.uncertainty()`), `<h>_qcdScaleUp/Down` (μR/μF envelope over all 9
  points, per-bin max/min incl. the nominal) and `<h>_alphaSUp/Down` (the
  α_s 0.119/0.117 member templates) — `run_lhe_updown.sh`, env `lhe_env.sh`;
  LHAPDF 6.5.6 built under `~/local/lhapdf` for the PyROOT python 3.12. User
  decisions: raw variation (no gen-level division), 107 EPPS21 members + 14
  scale/α_s weights stored (not 217), the official LHAPDF function for the
  nPDF combination, the all-points envelope for the scales, the two ±0.001
  members for α_s. **Combine inputs DONE the same day:** `mtandmet.C` /
  `dileptonpeak.C` carry `<process>_<syst>Up/Down` for `signal/z/ztau/wtau`
  (`signal/ztau/w/wtau` under the Z peak) into every region of
  `combine_input_W[_leppt_mt40].root` and `combine_input_Z.root` (nominal
  templates bit-identical) plus the `<input>_systs.txt` sidecar listing the
  systematics/processes written (see "Structured inputs"). **Fork side DONE the same evening (local checks only):**
  `make_pO_simfit_cards.sh` reads the sidecars → three `shape` rows + fifth
  `shapes` token + `lhe` group (`LHE_SYST=auto|off|list`; `off` ≡ baseline
  cards + 1 comment line; alignment verified 248/324 columns), driver copies
  the sidecars with the inputs and passes `lheSysts` to the extractor
  (pulls/constraints rows, Asimov closure θ=0, safety-net sweep of unreported
  floating params), `sync_lxplus.sh` uploads the sidecars and downloads the
  `fitDiagnostics_simfit_<B>.root`, `postfit_incl.C` reads `shapes_fit_s` when
  the fitted sidecar lists shape nuisances (loud approximate fallback).
  **FITTED on lxplus the same night (abcd mode, leppt_mt40, 19:53–19:57;
  downloaded incl. the fitDiagnostics files; observables + postfit_incl via
  `shapes_fit_s` rerun):** status 0 / covQual 3 both variants; Asimov closure
  EXACT (all 25 POIs = 1, all 12 scales = 1, every θ ≈ 1e-7). Pulls
  (lab | fb): **nPDF +0.32 ± 1.08 | +0.14 ± 1.43** (unpulled, unconstrained;
  the fb post-fit width > 1 = non-parabolic likelihood along θ, from the
  asymmetric templates), **alphaS −0.04 ± 0.99 | −0.11 ± 1.00**, **qcdScale
  +1.30 ± 0.85 | +1.71 ± 0.86** — the lepton-pT shape constrains μR/μF and
  the data prefer the upper envelope (χ²/ndf of the inclusive postfit stacks
  improved μ 5.96 → 3.07, e 13.2 → 6.7; the e tail excess and the 25–32
  turn-on deficit remain). QCD multipliers unchanged (μ 0.912/0.935, e
  0.837/0.656, θ_e⁻ still −2.7σ). r_Z 1.038 ± 0.065 (was 1.095 ± 0.053).
  **σ = r×σ_gen(nominal): W⁺ 55.24 ± 3.26 / W⁻ 44.63 ± 2.63 / W 99.87 ± 5.75
  nb (eff. r 1.110/1.127/1.118) vs 58.71 ± 2.03 / 47.44 ± 1.68 / 106.15 ±
  3.44 before — the errors grew 3.5% → 5.9% (theory ⊕ ≈ 4.7%) as intended,
  BUT the central values dropped 6%, which is NOT physics: with the RAW
  variation the pulled qcdScale supplies κ⁺^1.30 = 1.050 (×1.0078 nPDF) of
  normalization to every MC template (verified: postfit signal = r × prefit
  × 1.0575 in y5, = the product of the four κ^θ), r compensates, and
  σ_meas = N/(L·A·ε) × σ_gen(0)/σ_theory(θ̂) inherits the theory-normalization
  pull.** DO NOT quote these σ. The likelihood with the shape nuisances
  (templates, Poisson × Gaussian form, the κ(θ) normalization and the
  quadratic/linear shape morphing, the pull/constraint definitions, this
  fit's numbers and the acceptance-only remedy, Eq. by Eq.) is written up
  AN-style in `docs/AN_lhe_systematics.tex` (2026-09-08; fragment like the other
  `AN_*.tex`, `\input` it — `eq:infit` resolves inside the AN). **2026-09-09:
  + subsection "What a PDF member is allowed to change in the measured cross
  section"** — the step-by-step derivation that r_y×σ_gen,y ≡ (N_y−B_y)/(L·(A·ε)_y)
  (σ_th, A, N_gen cancel = the pPb corrected-count formula), the per-member
  ratios ρ_sel,y(k)/ρ_fid,y(k), the consistent propagation
  σ_meas(k) = σ_meas(0)·(A·ε)(0)/(A·ε)(k) (only A·ε survives — the quantity the
  pPb HIN-17-007 `MC_Syst_PDF` varies, 0.01–0.09%), what the raw variation
  makes the fit compute instead (r_y(θ) = r_y(0)/ρ_sel,y(θ) with σ_gen fixed ⇒
  the full theory-σ PDF uncertainty lands on r, not an oxygen-PDF property),
  the raw-vs-shape-only table (raw +2.5/−3.4% central, +5.7/−8.2% y11 vs
  shape-only residual 0.4–0.7% per pT bin, from area-normalizing the stored
  `_epps21` members before the LHAPDF combination) and the acceptance-only
  recipe X̃(k) = X(k)/ρ_fid,y(k) (needs the gen twins). **PDF:
  `docs/AN_lhe_systematics.pdf` (4 pp.), built by `AN_lhe_systematics_wrapper.tex`
  in `docs/`** (moved from the repo root 2026-09-27; `pdflatex AN_lhe_systematics_wrapper.tex` ×2 from
  `docs/`; provides amsmath, the CMS macro stubs \PW/\PZ/\Pgm/\Pe/\pt/\GeV/
  \pPb, a `\newlabel` stub for `eq:infit` and the HIN-17-007 bib entry).
  **DECIDED + APPLIED 2026-09-14 (user's final version, after the 09-10..14
  discussion that the r uncertainty must not carry the theory-σ uncertainty —
  a count-based σ = (N−B)/(L·Aε) never contains σ_th, and r×σ_gen equals it
  only if the same σ_th enters the template and σ_gen; "normalize to the
  reco nominal" = area-normalize, "to the gen nominal" = Way 1 with the gen
  twins, the difference being the Aε change, ≲0.1% for PDF/α_s, possibly ~1%
  for scales): every member template is area-normalized to the nominal
  integral before combining (all three families; the Aε part deferred),
  scales = 6-point envelope (antagonistic corners dropped), α_s symmetrized
  ±(N[0.119]−N[0.117])/2 — implemented in `lhe_updown.py` (`--norm reco`,
  `--scale-points 6`, `--alphas-mode symm` defaults), Up/Down rewritten in
  all 24 MC files, all Combine inputs regenerated (nominals bit-identical),
  AN paragraph "Implemented choice" + PDF rebuilt. Stored shifts now nPDF
  ±0.4% / scales ±0.3% / α_s 0. **NEXT = refit on lxplus** (`sync_lxplus.sh
  upload`, `./run_pO_fits.sh --disc leppt_mt40 --asimov`, download, then
  `run_observables.sh leppt_mt40` + `postfit_incl`): expect the θ pulls to
  move shapes only, σ back at the 08-24 level with a much smaller theory
  contribution, and the qcdScale pull to be re-examined (it can no longer
  buy normalization). **The ECAL-gap veto was then DONE the same day
  (2026-09-14, reco-only on `eleSCEta`, gen fiducial kept at |η| < 2.4 —
  user decision after the assessment; electron re-skim + full downstream
  re-run, see selection step 5): the electron Combine inputs the refit will
  use carry BOTH changes (area-normalized theory variations + the crack
  veto), so the next fit's e-channel r's and the crack-bin A·ε are not
  comparable one-to-one with the 09-07 fit.**
  **Earlier NEXT list (partly superseded by the above): (1) the
  gen-level twins in `gen_xsec.C` + the acceptance-only variation (divide each
  varied template by σ_gen,i(k)/σ_gen,i(0) so the nuisance moves only A·ε and
  shape) — or equivalently evaluate σ_gen at θ̂; (2) consider replacing the
  per-bin envelope by separate μR / μF nuisances (physical members idx 3/6,
  1/2; the envelope "Up" is a per-bin composite, so θ = +1.3 has no single
  scale interpretation); (3) then re-fit, plus `LHE_SYST=off` for the
  breakdown, impacts; (4) per-eigen-direction nPDF (collapse inflation
  ×1.27/1.16 measured by `syst_shapes.C`), pairing verification.**
- **IN-FIT ABCD QCD (2026-08-23, `QCD_MODE=abcd`)** — the leppt_mt40 QCD
  normalization moved into the simfit likelihood, per the Mattermost decision
  (Andre's (m_T × iso) region layout — pT is never an ABCD axis — executed
  with Combine's rateParam-formula ABCD so the EWK subtraction floats with the
  fitted r's, killing the prefit-r circularity; the SHAPE stays the anti-iso
  m_T>40 pT template). Implemented end-to-end: `qcd_abcd.C` (C40 region +
  `abcd_counts_*` export + in-fit window scan + reduced κ + the per-pT-bin
  fake-factor diagnostic `runFFCheck`), `mtandmet.C` (`qcd_abcd` 7th template
  + 6 CR dirs), fork (`QCD_MODE=abcd` cards + maps, extractor `AbcdScale`,
  closure extension, card-gen failure guard), `postfit_incl.C` (qcd_model
  aware). All local checks green (lnN/free cards byte-identical; Σqcd_abcd =
  A0; CR contents ≡ exports; observables chain regresses). **e-FF finding:
  F(pT) rises with pT, total +13% — the direct e iso-pT correlation → κ_e
  1.15.** **VALIDATED on lxplus 2026-08-24**: Asimov closure EXACT (<1e-7,
  all 25 POIs + 12 scales + θ's), fit status 0 / covQual 3; multipliers
  κ^θ·sB·sC/sD = μ 0.894±0.102 / 0.919±0.106 (healthy — driven by sB≈0.94,
  the floating subtraction at r̂>1; θ −0.4/−0.3σ), e 0.833±0.067 /
  **0.654±0.053 with θ_e⁻ = −2.6σ** (fb variant same pattern:
  0.904/0.913/0.734/0.565); r_Z 1.095±0.053. **Postfit diagnosis (2026-08-24,
  era-consistent required-QCD reconstruction)**: (i) e required/template
  ratio rises 0.4→3 across pT 25–100 for BOTH charges (μ control flat) —
  the flat-T shape is too soft = the FF finding in the fit; FF-weighted
  prediction reproduces the e⁺ 34–60 demand (338 vs 312; flat-T 244);
  (ii) the 25–34 deficit is FLAVOR-UNIVERSAL (all 4 channels) → the missing
  turn-on TnP SFs, not QCD — **HALF-SUPERSEDED 2026-09-15: after the muon
  ID+ISO+trigger SFs (09-14) the μ deficit SHRANK but did not close (−3.4%
  over 24–36 GeV vs the electron's −7.6%, and μ recovers one bin earlier),
  so it is no longer flavor-universal in magnitude; the leading μ suspect is
  now the lepton momentum scale/smearing, which is still not applied and
  migrates events across the 25 GeV cut — see the `postfit_incl.C` bullet
  for the per-bin table**; (iii) the +/− split = charge composition: data
  ratio 1.44±0.10 at [25,30) vs model 1.21, sideband predicts 0.97 (C40
  e⁻-rich at 2.3σ) — iso-pass fakes e⁺-rich, anti-iso e⁻-rich; (iv) the
  >60 GeV excess has W-like charge ratio ~1.4 → likely non-QCD (scale tail /
  Wγ/conversions). The e multipliers = the only per-(e,charge) norm dof →
  absorb all of this. **OPEN: (a) the lnN→abcd default flip for leppt_mt40
  (validated, decision pending); (b) the B/D boundary question — user will
  revisit whether to tile at m_T<40 (no buffer; moves A0 by only +2.0–2.6%:
  131.9→135.3, 126.6→129.1, 652.1→667.1, 670.7→688.1, B W-contamination
  +~50% rel.) vs the implemented m_T<30 + 30–40 buffer; (c) follow-ups:
  FF-shaped e template variant, TnP SFs at threshold, >60 GeV e excess
  composition.** AN tex + pdf updated 2026-08-24 (in-fit definitions +
  likelihood eqs, region map with the buffer stated explicitly, results
  table, FF figure; pdf rebuilt via a standalone wrapper — external \ref's
  to other AN sections render as ??). Details in the `qcd_abcd.C` bullet +
  "Downstream fit".
- **GRAND SIMULTANEOUS FIT — simfit (2026-08-04, NEW DEFAULT)** — the fit
  procedure changed completely: one likelihood per binning variant (lab/fb)
  with all 48 (flavour, charge, y) W channels + both Z peaks; 25 POIs
  (`r_<C>_y<i>` μ/e-shared + global `r_Z`); implemented end-to-end (fork:
  `make_pO_simfit_cards.sh`/`extract_pO_simfit.C`/driver `simfit` mode +
  `--asimov` closure; repo: covariance-aware `charge_asym`/`FBratio`,
  `observables_comb`, run_observables simfit chain). **Validation pending:**
  needs the first lxplus run (`./run_pO_fits.sh --asimov`) → check Asimov
  closure PASS, fit status 0/covQual 3, r's vs the legacy per-bin values, and
  `r_Z ≈ 1`. Legacy per-bin pipeline intentionally kept runnable ("wrong way"
  reference: it re-fits the same Z data per bin). Details in "Downstream fit".
  DY vetoes included) and both isolation studies now use the Δβ PU-corrected
  relIso; cuts 0.15/0.095 re-confirmed as optima. **Requires a full re-skim
  (all channels, all samples) + downstream re-run** (ABCD QCD, plots, Combine
  inputs) — existing `skim/rootfile/` outputs still carry the old definition.
- **Electron isolation + background** — working point confirmed
  (`eleMVAIdWP95` + Δβ relIso < 0.095); PU correction for electron isolation
  DONE via the Δβ switch above.
- **Tau channels as background** — sample enums (`kWptau`, `kWmtau`, `kDYtau`)
  and `Wptau/Wmtau/DYtau` filelists exist; selection code largely there but
  validation pending. Commit `447a2b5` ("before adding tau") marks the boundary.
- **Corrections still WIP** — MUON efficiency SFs (ID / ISO / trigger) ARE
  applied since 2026-09-14 (`skim/muon_sf.h`, muon MC only); electron SFs
  and the lepton momentum scale/smearing are not applied yet. Don't assume MC
  is fully corrected when reading skim output. (The per-event *generator* weight IS now
  applied — see "Generator event weight" below — but these reco-level
  corrections are separate and still missing.) **Recoil: checked
  2026-07-02 via the raw-recoil study (`correction/recoil_raw[_ele].C`)
  and looks stable for now — no recoil correction applied; `recoil_fit.C`
  deferred. The MET data/MC shape discrepancy is understood to come from
  the current MC samples being UNEMBEDDED (no pO underlying event embedded
  in the simulation), not from recoil miscalibration.**
- **ABCD QCD background (muon + electron, done 2026-06-23)** — replaces the
  Rayleigh shape-extrapolation. Both `skim_Wmu`/`skim_Wel` store the relIso × MET
  (and × m_T) 2D per charge (`h_iso_{met,mt}_{mu,ele}{Plus,Minus}`);
  `correction/qcd_abcd.C[+(true)]` does `N_A = N_B·N_C/N_D` on QCD-only
  (EWK-subtracted) counts and emits the iso-pass QCD template, wired into the
  `mtandmet.C` MET stacks for both channels. Electron QCD is ~10× the muon and
  is the dominant low-MET background there. Remaining: use these templates in
  the Combine fork datacards.

Future plans / on the TODO list (not yet started):
- **σ extraction as r × σ_gen-fid — DONE (2026-08-12)** — the measured cross
  section is now quoted as fitted-signal-strength × gen fiducial σ instead of
  N_fit/(L·A·ε); implemented in `xsec_fiducial_comb` (see the Stage-2 bullet).
  Algebraically identical when A·ε comes from the same MC
  (r×σ_gen = N_fit/L/(A·ε)_MC), but kA_O and kSigma cancel between r and σ_gen
  (result independent of the assumed theory σ and the A=16 scaling), and future
  SFs on the reco templates propagate into r automatically. Counts CSVs
  deliberately stay count-based (the fit's record); charge_asym/FBratio did
  too until 2026-09-22, when the user ruled that EVERY observable uses
  r × σ_gen (the "A×ε cancels in a ratio" argument holds for A_ch, not for
  R_FB — see the `fiducial_yields.C` bullet). Residual model dependence
  concentrates at the pT=25 threshold (lepton scale — the known e-scale shift
  matters here) and in-fiducial efficiency.
- **ECAL transition gap veto — DONE 2026-09-14 (TODO added 2026-08-12)** —
  1.4442 < |η_SC| < 1.566 excluded in the electron RECO selection, data and MC
  (skim_Wel leading electron + DY-veto legs, skim_Zee both legs, and the four
  replicating macros), see selection step 5. Decision taken: the electron GEN
  fiducial keeps the common |η| < 2.4 (reco-only veto; the gap extrapolation
  is an intra-bin shape effect, to be quoted as a small systematic). Electron
  re-skim + downstream re-run done the same day (numbers under step 5);
  fit NOT yet re-run.
- **Tau cross-section measurement** — currently τ is only a background to W/Z.
  Adding W→τν / Z→ττ as their own signal channels is on the wishlist (tagged
  as "probably just an attempt").

Known FIXMEs (mostly in the electron path):
- ~~`DrawDielectronPeak.C` — ECAL coverage comment says 2.4 but should be 2.5
  with the 1.4442–1.566 transition gap excluded~~ — RESOLVED 2026-09-14: the
  gap is vetoed per leg on `eleSCEta` in `skim_Zee` (and everywhere else an
  electron is ID'd); |η| < 2.4 is kept deliberately as the common μ/e fiducial
  edge (the legacy macro's comment is left as is).
- **W→eν electron ID: DONE (2026-06-24)** — the isolation/ID study
  (`correction/isolation_ele.C`, now with a continuous-relIso scan) showed the old
  `eleCutIdWP95` was the *worst* QCD rejector; `skim_Wel` now uses **`eleMVAIdWP95`**
  (leading lepton + DY veto), continuous relIso < 0.095 kept (Youden-J optimum). This
  cut the electron ABCD QCD ~64% (≈10037→3634 iso-pass; signal −6.5%).
  **Electron gates HARMONIZED (2026-07-02):** `skim_Zee` switched to `eleMVAIdWP95`
  too (both legs; Z-window data 276→247, −10.5%, consistent with the tighter ID),
  and the `skim_Wel` DY veto's iso gate switched from the integer `eleMVAIsoWP95`
  WP to the **continuous `RelIsoPF < 0.095`** (the study characterized the
  continuous cut; the MVA-Iso WP was shown not to track it). Every electron ID
  gate in the skim is now `eleMVAIdWP95` and every iso gate is continuous relIso.
  Requires an electron re-skim (Wel + Zee, all samples) + downstream re-run.
- **Electron ΔβPU correction: DONE (2026-07-06)** — all relIso (both flavours,
  skim + studies) now Δβ-corrected; the 0.095 optimum was re-confirmed under the
  new definition (J_MET optimum exactly 0.095, AUC_MET 0.909). As anticipated,
  the correction changes very little at pO pileup.

Not in repo: no tests, no CI, no config files. Many recent commit messages are
placeholder ("xx") so `git log` is not a reliable narrative — read the diffs.

## When working on this repo

- Don't expect type-checking or a test suite. Sanity check by re-running the
  relevant `skim/run_all.sh <channel>` invocation and reading `skim/logs/`.
- **Run stages through their logging wrapper, never the bare macro, whenever the
  console output is itself a deliverable.** Today that means
  `correction/run_qcd_abcd.sh` for the ABCD QCD estimate (NOT
  `root -l -b -q 'qcd_abcd.C+'`), `correction/run_charge_flip.sh` for the
  charge-flip study, `correction/run_isolation.sh` for the isolation scans,
  `correction/run_trig_eff_mb.sh` for the trigger turn-on,
  `plotting/run_combine_inputs.sh` for the Combine inputs (2026-09-21; NOT
  `root -l -b -q 'mtandmet.C+(false)'` — see the Stage-2 bullet),
  `plotting/run_syst_shapes.sh` for the shape-systematics diagnostics,
  `skim/run_lhe_updown.sh` for the theory Up/Down templates,
  and `skim/run_all.sh` for the skim. The ABCD
  report — region composition, closure, factorisation tests, window scan,
  two-plane transport + r-scan, systematic budget → lnN κ — is the source of the
  numbers quoted in `docs/AN_qcd_background.tex`, and ROOT prints it to stdout only,
  so a bare run silently discards the record. If a new macro starts producing
  numbers that get quoted anywhere, give it a wrapper too (`mkdir -p logs`,
  pre-build once to avoid the ACLiC race, tee per job, echo the headlines).
  **Corollary (user, 2026-09-21): let the wrapper name its own log.** Every one
  derives the path from its arguments (`logs/skim_${ch}_${s}.log`,
  `logs/qcd_abcd_${chan}.log`, `logs/syst_shapes_${d}.log`, …), which is what
  makes any stage's record findable months later without reading the code that
  wrote it. A one-off `> some.log` redirect invented inside a driver breaks
  that for no gain and is orphaned when the driver is — drivers written during
  a session live in a scratchpad, not the repo. A macro with no wrapper should
  GET one rather than an ad-hoc name; that is exactly why
  `run_combine_inputs.sh` exists.
- After the skim refactor, common scaffolding lives in `skim_common.h` and
  the four channel selections live as functions in `skim/skim.C`. Most
  cross-channel changes (new cut, new histogram) can be made once in the
  header; only physics differences need to go in the per-function bodies.
  When in doubt, diff against `skim/legacy/` — that's the pre-refactor
  reference.
- Histogram naming is load-bearing — `analysis/charge_asym.C` and
  `analysis/FBratio.C` read by name, and so do the Combine input scripts.
  Rename in skim → must rename in analysis AND in the Combine fork.
  **Naming audit (2026-07-30):** names must state the quantity they hold.
  Fixed then: the ABCD region printout said "MET-high/low" even on the m_T
  plane (`ABCDConfig::metCut` → `yCut`, labels now derived from the plane);
  the fitted-yield containers were called `h_mt_*` although the fit
  discriminant is PF MET → the fork now writes **`h_yield_W{p,m}_y*(_FB)`**
  as primary with `h_mt_*` kept as a deprecated alias, and
  `charge_asym.C`/`FBratio.C` prefer `h_yield_*` and fall back (so old and new
  files both work in both directions). Same dual-name scheme for the output
  graphs: **`g_chargeAsym`, `g_RFB_{sum,Wp,Wm}`** primary + `_mt` alias;
  `observables.C` resolves via a `GetGraph(file, name, legacy)` helper and its
  plots are now `chargeAsym.png` / `RFB_{sum,Wp,Wm}.png`. `useMT` still
  genuinely selects m_T vs MET **for raw skim files only** — it is ignored when
  `h_yield_*` is present. Also fixed: a real crossed assignment in
  `correction/PlotsIsoROC.C` (the SS and MET efficiency graphs were drawn under
  each other's legend entry).
  **Follow-up (2026-09-21): the untagged plot outputs that naming audit
  orphaned were DELETED** — `plots[/Elec]/charge_asym/chargeAsym_mt.{png,pdf}`,
  `plots[/Elec]/FBratio/RFB_mt_{sum,Wp,Wm}.*`, `plots/xsec/W_{fiducial,
  dsigma_deta_mu,dsigma_deta_ele}.*`, `plots/merged/{chargeAsym,RFB_*}_overlay.*`
  (30 files, all 2026-07-06) plus the orphan `plotting/RpO_rootfile/
  RpO_graphs.root` (only `RpO_FB_graphs.root` is ever read). They were
  **doubly obsolete**: written before the 2026-08-03 disc-tag scheme, so no
  live code path could ever overwrite them (every writer targets
  `<disc>/`), AND named after the retired `_mt` discriminant convention. If
  you find a reference to one of those paths anywhere, it is stale.
- Output paths in the skim functions are referenced again in plotting and in
  the Combine input scripts (`make_combine_input*.C`). Moving files mid-pipeline
  silently breaks downstream stages.
- `.so` / `_ACLiC_dict_rdict.pcm` files in `skim/` are ROOT's compiled-macro
  cache (now gitignored). Delete them if the macro signature changes and
  ROOT picks up the stale build.
- **NEVER hardcode an input ntuple path — every macro that reads the raw
  ntuples must take it from `skim/skim_common.h`.** `pOSkim::kDefaultDataFile`
  (data) and `pOSkim::ResolveMCSample()` (MC) are the SINGLE SOURCE OF TRUTH for
  the production in use. The user repoints them when a new production lands, and
  the whole repo — `skim/` **and** `correction/` — must follow with no further
  edits, so a macro's file argument must DEFAULT to those symbols (an explicit
  path stays available as an override argument). Enforced 2026-08-20 after an
  audit found four stragglers still on dead defaults:
  `correction/isolation_mu_tight.C` (May-26 local), `correction/isolation_ele.C`
  and `correction/isolation.C` (`pO_2025.root` on EOS — long gone) and
  `skim/PrintHiForestStructure.C`; all four now `#include` the header and
  default to `pOSkim::kDefaultDataFile`. **This is why the isolation working
  points had been quoted on May-26 data while everything else ran on July-29** —
  a silent version skew that no error message would ever have surfaced, since a
  stale-but-existing file just opens and runs. Macros reading the SKIM outputs
  (`../skim/rootfile/…`: `qcd_abcd.C`, `ptmt_scan.C`, `dataMC_kinematics.C`,
  `recoil_raw*.C`, all of `plotting/` and `analysis/`) follow the production
  automatically and need nothing — the rule is about the ones that loop the
  ntuples themselves (today: `skim.C`, `count_ngen.C`, `gen_xsec.C`,
  `njet_WZ.C`, `charge_flip.C`, the three isolation studies,
  `PrintHiForestStructure.C`).
  When adding such a macro, check it against this list. Two practical notes
  from doing this: (1) pulling `skim_common.h` into a macro puts
  `pOSkim::ComputePFMET` etc. in scope for the WHOLE ROOT session, so a macro
  with its own same-named helper will break a *different* macro that does
  `using namespace pOSkim` (this bit `njet_WZ.C` — hence the
  `ComputePFMET_ele`/`_isoleg` renames); give session-visible helpers unique
  names. (2) Verify a repoint by grepping the run log for the entry count, not
  by trusting the default: `Processing entries: 3890722` = July-29,
  `2425836` = May-26.
- **Always create output directories before writing.** Output dirs (`rootfile/`,
  `plots/`, `output/`, `logs/`) are gitignored and absent on a fresh checkout
  (e.g. a pull on lxplus), so any code that writes must `gSystem->mkdir(dir,
  kTRUE)` (C macros) or `mkdir -p` (shell) first. This is already done in
  `skim/skim.C` (creates `rootfile/`), `skim/run_all.sh` (`rootfile output
  logs`), and every `correction/` macro (their `plots/` / `rootfile/`). When
  adding a new writer, do the same — and zombie-check input files you open.
- **Never run two `run_all.sh` invocations in parallel from `skim/`.** They share
  the ACLiC build artifacts (`skim_C.so`, `skim_C_ACLiC_dict.*`) in that one
  directory, so two simultaneous first-jobs race to rebuild them and one dies
  with `Error in <ACLiC>: Executing 'rootcling ...' failed!`. The driver reports
  it only as `failures=1` in its summary line, and because the compile (not the
  physics) failed, the previous output file is silently left in place — easy to
  mistake for success. Either run the channels sequentially, or pre-build once
  (`root -l -b -q -e '.L skim.C+'`) before launching parallel jobs. Hit
  2026-07-30 with `run_all.sh Wmu` + `run_all.sh Wel` started together.
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
`pthat`, `ProcessID`, and `ttbar_w` = the 217 LHE systematic weights — scale
variations + nPDF error sets, decoded by `skim/lhe_weights.C`, see Stage 1).
All four channel functions:

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
overshoot: initially suspected uncalibrated MET/recoil, but the recoil check
(2026-07-02) came out stable — now attributed to the **unembedded MC** (no pO
underlying event in the simulation; NOT efficiency). The Z mass peak (no MET)
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
- **Isolation cut**: 0.15 in `skim_Wmu`, 0.095 in `skim_Wel` — **intentional**;
  both re-confirmed as Youden-J optima under the Δβ-corrected definition
  (2026-07-06), so the values are now settled.
- **Z-channel isolation cut = the W channel's since 2026-09-14** (user
  decision while applying the muon SFs): `skim_Zmm` moved from 0.20 (the POG
  Medium PF-iso WP, for which the pp SF file has no TightID-denominator table)
  to **0.15**, so the one `NUM_TightPFIso_DEN_TightID` SF serves the W and
  both Z legs; `skim_Zee` was already at 0.095 = `skim_Wel`. `correction/njet_WZ.C`
  replicates the value (`isoMax = isMu ? 0.15 : 0.095` in `RunZ`). Every
  Z→μμ number quoted before that date (388 Z-window data events, MC 372.1,
  the recoil/kinematics plots) is at the 0.20 cut.
- **DY-veto gate type**: ~~continuous PF relIso (muon) vs integer
  `eleCutIdWP95` / `eleMVAIsoWP95` (electron)~~ — was pre-existing; **harmonized
  2026-07-02**: both flavours now gate veto legs on continuous PF relIso
  (μ < 0.15, e < 0.095) with the flavour's tight ID (`eleMVAIdWP95` for e).
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
