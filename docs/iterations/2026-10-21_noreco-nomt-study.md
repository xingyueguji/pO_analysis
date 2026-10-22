# Iteration 2026-10-21 — electron Reco SF → one ±10% nuisance; inputs + multijet study for the no-m_T lepton-pT fit

**Status:** step 1 built and propagated to the Combine inputs (not fitted);
step 2 = skim inputs + a design study. **The user chose scheme C (§5); it is
built in `2026-10-21_leppt-nomt-schemeC.md`.** The last fitted
iteration is `2026-10-19_abcdy-mode.md`; its results stay valid for that model.
Backups: tag `pre_noreco`, §7.

## 1. What was asked

1. *"Remove the reco SFs applied on electron. Instead, we will only use iso+ID
   SF, and assign a ~10% nuisance on top of it"* — the pO data support a Reco
   SF ≈ 1 against ≈ 0.93 for pp, so the fit should tell; one inclusive
   nuisance, no rapidity or pT dependence; electrons only.
2. *"Re-use the pT fit without mT cut (but keep the current pT with mT > 40
   there for future use). I would quote the current nominal fit to be the pT
   fit, without mT cuts"* — because of the known electron PF-MET issue. The
   fiducial volume and, most importantly, the in-fit ABCD have to change; the
   user asked whether a per-rapidity-bin multijet shape with a pT-dependent
   transfer factor works, mainly for electrons, and to confirm the result first.
3. (mid-task) *Can the QCD estimate use charged PF candidates only for the MET
   (less biased; e⁺ and e⁻ behave differently in the ABCD plane), since this
   is the only place m_T / MET is used — or does it not matter?*

## 2. Step 1 — electron Reco SF

`skim/electron_sf.h`: new `kEleRecoSource = kRecoUnit` = factor 1 for every
electron, ±`kEleRecoUnitUnc` = 10% per electron (Z→ee legs coherently). The
Reco source is left out of the combined `eleSF` (new `combMask` in
`skim/muon_sf.h` `BookSFTwins`/`FinalizeSFTwins`/`SFStats::Print`; the muon
default = all sources, output unchanged) and its per-source twins
`<h>_eleRecoUp/Down` are shipped as a second nuisance, `kEleSFSystNames =
{eleSF, eleReco}` → sidecars → fork `shape` row in group `lepsf` (no fork
change needed). ID+iso (EGM wp90iso) and the trigger SF (MB study, 1.0010)
are kept.

| (skim logs) | before (EGM Reco) | now |
|---|---|---|
| ⟨SF⟩ W⁺ / W⁻ electron | 0.8865 / 0.8876 | 0.9548 / 0.9555 |
| ⟨SF⟩ Z→ee pair | 0.7879 | 0.9110 |
| eleSF on the W templates | +1.23 / −1.25% | +1.03 / −1.04% (ID+iso ⊕ trigger) |
| eleSF on the Z peak | +2.67 / −2.64% | +2.38 / −2.35% |
| eleReco | (inside eleSF) | W ±10.0%, Z +21.0 / −19.0% |
| prefit r (leppt_mt40, display QCD), e W incl | 1.163 ± 0.025 | 1.088 ± 0.023 |
| prefit e/μ W incl (W⁺ / W⁻) | 0.980 (1.031 / 0.918) | 0.917 (0.964 / 0.859) |

Chain: full re-skim (32 jobs, 0 failures) → `run_lhe_updown.sh all`
(member 0 ≡ nominal everywhere) → `run_qcd_abcd.sh both` (electron residual
κ 4.8% → still 1.05; reduced κ 1.12 / 1.14) → `run_combine_inputs.sh both`.
Regression: muon W/Z skim files, the electron data files and the Z→ee data
file bit-identical apart from the new objects (`skim/logs/compare_noreco.log`);
the four muon Combine inputs bit-identical; electron inputs: every MC template
moved, + 510 (W leppt_mt40) / 10 (Z) `eleReco` Up/Down. Fork dry-run into a
scratch dir (`QCD_MODE=abcdy both simfit`): `eleReco` has 125 entries, all on
electron MC columns (24 W bins × 5 + the Z peak), 0 on muon columns;
`lepsf group = muSF eleSF eleReco`. Expected in the grand fit: the muons pin
r, so the eleReco pull measures the electron Reco efficiency; in the e-only
fit it is degenerate with r (σ_e gains ±10%, σ_Z,e ±21%).

## 3. Step 2 — skim inputs for the no-m_T fit

`skim/skim.C` (both W channels):
- LHE member twins + lepton-SF twins of the no-m_T templates `h_leppt_W{p,m}_y*[_FB]`
  (the fit templates of the no-m_T fit; LHE Up/Down written by `run_lhe_updown.sh`);
- per-rapidity cubes `h_iso_pt_mt_<lep><chg>_y<i>[_FB]` = (iso axis, pT 2 GeV,
  m_T) with coarse iso/m_T edges on every boundary in use, and their twins
  `h_iso_pt_mttrk_*` with m_T from the **charged-only MET**
  (`skim_common.h ComputeChargedPFMET`: PF types 1/2/3, lepton included, no
  η cut, no vertex association).

W MC files 6 → ~24 MB. Every pre-existing histogram bit-identical (§2).

## 4. The multijet study (`correction/qcd_nomt_study.C`, log `correction/logs/qcd_nomt_study.log`)

SR = iso-pass, pT 24–100 (the fit window), all m_T; transfer factor F(pT) =
QCD(iso-pass)/QCD(anti-iso) at m_T < 30 (profile likelihood, EWK at r = 1);
three in-fit schemes compared by Asimov Fisher matrices and data stand-in fits
(stat only, no μ/e sharing, Z peak not included):
**A** = all-m_T SR, level T_g free per group (no m_T anywhere);
**C** = the same SR split at m_T = 30 into two channels sharing r (the low-m_T
channel is region B); **R** = today's leppt_mt40 in-fit ABCD per group.

**Composition (PF MET).** Multijet fraction of the SR: μ⁺/μ⁻ 15% / 17%
(m_T > 40: 5%), e⁺/e⁻ 38% / 47% (m_T > 40: 19% / 27%); e⁺ y10 73%. Half of
it (μ 57–58%, e 45–47%) sits at m_T < 30, i.e. in region B itself.

**F(pT) (PF MET, m_T < 30), F = F0 exp(s (pT − 40)/10).** μ flat
(s = +0.08 / −0.10, flatness p 0.48 / 0.42); e rises: ×1.21 (e⁺, 2.4σ) /
×1.34 (e⁻, 4.1σ) per 10 GeV.
The **level** depends strongly on the group: e F0_fwd,bwd / F0_cen = 2.2–2.9
(p ~ 1e-17), μ⁺ 1.8 (p 5e-4), μ⁻ 1.1–1.5 (p 0.3). The **slope** is common to
the groups for e⁺ (p 0.92) and μ⁻, not for e⁻ (p 0.017: fwd 0.60, cen 0.08)
and μ⁺ (p 0.018). F at m_T 30–40 over m_T < 30 (pT-integrated): e 1.15 / 1.10,
μ 1.04 / 1.04.

**Per-bin anti-iso shapes.** vs the rest of the acceptance 5–6 of 12 bins per
channel at p < 0.05 (the edges); vs the rest of their own group 1–5 of 12 (μ
1–2, e 4–5, mostly 0.002–0.05).

**Expected precision on ⟨r⟩_W (Asimov, stat + the slope prior).**

| | A | C | R |
|---|---|---|---|
| μ⁺ / μ⁻ | 3.6 / 3.9% | 2.4 / 2.8% | 2.4 / 2.7% |
| e⁺ / e⁻ | 4.6 / 5.1% | 3.4 / 4.0% | 3.3 / 3.8% |

Scheme A: correlation r–T_g −0.3…−0.6.

**Data stand-in fits, SR multijet postfit / prefit (the B/D prediction), PF MET.**

| | A (per-y shape) | C (per-y shape) | C (group shape) |
|---|---|---|---|
| μ⁺ fwd/bwd/cen | 0.87 / 0.88 / 0.65 | 0.93 / 0.86 / 0.88 | 0.96 / 0.93 / 0.94 |
| μ⁻ | 0.71 / 1.31 / 0.69 | 0.89 / 1.02 / 0.86 | 0.95 / 1.07 / 0.94 |
| e⁺ | 0.74 / 0.53 / 0.76 (−1.4…−2.9σ) | 0.80 / 0.74 / 0.81 | 0.85 / 0.78 / 0.96 (−0.3…−1.5σ) |
| e⁻ | 0.54 / 0.62 / 0.65 (−2.2…−3.4σ) | 0.72 / 0.74 / 0.72 | 0.79 / 0.76 / 0.84 (−1.2…−2.0σ) |
| χ²_BC/ndf e⁺ / e⁻ | 496/428, 522/428 | 873/872, 1003/872 | 812/872, 944/872 |

R gives 0.70–0.88 for electrons (cf. the real abcdy fit's 0.76–0.84). An
independent slope for the m_T ≥ 30 channel (C_2s) changes nothing (fitted
slopes lo / hi similar), so the deficit is a level, not a pT-factorization
failure. In A the electron edge-bin r swing with the shape choice (e⁻ y11
2.53 / 2.41 / 2.19 per-y / group / inclusive shape).

**e/μ of r (stand-in).** C (group shape): χ²(r_μ = r_e) 6.9/12 (W⁺),
6.5/12 (W⁻); W⁻ edges y0 1.14 ± 0.22, y11 1.09 ± 0.21. R (m_T > 40):
W⁻ edges 0.90 ± 0.18, 0.98 ± 0.18. **Without the m_T cut the electron W⁻ edge
deficit is gone** — consistent with the 2026-10-20 edge cutflow (deficit
present only with m_T > 40 / MET ≥ 30).

**Charged-only MET (user question).** Worse on every measure: the multijet
moves to higher m_T (SR multijet at m_T ≥ 30: 71–74% vs 42–55%), so region B
is smaller and less useful; the electron F(pT) becomes much steeper (×1.5–1.6
per 10 GeV, 6–7σ); the muon F(30–40)/F(<30) degrades (1.33 / 1.46 vs 1.04);
expected ⟨r⟩ precision worse (C: μ⁺ 2.8 vs 2.4%, e⁺ 4.1 vs 3.4%); stand-in
closure worse (C e: 0.46–0.68 of prefit, χ² e⁺ 1054/872 vs 873/872). e⁺ vs
e⁻ (the [CHG] table): with PF MET T⁺/T⁻ is consistent with 1 in every group
(0.83–1.08, ≤ 1.4σ), the slopes differ only in fwd (2.3σ); with charged-only
MET the charges differ more (T⁺/T⁻ fwd 0.76, cen 1.42, ~1.9σ each). In the
no-m_T fit m_T enters only through the multijet estimate (and the m_T > 40
purity selection of the trigger-SF measurement, whose plain/mt40 SFs agree:
e 0.9974 / 1.0010, μ 0.9962 / 0.9971). **Recommendation: keep PF MET.**

## 5. Proposal (for the user to confirm)

Scheme C with the group anti-iso shape:
- SR per (flavour, charge, rapidity bin) = iso-pass, pT 24–100, **no m_T cut**,
  split into two pT channels at m_T = 30 sharing r (all events used once; the
  sum is the no-m_T pT fit; fiducial = the gen set without m_T,
  `gen_xsec_fid.txt`).
- Multijet per bin = T_g (free per group) × exp(s (pT − 40)/10) × the group's
  anti-iso pT shape (per m_T channel) × the bin's anti-iso count (= the
  per-rapidity-bin normalization); the low-m_T channel IS region B and fixes
  T_g in-fit; anchors D_g / C_hi,g as 1-bin CRs; slope per charge (electron:
  free or constrained; muon: flat), the per-bin-vs-group shape difference as a
  shape systematic.
- MET enters only through how the W MC splits at m_T = 30: 2.0% (e) / 2.2% (μ)
  of the SR W sits below 30 → a MET mis-modelling moves ≲ 1% of W
  (estimate ~1% on r via the B-level), to be covered by a nuisance on the W
  lo/hi split — vs 9–15% on the m_T > 40 acceptance of endcap electrons.
- Open: the electron multijet level lands 15–25% below the B/D prediction
  (1–2σ per group) in every scheme including today's; with the no-m_T SR it
  has twice the leverage on r. To decide whether it stays a fitted effect or
  gets a dedicated study (the electron pT-shape mis-modelling).

## 6. Not done

- No fit (lxplus). The no-m_T fit is not wired: needs the `leppt` variant
  back in `mtandmet.C` (its CRs/templates per §5), the two-m_T-channel
  templates with LHE/SF twins (another skim addition), the fork's card
  generator, `disc_variants.h` / `run_observables.sh`, and the AN.
- `run_syst_shapes.sh` not re-run (diagnostic).

## 7. Backups (tag `pre_noreco`, the state before this change)

`backups/pre_noreco/skim/{rootfile,logs,output}/`,
`backups/pre_noreco/correction/{rootfile,logs}/`,
`backups/pre_noreco/correction/plots/{qcd_abcd_mu,qcd_abcd_ele}/`,
`backups/pre_noreco/plotting/{plots,logs}/`,
`backups/pre_noreco/analysis/logs/`. The fork's fit trees were not touched.
