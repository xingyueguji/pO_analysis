# Analysis Note outline — W/Z production in pO at √s_NN = 9.62 TeV

Checklist of content for the AN, mapped to where each ingredient already exists
in this repo (details live in the code / CLAUDE.md / README). Structure loosely
mirrors AN-2017/058 (W in pPb 8.16 TeV), already our recoil reference.
Items marked **[GAP]** don't exist yet and will need work before (or while)
the note is finalized — the note should at minimum state them honestly.

## 1. Introduction & motivation
- [ ] Physics case: W/Z in pO; nPDF sensitivity; the observables:
      A_ch(η_CM), R_FB(|η_CM|) (sum, W⁺, W⁻), σ_fid (W⁺/W⁻/W incl.),
      dσ_fid/dη_CM; Z as validation/standard candle
- [ ] Asymmetric collision system: define lab vs CM frame ONCE, incl. the
      rapidity shift Δy = 0.3466 and the sign convention
      (`analysis/analysis_helpers.h::kPORapidityShift`)
- [ ] Four channels (W→μν, W→eν, Z→μμ, Z→ee); τ treated as background
- [ ] nPDF comparison sets: EPPS21, nCTEQ15HQ, nNNPDF3.0, TUJU21nlo
      (MCFM-based R_pO curves; provenance = `plotting/plotRpOtheory.C` +
      `filelist_theory.txt`)

## 2. Data sample, trigger, luminosity
- [ ] 2025 pO run; L_int = 46.5 nb⁻¹ (`skim/mc_norm.h`) — document the lumi
      source, whether it is per-NN, and its uncertainty **[GAP: lumi unc.]**
- [ ] Certification / event quality: which pO filters are applied
      (auto-detected in `skim/skim_common.h`) + golden-JSON equivalent story
- [ ] Triggers: `HLT_PAL3Mu12`, `HLT_PAL3Ele12`; trigger matching to the
      leading lepton; trigger efficiency **[GAP: not measured]**
- [ ] Input datasets / ntuple production: ggHiNtuplizer EventTree, July-29
      production, EOS paths (`merge_rootfile/` filelists), ~2.2M events

## 3. Simulated samples & normalization
- [ ] POWHEG W⁺/W⁻/DY(m>50) per lepton flavour incl. τ (9 samples; no
      low-mass DY)
- [ ] Per-nucleon σ read from generator weights (⟨w⟩ = σ): W⁺ 6.376 /
      W⁻ 5.464 / DY 1.175 nb; A-scaling σ_pO = A·σ_NN, A = 16
- [ ] Absolute normalization k_s = A·σ·L/N_gen; N_gen = Σ gen weight with no
      selection (`skim/count_ngen.C`; table from `skim/output/ngen.txt` →
      appendix), ~0.9% NLO negative weights; data weight ≡ 1
- [ ] **State clearly: MC is UNEMBEDDED (no pO underlying event)** — its
      consequences (MET shape disagreement) and how the analysis copes
      (discriminant choice, systematics) **[GAP: quantify or embed]**

## 4. Lepton selection & working points
- [ ] Muons: Tight ID + Δβ-corrected PF relIso < 0.15 (MuonPOG convention);
      optimization study (`correction/isolation_mu_tight.C`: AUC_MET 0.943,
      Youden-J optimum 0.156 ⇒ 0.15 confirmed)
- [ ] Electrons: `eleMVAIdWP95` + relIso < 0.095; include WHY (CutWP95 = worst
      QCD rejector; ID switch cut electron QCD ~64%; J-optimum exactly 0.095)
- [ ] relIso definition: Δβ PU correction formula, validation against the
      ntuplizer variable, note pO pileup is tiny
- [ ] ROC curves / WP summary figures (`correction/plotsROC[_ele]/`,
      `summary_*.png`)
- [ ] Electron acceptance statement: 2.5 with 1.4442–1.566 transition gap
      (fix the known 2.4 comment FIXME when writing this)

## 5. Event selection
- [ ] W selection: the 8-step cutflow incl. exact DY-veto definition (OS pair,
      ID+iso legs, mll ∈ (80,110); veto-leg pT 15 (μ) / 10 (e) — intentional)
- [ ] Z selection: OS pair, both legs ID+iso, mass window [60,120]
- [ ] Cutflow tables per channel & sample (`skim/output/*.txt`)
- [ ] Rapidity binning: 12 lab bins in ±2.4 (A_ch); F/B edges symmetric
      around Δy so CM bins pair cleanly (`analysis/FBratio.C` edge list)
- [ ] Nominal W selection has NO MET/m_T cut (MET is the fit discriminant);
      the pT>25 && m_T>40 variant selection for the lepton-pT fits

## 6. Data/MC validation & corrections
- [ ] Z→ℓℓ kinematics data/MC (q_T, y, φ, lepton pT/η/φ)
      (`correction/dataMC_kinematics.C`, both flavours)
- [ ] W leading-lepton kinematics data/MC (+ _mt40 variants; state the QCD
      fraction in the plotted sample)
- [ ] Z mass peak = absolute-normalization + lepton-eff cross-check (MC within
      ~4% of data) (`plotting/dileptonpeak.C`)
- [ ] Hadronic recoil u_∥/u_⊥ vs q_T for Zμμ + Zee: checked stable, NO recoil
      correction applied (`correction/recoil_raw[_ele].C`); MET data/MC
      discrepancy attributed to unembedded MC, not recoil
- [ ] **Corrections NOT applied — list explicitly + plan**: lepton eff SFs
      (tag-and-probe on the Z sample) **[GAP]**, momentum/energy scale &
      smearing **[GAP]**, recoil correction (deferred), trigger eff **[GAP]**

## 7. Background estimation
- [ ] QCD (data-driven ABCD): relIso × PF MET plane per charge; region
      definitions (all boundaries on histogram bin EDGES — cite the 0.005
      binning rule), EWK subtraction with absolute MC
      (`correction/qcd_abcd.C`)
- [ ] Transfer factors + iso-pass QCD yields: μ T≈0.18–0.20 (~3% of high-MET
      signal region); e T≈0.32–0.38 (~10–11% after the ID switch)
- [ ] Closure tests (data vs EWK+ABCD QCD overlays); anti-iso pT-shape
      stability check (electron tilt = the systematic knob)
- [ ] QCD template per discriminant: MET/m_T from the plane; lepton-pT =
      anti-iso shape × MET-plane T (why the pT plane is never counted 2×2);
      m_T>40 keeps ~25% (μ) / ~33% (e) of QCD
- [ ] EWK/τ backgrounds from MC at fixed relative norm: DY→ℓℓ, W→τν, DY→ττ
      (+ W under the Z peak). τ selection validation status **[GAP: pending]**
- [ ] (Superseded Rayleigh sideband-fit method → one line or appendix)

## 8. Signal extraction (fit)
- [ ] Inputs: structured absolute templates per fit region
      (`combine_input_W[_leppt[_mt40]].root`, `combine_input_Z.root`;
      data_obs/signal/z/ztau/wtau/qcd)
- [ ] Fit model (two-parameter): POI r = all W-related MC; dy_norm = all
      DY-related, pinned by the simultaneous Z_incl channel in every per-bin
      card; free qcd_norm; standalone Z_incl card flips roles (r = DY)
- [ ] Discriminants: PF MET shape (nominal) + lepton-pT and
      lepton-pT (m_T>40) variants — motivation (recoil/MET-free) and the
      (pT × m_T) cut-pair optimization (`correction/ptmt_scan.C`: pT>25
      always; m_T plateau 40–45; channel-uniform (25,45))
- [ ] qcd_norm free for all three — note the weak constraint in the pT
      variants and the decision (or ABCD-constrain) **[decision to revisit]**
- [ ] Combine technicalities: FitDiagnostics, per-region cards, workspace;
      fitted yield = r × prefit signal (`extract_pO_yields.C` → `h_yield_*`)
- [ ] Postfit plots with pull pads; fit quality (χ²/GoF, status/covQ),
      dy_norm ≈ 1 / w_norm ≈ 1 sanity checks

## 9. Results
- [ ] A_ch(η_CM) per channel + μ+e overlay
- [ ] R_FB(|η_CM|) sum / W⁺ / W⁻ with the 4 nPDF bands; document the
      count-weighted sum-theory construction (abundance-weighted mean, its
      approximation vs the exact identity)
- [ ] σ_fid (W⁺/W⁻/W incl.) μ vs e — state the σ_fid×ε caveat (gap = channel
      efficiency until SFs exist); dσ_fid/dη_CM per channel
- [ ] Z→ℓℓ result (peak, dy_norm interpretation; quote a Z fiducial σ?)
- [ ] Discriminant stability: met vs leppt vs leppt_mt40 side-by-side
      (per-disc folders from `analysis/run_observables.sh all`)

## 10. Systematic uncertainties (currently stat-only — program to define)
- [ ] QCD ABCD: anti-iso window choice, T stability, iso–MET correlation,
      anti-iso pT shape (e)
- [ ] MC normalization: σ_NN theory unc., A-scaling assumption, lumi
- [ ] Unembedded MC / UE-in-MET modeling (MET-fit specific); recoil
- [ ] Lepton ID/iso/trigger efficiencies; momentum/energy scale & resolution
- [ ] τ background modeling
- [ ] Frame/binning: Δy value, F/B pairing, bin migrations
- [ ] Theory/PDF where acceptance enters (only if extrapolating beyond
      fiducial)

## 11. Summary

## Appendices (candidates)
- [ ] Full cutflow tables; N_gen/⟨w⟩ table; ABCD region details & scans;
      per-ID iso ROC gallery; ptmt-scan maps/ROCs; per-bin postfit gallery;
      theory-curve provenance; reproducibility: repo + fork branch/commits and
      the README end-to-end runbook

## Cross-cutting gaps the note will force (schedule these)
- lepton eff / trigger SFs (tag-and-probe), lumi uncertainty number,
  acceptance if any boson-level comparison is wanted, τ validation,
  the systematics program above, possibly an embedded-MC request.
