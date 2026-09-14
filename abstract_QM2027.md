# QM 2027 abstracts (shortened 2026-09-13; plain LaTeX math, no CMS macros)

Writing rule (user, 2026-09-13): the abstract should poke interest, not report
results. No numerical results (shifts, significances, x-ranges, cross sections);
dataset identifiers (energy, luminosity) are fine. Both texts below follow it.

## pO — Measurement of W and Z boson production in proton-oxygen collisions at $\sqrt{s_{\mathrm{NN}}} = 9.62$~TeV with CMS

(161 words, −47% vs the 09-13 refined version)

CMS W boson measurements in pPb collisions at 5.02 and 8.16~TeV revealed nuclear modifications of the lead parton densities and now anchor global nuclear PDF fits. The 2025 LHC proton-oxygen run extends these clean probes to a light nucleus bridging the pp and pPb systems, whose isoscalar nature also frees the W charge asymmetry from the neutron-excess isospin effects that dominate it in pPb. We present the first measurement of W and Z boson production in pO collisions at $\sqrt{s_{\mathrm{NN}}} = 9.62$~TeV recorded by CMS. Fiducial $\mathrm{W}^{\pm}$ and Z cross sections are extracted from a simultaneous fit of the muon and electron channels with EPPS21-based templates; the W cross sections, charge asymmetry and forward-backward ratios versus lepton pseudorapidity probe the oxygen parton densities from the shadowing to the antishadowing region. Comparisons with free-proton and nuclear PDF predictions give the first constraints on the partonic structure of a light nucleus at LHC energies and a baseline for hard probes in OO collisions.

Cut to reach the length, restorable: the luminosity ("with 46.5~nb$^{-1}$ of
CMS data" in place of "recorded by CMS", +3 words), the explicit x-range
$10^{-3} \lesssim x \lesssim 10^{-1}$ (now "from the shadowing to the
antishadowing region"), the fiducial cuts, the list of the four nPDF sets
(EPPS21, nCTEQ15HQ, nNNPDF3.0, TUJU21 — now "nuclear PDF predictions").

## Z — mass shift as a probe of the magnetic field in PbPb

(189 words, −30% vs the 09-13 refined version)

Heavy ion collisions are expected to generate the strongest magnetic fields in nature, of order $10^{15}$~T and a prerequisite for the chiral magnetic effect, yet their magnitude and especially their time evolution remain poorly constrained. Z bosons offer a clean probe: decaying within $c\tau \approx 0.08$~fm of the initial hard scattering, their muons traverse the field while it is still strong, and the Lorentz force is predicted to lower and broaden the dimuon mass peak. We present a search for this signature in $\mathrm{Z}\to\mu^{+}\mu^{-}$ decays in 1.8~nb$^{-1}$ of PbPb collisions at $\sqrt{s_{\mathrm{NN}}} = 5.02$~TeV recorded by CMS, against a 13~TeV pp reference reweighted to the PbPb Z kinematics. Peak position and width are extracted with three complementary methods, with detector effects controlled using the $\mathrm{J}/\psi$, and the PbPb$-$pp differences are studied versus centrality. All three indicate a downward mass shift and a broadening in PbPb relative to pp, and a comparison with model predictions spanning the field magnitude, lifetime, and decay profile excludes the largest predicted shifts, disfavors both very short- and very long-lived fields, and provides the first constraints on the field's time evolution from Z boson decays.

Cut, restorable: "about 19,000 candidates", "in 2018", "taken under the same
detector conditions", the names of the three methods (window counting,
Breit--Wigner ⊗ double-sided Crystal Ball fit, template fit to reweighted NLO
simulation), and all numbers ($\Delta\langle M\rangle$, $\Delta\sigma$, the
125 models, the $5\sigma$, the 300--400~MeV shifts).

## Where the pO numbers come from / verify before submission

- √s_NN = 9.62 TeV, L = 46.5 nb⁻¹: `skim/mc_norm.h` (`kLumi_invnb`), `plotting/CMS_lumi.h`. Check 46.5 against the official CMS pO luminosity if it goes back in.
- Fiducial W: lepton pT > 25 GeV, |η_lab| < 2.4 (`skim.C`, `gen_xsec.C::kFidPtMin`). Z: 60 < m_ℓℓ < 120 GeV, leading/subleading lepton pT > 15/10 GeV (`skim_Zmm`/`skim_Zee`).
- CM frame: y_CM = y_lab − 0.3466 (= ½ ln 2), proton-going = forward (`analysis/analysis_helpers.h`). Lab |η| < 2.4 ⇒ η_CM ∈ [−2.75, 2.05]; F/B pairs cover |η_CM| < 2.05 (`kYEdgesFB`).
- x range: x_O ≈ (m_W/√s_NN)·exp(∓y_CM) with m_W/√s_NN = 8.4×10⁻³ ⇒ ≈1×10⁻³ (forward) to ≈0.13 (backward) — the "shadowing to antishadowing" statement.
- nPDF sets: EPPS21, nCTEQ15HQ, nNNPDF3.0, TUJU21 = the MCFM R_pO curves in `plotting/pQCDLightIon/` (`plotRpOtheory.C`); MC templates = POWHEG CT18ANLO × EPPS21 (O16).
- pPb context (from memory, verify before citing): W→μν at 5.02 TeV, 34.6 nb⁻¹ (PLB 750 (2015) 565); W→μν at 8.16 TeV, 173.4 nb⁻¹, first observation of nuclear modifications (PLB 800 (2020) 135048); Drell-Yan at 8.16 TeV (JHEP 05 (2021) 182).
- "First measurement": confirm no other LHC experiment shows W/Z in pO before QM27; otherwise "first CMS measurement".
- "Z cross section": the fit currently returns r_Z; drop the Z σ if only the W σ will be quoted.
