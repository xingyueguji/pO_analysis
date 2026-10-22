# Edge study: the W⁻ → eν deficit at |η| > 2.0 (2026-10-20)

Scratch studies made while chasing why the electron W⁻ yield is low in the outermost lab bins
(y0, y11) relative to muons. They are not part of the pipeline. The macros include the repo
headers by absolute path (`/Users/zhenghuang/pO_analysis/skim/...`), write their outputs next to
themselves and are run from this folder:

```bash
cd correction/edge_study && DEVELOPER_DIR=/Library/Developer/CommandLineTools root -l -b -q 'edge_cutflow.C+("loop")'
```

## What is here

| study | macro | record | plots |
|---|---|---|---|
| per-charge cutflow at the edges, φ and run | `edge_cutflow.C` (`"loop"`, then `"report"`) | `loop.log`, `report.log`, `edge_cutflow.csv`, `edge_cutflow.root` | `plots/phi_{e,mu}_{edges,central}` |
| per-charge m_T fit with each charge's own anti-iso multijet shape | `fit_mt_charge.C` (reads the production skim) | `fit_mt_charge.log` | — |
| E/p, track pT, PF accounting, recoil along the electron | `ep_charge.C` (`"loop"`, then `"report"`) | `ep_loop.log`, `ep_report.log`, `ep_charge.root` | — |
| PF MET with the calibrated electron instead of the PF electron | `met_corr.C` | `met_corr.log` | — |
| e⁺ / e⁻ MET and m_T distributions, no MET/m_T cut | `dist_charge.C` (`dist_charge_replot()` redraws) | `dist_charge.log`, `dist_charge.root` | `plots/dist/` |
| gen-level e vs μ W samples, MC A·ε per charge | `gen_flav/gen_flav.C` | `gen_flav/gen_flav.log` | — |
| prefit r with pO charge-dependent electron SFs (sandbox re-skim) | `chargesf_sandbox/` (see its README) | `chargesf_sandbox/compare_chargesf.log` | `chargesf_sandbox/plots/` |
| the skim's own cutflow split by charge (sandbox run of the new skim code) | `skim/skim.C` + `skim_common.h` | `cutflow_by_charge/*.txt` | — |
| which energy the PF electron carries (PF / elePt, eleRawPt, SC energy) | `pf_raw_check.C` | `pf_raw_check.log` | — |
| the same for muons (PF / muPt) | `pf_raw_check_mu.C` | `pf_raw_check_mu.log` | — |

All selections replicate the skim's W selection exactly (data: 6182 electron / 6055 muon W
candidates, as in the skim).

## How to read the plots

### `plots/dist/<var>_<region>_<charge>`: per-charge stacks

- Events: the skim's W selection (isolated, ID'd, trigger-matched leading electron) with **no MET
  or m_T cut**, electron in one region (`edges` |η| > 2.0, `mid` 1.566–2.0, `cen` barrel), one
  charge (`em` = e⁻, `ep` = e⁺).
- `<var>`: `met` = PF MET, `metcorr` = MET with the PF electron replaced by `elePt`
  (see `met_corr.C`), `mt` / `mtcorr` = m_T with each MET.
- Black points: data. Stack: W→eν MC (blue) and DY/τ/tt̄ MC (orange), both × 1.064 (= the
  muon r 1.2 × the electron SF 0.887; no SFs in these loops), and the multijet (green) = that
  charge's anti-isolated data minus the electroweak MC there, normalized so that data = model below
  20 GeV (MET) / 30 GeV (m_T). Bottom pad: (data − model)/σ.
- Read the region above ~25 GeV, which is mostly W: the low region agrees by construction.
- Limitation: because the multijet is normalized at low MET, extra events at low MET are absorbed
  as "more multijet". The stack cannot tell "lost" from "moved to low MET".

### `plots/dist/diff_<var>_<region>`: the charge difference (no multijet model needed)

- Points: data(e⁺) − data(e⁻) per bin, error √(N⁺ + N⁻). Blue: the same difference in MC (× 1.064),
  positive because W⁺ > W⁻.
- Multijet fakes are charge-symmetric (checked: anti-iso e⁻/e⁺ = 1.05 ± 0.09 at low MET at the
  edges), so they cancel. If data behaved like MC the points would follow the blue shape.
- Box: the sums below / above 30 GeV (MET) or 40 GeV (m_T), data and MC, and the χ² of points vs blue.
- Control: `diff_met_cen` (barrel), points on the blue, χ² 17/20.
- At the edges: points below the blue at 5–20 GeV and above it at 25–40 GeV, totals equal
  (+84 data vs +80 MC) → e⁻ events at lower MET than in MC, not missing.
- PF MET vs `metcorr`: the MC barely moves (its PF electron is calibrated, −0.2 GeV), the data move
  (+1.5 GeV at the edges), so comparing `diff_met_*` with `diff_metcorr_*` shows the effect of the
  electron calibration on data only.

### `plots/phi_{e,mu}_{edges,central}`

(data/MC) per charge vs the lepton φ at the W selection + m_T > 40, each normalized to its φ
average. Flat = no detector sector involved (e edges χ² 11.2/11).

### `chargesf_sandbox/plots/`

`rpre_<C>`: prefit r = (data − DY − QCD)/W per lab bin for μ, e with the pp EGM SFs (open) and e
with the pO charge-dependent SFs from the Z→ee tag-and-probe split by charge (filled). `emu_<C>`:
r_e / r_μ for both electron SF choices.

## Findings (statistical errors; the edge bins were chosen after seeing the data)

1. Generator and MC detector response are charge-symmetric: the e and μ W samples agree at gen
   level (e/μ 0.97–0.99 flat in y = bare-electron FSR), the MC A·ε e⁻/e⁺ agree per bin.
2. The e⁻ deficit at |η| > 2.0 appears only once MET or m_T enters (cutflow, per-charge cutflow);
   trigger matching is charge-symmetric, φ and run are flat.
3. No charge dependence in E/p, track pT / elePt, PF accounting, PF-electron fraction.
4. (On the TODO list in CLAUDE.md; the analysis keeps the raw PF MET for now — user decision
   2026-10-20.) A charge-independent data-only bias: PF pT (ΔR < 0.1) / elePt = 0.916–0.92 (|η| > 1.566),
   0.98 (barrel) in data vs 1.01 in MC = 1/(the forest's data-only electron correction). The PF MET
   in data carries the uncorrected electron: per electron the matched PF electron / `eleRawPt` =
   1.000 (median) in data and MC (`pf_raw_check.C`), while `elePt/eleRawPt` = 1.094 (endcap) /
   1.021 (barrel) in data and 0.998 in MC. The supercluster energy is unrelated (PF / E_SC,T ≈ 1.13
   in data and MC alike). Correcting it raises the data electron W yields at
   m_T > 40 / MET ≥ 30 by 9–15% at |η| > 1.566. In Z→ee this appears as extra recoil width, not a
   shift (two back-to-back electrons), which is how `correction/recoil_res.C` saw it (k_e = 1.27).
5. The charge part survives the MET correction: at |η| > 2.0 e⁻ W candidates sit at lower MET
   than in MC relative to e⁺ (QCD-free difference: below MET 30 −9 ± 31 vs +24, above +93 ± 19
   vs +56). Mechanism not identified.
