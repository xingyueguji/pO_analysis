#!/usr/bin/env python3
"""
skim/lhe_updown.py -- theory Up/Down templates from the skim's LHE member twins.

skim.C stores, for every fit-template histogram <hist> of an MC file, three
TH2D twins (x = the fit discriminant, y = member index; filled with the
per-event generator weights w * ttbar_w[i]/ttbar_w[0], member 0 == nominal;
layout in skim/lhe_index.h and in each twin's title). This script turns them
into Combine-style Up/Down TH1Ds, written back INTO the same file (UPDATE,
overwrite; nominal bin errors kept, under/overflow left at the nominal):

  STEP 0    (--norm reco, DEFAULT since 2026-09-14, all three families): every
            member template is first AREA-NORMALIZED to the nominal integral
            over the fit bins 1..nx, i.e. member k is rescaled by I_0/I_k.
            Why: the templates are absolutely normalized with the theory cross
            section, and a member changes that cross section (nPDF +-2.5%
            centrally, +5.7/-8.2% forward; scales +4/-6.6%; alpha_s +-0.5%).
            The measured sigma = r x sigma_gen does not depend on the theory
            cross section (it cancels between r and sigma_gen), so a nuisance
            that moves the template normalization only makes r absorb a shift
            that sigma_gen(nominal) never compensates (the 2026-09-07 fit: the
            +1.3 sigma qcdScale pull moved sigma by -6%). Normalizing to the
            reco nominal removes the whole normalization change and keeps the
            shape; the small acceptance x efficiency part of the change (which
            would need dividing by the GEN-level member integral instead --
            the gen twins of gen_xsec.C) is deliberately deferred (user
            decision 2026-09-14). --norm none = the pre-09-14 raw variation.
  nPDF      <hist>_epps21 (107 members = the LHAPDF set EPPS21nlo_CT18Anlo_O16,
            member for member) -> bin by bin, LHAPDF's official
            PDFSet.uncertainty(), exactly the `_nPDF` branch of the pPb
            prepareSystVariation():
                val = [nominal] + [member 1..106]        (107 numbers, normalized)
                err = pset.uncertainty(val)   # asymmetric Hessian over the
                                              # pairs (1,2),(3,4),...; the set
                                              # is 90% CL, LHAPDF rescales to
                                              # 68.27% (default cl)
                Up = max(val[0] + err.errplus, 0), Down = max(val[0] - err.errminus, 0)
            The member templates are the irreducible input (the Hessian is
            nonlinear in the per-member bin sums), no per-event "up weight" exists.
  qcdScale  <hist>_scale (9 members = ttbar_w[0..8], the (muR,muF) grid) ->
            ENVELOPE: bin by bin, Up = max and Down = min over the members,
            member 0 (nominal) included so Up >= nominal >= Down. Default =
            the 6 variation points that move muR and muF in the same direction
            or one at a time (user decision 2026-09-14); --scale-points 8 adds
            the antagonistic corners (2,0.5) and (0.5,2) back (the 2026-09-07
            "all 9 points" choice).
  alphaS    <hist>_alphas (5 members: 0.118 nominal, 0.116, 0.117, 0.119, 0.120)
            -> SYMMETRIZED (user decision 2026-09-14): per bin
                err  = (N(alpha_s=0.119) - N(alpha_s=0.117)) / 2
                Up   = nominal + err,  Down = nominal - err   (floored at 0)
            i.e. +-0.001 = the PDF4LHC21 68% CL shift, one symmetric template
            pair (no one-sided bins by construction). --alphas-mode pick
            restores the 2026-09-07 recipe (the two member templates as they are).

Member 0 of every twin IS the nominal histogram; the script checks that bin by
bin (a wiring guard, printed as maxdev0 -- must be 0). With --norm reco the
per-template log line also prints the normalization factors that were divided
out (min..max of I_k/I_0 over the members used).

Usage (env from lhe_env.sh: PyROOT python + the lhapdf module + the set):
    python3 lhe_updown.py [--systs nPDF,qcdScale,alphaS] [--set EPPS21nlo_CT18Anlo_O16]
                          [--norm reco|none] [--scale-points 6|8] [--alphas-mode symm|pick]
                          [--alphas-up 3 --alphas-down 2] [--dry-run] FILE.root ...
Wrapper with logging: ./run_lhe_updown.sh [channel|all]  (re-run after every re-skim).
"""
import argparse
import sys

import ROOT  # PyROOT (Homebrew python 3.12 on this Mac; see lhe_env.sh)

# twin suffix (skim/lhe_index.h kFamilySuffix) -> nuisance name written as <hist>_<name>Up/Down
FAMILIES = {"_epps21": "nPDF", "_scale": "qcdScale", "_alphas": "alphaS"}

# ttbar_w[0..8] = (muR,muF), inner loop muF (lhe_index.h kScaleLabel)
SCALE_LABELS = ["(1,1)", "(1,2)", "(1,0.5)", "(2,1)", "(2,2)", "(2,0.5)", "(0.5,1)", "(0.5,2)", "(0.5,0.5)"]
# keyed by the number of VARIATION points; member 0 (nominal) is always part of the max/min
SCALE_MEMBERS = {6: [0, 1, 2, 3, 4, 6, 8],  # drops the antagonistic corners (2,0.5) [5] and (0.5,2) [7]
                 8: list(range(9))}         # all 8 variations
ALPHAS_LABELS = ["0.118 (nominal)", "0.116", "0.117", "0.119", "0.120"]  # members 0..4
# NB no ';' in any title note: TH1::SetTitle splits "title;xaxis;yaxis" at semicolons
NORM_NOTE = {"reco": " -- members area-normalized to the nominal integral first", "none": ""}


def member_matrix(h2, nx, members):
    """Bin contents [m][ix-1] of the twin for ix = 1..nx (Combine reads bins 1..nx)."""
    return [[h2.GetBinContent(ix, m + 1) for ix in range(1, nx + 1)] for m in members]


def normalize_members(vals):
    """Rescale every member to the nominal (first row) integral over the fit bins.

    This divides out the member's change of the theory cross section and keeps
    only its shape ("normalize to the reco nominal"). Returns (normalized rows,
    factors) with factors[m] = I_m / I_0 = what was removed; a member with a
    non-positive integral is left untouched (factor 1)."""
    i0 = sum(vals[0])
    out, fac = [], []
    for v in vals:
        im = sum(v)
        r = im / i0 if (i0 > 0.0 and im > 0.0) else 1.0
        fac.append(r)
        out.append([x / r for x in v] if r != 1.0 else list(v))
    return out, fac


def prepare(nominal, h2, members, norm):
    """Read the member rows, check member 0 == nominal, optionally normalize."""
    nx = nominal.GetNbinsX()
    vals = member_matrix(h2, nx, members)
    maxdev0 = max(abs(vals[0][ix - 1] - nominal.GetBinContent(ix)) for ix in range(1, nx + 1))
    fac = [1.0] * len(members)
    if norm == "reco":
        vals, fac = normalize_members(vals)
    return nx, vals, fac, maxdev0


def clone_pair(nominal, syst, note):
    up = nominal.Clone(f"{nominal.GetName()}_{syst}Up")
    dn = nominal.Clone(f"{nominal.GetName()}_{syst}Down")
    up.SetDirectory(ROOT.nullptr)
    dn.SetDirectory(ROOT.nullptr)
    up.SetTitle(f"{nominal.GetTitle()} {syst} up ({note})")
    dn.SetTitle(f"{nominal.GetTitle()} {syst} down ({note})")
    return up, dn


def check_twin(nominal, h2, nmembers):
    if h2.GetNbinsY() != nmembers:
        raise RuntimeError(f"{h2.GetName()}: {h2.GetNbinsY()} members stored, expected {nmembers}")
    if h2.GetNbinsX() != nominal.GetNbinsX():
        raise RuntimeError(f"{h2.GetName()}: x binning differs from {nominal.GetName()}")


def combine_lhapdf(nominal, h2, pset, syst, norm):
    """nPDF: per bin, PDFSet.uncertainty() of the 107-member column."""
    check_twin(nominal, h2, pset.size)
    nx, vals, fac, maxdev0 = prepare(nominal, h2, range(pset.size), norm)
    up, dn = clone_pair(nominal, syst, f"LHAPDF {pset.name} uncertainty, 68.27% CL{NORM_NOTE[norm]}")
    for ix in range(1, nx + 1):  # Combine reads bins 1..nx
        col = [vals[m][ix - 1] for m in range(pset.size)]
        err = pset.uncertainty(col)  # PDFUncertainty: central, errplus, errminus, errsymm
        up.SetBinContent(ix, max(col[0] + err.errplus, 0.0))
        dn.SetBinContent(ix, max(col[0] - err.errminus, 0.0))
    return up, dn, maxdev0, fac


def combine_envelope(nominal, h2, members, syst, note, norm):
    """qcdScale: per bin, Up = max / Down = min over `members` (member 0 included)."""
    check_twin(nominal, h2, 9)
    nx, vals, fac, maxdev0 = prepare(nominal, h2, members, norm)  # members[0] must be 0
    up, dn = clone_pair(nominal, syst, note + NORM_NOTE[norm])
    for ix in range(1, nx + 1):
        col = [v[ix - 1] for v in vals]
        up.SetBinContent(ix, max(max(col), 0.0))
        dn.SetBinContent(ix, max(min(col), 0.0))
    return up, dn, maxdev0, fac


def combine_pick(nominal, h2, m_up, m_dn, syst, note, norm):
    """alphaS (pick mode): Up/Down = the member templates m_up / m_dn themselves."""
    check_twin(nominal, h2, 5)
    nx, vals, fac, maxdev0 = prepare(nominal, h2, range(5), norm)
    up, dn = clone_pair(nominal, syst, note + NORM_NOTE[norm])
    for ix in range(1, nx + 1):
        up.SetBinContent(ix, max(vals[m_up][ix - 1], 0.0))
        dn.SetBinContent(ix, max(vals[m_dn][ix - 1], 0.0))
    return up, dn, maxdev0, [fac[0], fac[m_up], fac[m_dn]]


def combine_symm(nominal, h2, m_up, m_dn, syst, note, norm):
    """alphaS (symm mode): err = (N_up - N_dn)/2 per bin; Up = nominal + err, Down = nominal - err."""
    check_twin(nominal, h2, 5)
    nx, vals, fac, maxdev0 = prepare(nominal, h2, range(5), norm)
    up, dn = clone_pair(nominal, syst, note + NORM_NOTE[norm])
    for ix in range(1, nx + 1):
        e = 0.5 * (vals[m_up][ix - 1] - vals[m_dn][ix - 1])
        up.SetBinContent(ix, max(vals[0][ix - 1] + e, 0.0))
        dn.SetBinContent(ix, max(vals[0][ix - 1] - e, 0.0))
    return up, dn, maxdev0, [fac[0], fac[m_up], fac[m_dn]]


def process_file(path, cfg, dry_run):
    mode = "READ" if dry_run else "UPDATE"
    f = ROOT.TFile.Open(path, mode)
    if not f or f.IsZombie():
        print(f"[ERR] cannot open {path}")
        return 1
    wanted = {suf: nm for suf, nm in FAMILIES.items() if nm in cfg["systs"]}
    twins = []  # (suffix, twin name)
    for k in f.GetListOfKeys():
        if k.GetClassName() != "TH2D":
            continue
        for suf in wanted:
            if k.GetName().endswith(suf):
                twins.append((suf, k.GetName()))
    twins = sorted(set(twins), key=lambda t: (t[1][: -len(t[0])], t[0]))
    if not twins:
        print(f"[WARN] {path}: no member twins found (data file, or skimmed before the twins existed)")
        f.Close()
        return 0
    print(f"[FILE] {path}: {len(twins)} member twins")
    for suf in wanted:
        first = next((nm for s, nm in twins if s == suf), None)
        if first:
            print(f"[LAYOUT] {f.Get(first).GetTitle()}")
    worst, counts = 0.0, {nm: 0 for nm in wanted.values()}
    for suf, nm in twins:
        h2 = f.Get(nm)
        nom = f.Get(nm[: -len(suf)])
        if not nom:
            print(f"[ERR] {path}: nominal {nm[:-len(suf)]} not found for {nm}")
            return 1
        syst = wanted[suf]
        norm = cfg["norm"]
        if suf == "_epps21":
            up, dn, maxdev0, fac = combine_lhapdf(nom, h2, cfg["pset"], syst, norm)
        elif suf == "_scale":
            up, dn, maxdev0, fac = combine_envelope(nom, h2, cfg["scale_members"], syst, cfg["scale_note"], norm)
        elif cfg["alphas_mode"] == "symm":
            up, dn, maxdev0, fac = combine_symm(nom, h2, cfg["alphas_up"], cfg["alphas_dn"], syst, cfg["alphas_note"], norm)
        else:
            up, dn, maxdev0, fac = combine_pick(nom, h2, cfg["alphas_up"], cfg["alphas_dn"], syst, cfg["alphas_note"], norm)
        worst = max(worst, maxdev0)
        counts[syst] += 1
        i0 = nom.Integral(1, nom.GetNbinsX())
        ru = (up.Integral(1, up.GetNbinsX()) / i0 - 1.0) * 100.0 if i0 > 0 else 0.0
        rd = (dn.Integral(1, dn.GetNbinsX()) / i0 - 1.0) * 100.0 if i0 > 0 else 0.0
        # the normalization change of the members that the recipe divided out (or, with --norm none, kept)
        normtxt = f"  I_k/I_0 {min(fac):.4f}..{max(fac):.4f}{' removed' if norm == 'reco' else ' kept'}"
        print(f"  {nom.GetName():34s} {syst:8s} nominal {i0:12.2f}  up {ru:+6.2f}%  down {rd:+6.2f}%  maxdev0 {maxdev0:.3g}{normtxt}")
        if not dry_run:
            f.cd()
            up.Write(up.GetName(), ROOT.TObject.kOverwrite)
            dn.Write(dn.GetName(), ROOT.TObject.kOverwrite)
    if worst != 0.0:
        print(f"[WARN] {path}: member 0 differs from the nominal histogram (max |diff| = {worst:.3g}) -- wiring problem?")
    done = ", ".join(f"{nm}Up/Down ({n})" for nm, n in counts.items())
    print(f"[SUMMARY] {path}: {done}{' (dry run, nothing written)' if dry_run else ' written'};"
          f" max|member0 - nominal| = {worst:.3g}")
    f.Close()
    return 0


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="MC skim files (rootfile/*_hist.root) holding the member twins")
    ap.add_argument("--systs", default="nPDF,qcdScale,alphaS",
                    help="comma list of the systematics to build (nPDF, qcdScale, alphaS)")
    ap.add_argument("--set", default="EPPS21nlo_CT18Anlo_O16", help="LHAPDF set whose .info defines the nPDF combination")
    ap.add_argument("--norm", default="reco", choices=("reco", "none"),
                    help="reco (default): area-normalize every member to the nominal integral before combining"
                         " (removes the theory-cross-section change, keeps the shape); none: raw variation")
    ap.add_argument("--scale-points", type=int, default=6, choices=(6, 8),
                    help="qcdScale envelope over 6 variation points (default; drops the antagonistic (2,0.5)/(0.5,2))"
                         " or all 8; the nominal is always included in the per-bin max/min")
    ap.add_argument("--alphas-mode", default="symm", choices=("symm", "pick"),
                    help="symm (default): Up/Down = nominal +/- (N_up - N_down)/2 per bin; pick: the two member templates as they are")
    ap.add_argument("--alphas-up", type=int, default=3, help="alphaS 'up' member of <hist>_alphas (default 3 = 0.119)")
    ap.add_argument("--alphas-down", type=int, default=2, help="alphaS 'down' member of <hist>_alphas (default 2 = 0.117)")
    ap.add_argument("--dry-run", action="store_true", help="compute and print, write nothing")
    args = ap.parse_args(argv)

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    systs = [s.strip() for s in args.systs.split(",") if s.strip()]
    bad = [s for s in systs if s not in FAMILIES.values()]
    if bad:
        print(f"[ERR] unknown systematic(s) {bad}; known: {list(FAMILIES.values())}")
        return 2
    cfg = {"systs": systs, "pset": None, "norm": args.norm}
    print(f"[RECIPE] norm     : {args.norm} -> "
          + ("every member is area-normalized to the nominal integral over the fit bins before combining"
             " (the theory-cross-section change is divided out, only the shape is kept; the acceptance x efficiency"
             " part is deferred -- would need the gen-level member integrals)" if args.norm == "reco"
             else "raw variation, the members keep their normalization change"))

    if "nPDF" in systs:
        import lhapdf
        lhapdf.setVerbosity(0)
        pset = lhapdf.getPDFSet(args.set)
        if pset.errorType != "hessian":
            print(f"[ERR] {pset.name} is not a hessian set ({pset.errorType}); the member layout assumes pairs")
            return 2
        cfg["pset"] = pset
        print(f"[LHAPDF] {lhapdf.__version__}: set {pset.name}: {pset.size} members, ErrorType {pset.errorType},"
              f" ErrorConfLevel {pset.errorConfLevel}")
        print(f"[RECIPE] nPDF     : <hist>_epps21 ({pset.size} members) -> per bin PDFSet.uncertainty() = asymmetric Hessian"
              f" over consecutive pairs, rescaled to 68.27% CL (default cl); Up/Down = nominal +/- errplus/errminus")
    if "qcdScale" in systs:
        cfg["scale_members"] = SCALE_MEMBERS[args.scale_points]
        cfg["scale_note"] = (f"(muR,muF) envelope over {args.scale_points} variation points + nominal: per-bin max/min over "
                             + " ".join(SCALE_LABELS[m] for m in cfg["scale_members"]))
        print(f"[RECIPE] qcdScale : <hist>_scale (9 members) -> {cfg['scale_note']}")
    if "alphaS" in systs:
        for m in (args.alphas_up, args.alphas_down):
            if not 0 <= m <= 4:
                print(f"[ERR] alphaS member {m} out of range 0..4")
                return 2
        cfg["alphas_up"], cfg["alphas_dn"] = args.alphas_up, args.alphas_down
        cfg["alphas_mode"] = args.alphas_mode
        if args.alphas_mode == "symm":
            cfg["alphas_note"] = (f"symmetrized: err = (N[alpha_s {ALPHAS_LABELS[args.alphas_up]}] - N[alpha_s"
                                  f" {ALPHAS_LABELS[args.alphas_down]}])/2 per bin, Up/Down = nominal +/- err")
        else:
            cfg["alphas_note"] = (f"Up = member {args.alphas_up} (alpha_s {ALPHAS_LABELS[args.alphas_up]}),"
                                  f" Down = member {args.alphas_down} (alpha_s {ALPHAS_LABELS[args.alphas_down]}) templates")
        print(f"[RECIPE] alphaS   : <hist>_alphas (5 members) -> {cfg['alphas_note']}")

    rc = 0
    for path in args.files:
        rc |= process_file(path, cfg, args.dry_run)
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
