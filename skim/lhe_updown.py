#!/usr/bin/env python3
"""
skim/lhe_updown.py -- theory Up/Down templates from the skim's LHE member twins.

skim.C stores, for every fit-template histogram <hist> of an MC file, three
TH2D twins (x = the fit discriminant, y = member index; filled with the
per-event generator weights w * ttbar_w[i]/ttbar_w[0], member 0 == nominal;
layout in skim/lhe_index.h and in each twin's title). This script turns them
into Combine-style Up/Down TH1Ds, written back INTO the same file (UPDATE,
overwrite; nominal bin errors kept, under/overflow left at the nominal):

  nPDF      <hist>_epps21 (107 members = the LHAPDF set EPPS21nlo_CT18Anlo_O16,
            member for member) -> bin by bin, LHAPDF's official
            PDFSet.uncertainty(), exactly the `_nPDF` branch of the pPb
            prepareSystVariation():
                val = [nominal] + [member 1..106]        (107 numbers)
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
            ALL 9 points (user decision 2026-09-07); --scale-points 7 drops the
            (2,0.5) and (0.5,2) corners (the usual 7-point convention).
  alphaS    <hist>_alphas (5 members: 0.118 nominal, 0.116, 0.117, 0.119, 0.120)
            -> Up = the member-3 template (0.119), Down = the member-2 template
            (0.117), i.e. +-0.001 taken as the templates themselves (user
            decision 2026-09-07; PDF4LHC21 quotes +-0.001 at 68% CL).

Member 0 of every twin IS the nominal histogram; the script checks that bin by
bin (a wiring guard, printed as maxdev0 -- must be 0).

Usage (env from lhe_env.sh: PyROOT python + the lhapdf module + the set):
    python3 lhe_updown.py [--systs nPDF,qcdScale,alphaS] [--set EPPS21nlo_CT18Anlo_O16]
                          [--scale-points 9|7] [--alphas-up 3 --alphas-down 2] [--dry-run] FILE.root ...
Wrapper with logging: ./run_lhe_updown.sh [channel|all]  (re-run after every re-skim).
"""
import argparse
import sys

import ROOT  # PyROOT (Homebrew python 3.12 on this Mac; see lhe_env.sh)

# twin suffix (skim/lhe_index.h kFamilySuffix) -> nuisance name written as <hist>_<name>Up/Down
FAMILIES = {"_epps21": "nPDF", "_scale": "qcdScale", "_alphas": "alphaS"}

# ttbar_w[0..8] = (muR,muF), inner loop muF (lhe_index.h kScaleLabel)
SCALE_LABELS = ["(1,1)", "(1,2)", "(1,0.5)", "(2,1)", "(2,2)", "(2,0.5)", "(0.5,1)", "(0.5,2)", "(0.5,0.5)"]
SCALE_MEMBERS = {7: [0, 1, 2, 3, 4, 6, 8],  # 7-point: no (2,0.5) [5], no (0.5,2) [7]
                 9: list(range(9))}
ALPHAS_LABELS = ["0.118 (nominal)", "0.116", "0.117", "0.119", "0.120"]  # members 0..4


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


def combine_lhapdf(nominal, h2, pset, syst):
    """nPDF: per bin, PDFSet.uncertainty() of the 107-member column."""
    check_twin(nominal, h2, pset.size)
    up, dn = clone_pair(nominal, syst, f"LHAPDF {pset.name} uncertainty, 68.27% CL")
    maxdev0 = 0.0
    for ix in range(1, nominal.GetNbinsX() + 1):  # Combine reads bins 1..nx
        val = [h2.GetBinContent(ix, m + 1) for m in range(pset.size)]
        maxdev0 = max(maxdev0, abs(val[0] - nominal.GetBinContent(ix)))
        err = pset.uncertainty(val)  # PDFUncertainty: central, errplus, errminus, errsymm
        up.SetBinContent(ix, max(val[0] + err.errplus, 0.0))
        dn.SetBinContent(ix, max(val[0] - err.errminus, 0.0))
    return up, dn, maxdev0


def combine_envelope(nominal, h2, members, syst, note):
    """qcdScale: per bin, Up = max / Down = min over `members` (member 0 included)."""
    check_twin(nominal, h2, 9)
    up, dn = clone_pair(nominal, syst, note)
    maxdev0 = 0.0
    for ix in range(1, nominal.GetNbinsX() + 1):
        val = [h2.GetBinContent(ix, m + 1) for m in members]
        maxdev0 = max(maxdev0, abs(h2.GetBinContent(ix, 1) - nominal.GetBinContent(ix)))
        up.SetBinContent(ix, max(max(val), 0.0))
        dn.SetBinContent(ix, max(min(val), 0.0))
    return up, dn, maxdev0


def combine_pick(nominal, h2, m_up, m_dn, syst, note):
    """alphaS: Up/Down = the member templates m_up / m_dn themselves."""
    check_twin(nominal, h2, 5)
    up, dn = clone_pair(nominal, syst, note)
    maxdev0 = 0.0
    for ix in range(1, nominal.GetNbinsX() + 1):
        maxdev0 = max(maxdev0, abs(h2.GetBinContent(ix, 1) - nominal.GetBinContent(ix)))
        up.SetBinContent(ix, max(h2.GetBinContent(ix, m_up + 1), 0.0))
        dn.SetBinContent(ix, max(h2.GetBinContent(ix, m_dn + 1), 0.0))
    return up, dn, maxdev0


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
        if suf == "_epps21":
            up, dn, maxdev0 = combine_lhapdf(nom, h2, cfg["pset"], syst)
        elif suf == "_scale":
            up, dn, maxdev0 = combine_envelope(nom, h2, cfg["scale_members"], syst, cfg["scale_note"])
        else:
            up, dn, maxdev0 = combine_pick(nom, h2, cfg["alphas_up"], cfg["alphas_dn"], syst, cfg["alphas_note"])
        worst = max(worst, maxdev0)
        counts[syst] += 1
        i0 = nom.Integral(1, nom.GetNbinsX())
        ru = (up.Integral(1, up.GetNbinsX()) / i0 - 1.0) * 100.0 if i0 > 0 else 0.0
        rd = (dn.Integral(1, dn.GetNbinsX()) / i0 - 1.0) * 100.0 if i0 > 0 else 0.0
        print(f"  {nom.GetName():34s} {syst:8s} nominal {i0:12.2f}  up {ru:+6.2f}%  down {rd:+6.2f}%  maxdev0 {maxdev0:.3g}")
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
    ap.add_argument("--scale-points", type=int, default=9, choices=(7, 9), help="qcdScale envelope: all 9 points (default) or the 7-point set")
    ap.add_argument("--alphas-up", type=int, default=3, help="alphaS Up = this member of <hist>_alphas (default 3 = 0.119)")
    ap.add_argument("--alphas-down", type=int, default=2, help="alphaS Down = this member of <hist>_alphas (default 2 = 0.117)")
    ap.add_argument("--dry-run", action="store_true", help="compute and print, write nothing")
    args = ap.parse_args(argv)

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    systs = [s.strip() for s in args.systs.split(",") if s.strip()]
    bad = [s for s in systs if s not in FAMILIES.values()]
    if bad:
        print(f"[ERR] unknown systematic(s) {bad}; known: {list(FAMILIES.values())}")
        return 2
    cfg = {"systs": systs, "pset": None}

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
        cfg["scale_note"] = (f"{args.scale_points}-point (muR,muF) envelope: per-bin max/min over "
                             + " ".join(SCALE_LABELS[m] for m in cfg["scale_members"]))
        print(f"[RECIPE] qcdScale : <hist>_scale (9 members) -> {cfg['scale_note']}")
    if "alphaS" in systs:
        for m in (args.alphas_up, args.alphas_down):
            if not 0 <= m <= 4:
                print(f"[ERR] alphaS member {m} out of range 0..4")
                return 2
        cfg["alphas_up"], cfg["alphas_dn"] = args.alphas_up, args.alphas_down
        cfg["alphas_note"] = (f"Up = member {args.alphas_up} (alpha_s {ALPHAS_LABELS[args.alphas_up]}),"
                              f" Down = member {args.alphas_down} (alpha_s {ALPHAS_LABELS[args.alphas_down]}) templates")
        print(f"[RECIPE] alphaS   : <hist>_alphas (5 members) -> {cfg['alphas_note']}")

    rc = 0
    for path in args.files:
        rc |= process_file(path, cfg, args.dry_run)
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
