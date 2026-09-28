#!/usr/bin/env python3
"""skim/sf/extract_electron_sf.py -- flatten one or more working points of the
EGM electron ID-SF correctionlib file into small provenance-stamped CSV tables.

    python3 extract_electron_sf.py [electron.json.gz] [WP ...]

Defaults: skim/sf/EGM/electron.json.gz (the 2025Prompt file the user added on
2026-09-22), working points wp90iso and wp90noiso. One CSV per working point is
written next to this script:

    electron_sf_2025Prompt_<WP>.csv

Consumer: correction/idiso_sf_plots.C (the ID+iso SF cross-check), which
overlays the EGM fine-binned SF on the SF we measure on W events and averages
it over our coarse bins. Nothing in the nominal skim reads these files.

Layout of the source (correction "Electron-ID-SF", schema v2):
    category year (2025Prompt) -> category ValType -> category WorkingPoint
    -> multibinning [eta = SUPERCLUSTER eta, signed ; pt], flow = error
Every ValType of one working point shares the same (eta, pt) grid; the CSV has
one row per cell and one column per ValType, plus the cell edges.

Three facts about the file that the CSV header records (checked here, per cell):
  * sfup - sf = sqrt(err_stat^2 + err_syst^2), i.e. the TOTAL error, although
    the file's ValType description says "sfup = sf + syst";
  * the crack cells 1.444 < |eta| < 1.566 hold placeholders (sf = 1);
  * sf = effData/effMC per cell only below 300 GeV -- above it EGM quotes one
    merged SF per eta row while effData/effMC stay per cell (irrelevant here:
    the W electrons of this analysis live at 25-100 GeV).
"""
import csv
import gzip
import hashlib
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SRC = os.path.join(HERE, "EGM", "electron.json.gz")
DEFAULT_WPS = ["wp90iso", "wp90noiso"]
CORRECTION = "Electron-ID-SF"
YEAR = "2025Prompt"
VALTYPES = ["sf", "sfup", "sfdown", "effData", "effMC", "err_stat", "err_syst"]
kRatioPtMax = 300.0  # sf == effData/effMC holds per cell only below this (merged SF above)


def categories(node):
    """{key: value} of a correctionlib category node."""
    if node.get("nodetype") != "category":
        raise ValueError(f"expected a category node, got {node.get('nodetype')}")
    return {c["key"]: c["value"] for c in node["content"]}


def main(argv):
    args = argv[1:]
    src = args[0] if args and args[0].endswith((".json", ".gz")) else DEFAULT_SRC
    wps = [a for a in args if a != src] or DEFAULT_WPS

    raw = open(src, "rb").read()
    md5 = hashlib.md5(raw).hexdigest()
    text = gzip.decompress(raw) if src.endswith(".gz") else raw
    d = json.loads(text)
    if d.get("schema_version") != 2:
        print(f"[ERR] schema_version {d.get('schema_version')} != 2")
        return 1
    corr = [c for c in d["corrections"] if c["name"] == CORRECTION]
    if not corr:
        print(f"[ERR] no correction '{CORRECTION}' in {src}")
        return 1
    corr = corr[0]
    years = categories(corr["data"])
    if YEAR not in years:
        print(f"[ERR] year '{YEAR}' not in {sorted(years)}")
        return 1
    vt = categories(years[YEAR])
    rel_src = os.path.relpath(src, os.path.dirname(os.path.dirname(HERE)))  # repo-relative

    rc = 0
    for wp in wps:
        tables = {}
        for v in VALTYPES:
            wpnodes = categories(vt[v])
            if wp not in wpnodes:
                print(f"[ERR] working point '{wp}' not in {sorted(wpnodes)}")
                return 1
            node = wpnodes[wp]
            if node.get("nodetype") != "multibinning" or node["inputs"] != ["eta", "pt"]:
                print(f"[ERR] {wp}/{v}: unexpected layout {node.get('nodetype')} {node.get('inputs')}")
                return 1
            tables[v] = node
        eta_e, pt_e = tables["sf"]["edges"]
        for v in VALTYPES:  # one grid for all ValTypes
            if tables[v]["edges"] != [eta_e, pt_e]:
                print(f"[ERR] {wp}: ValType {v} has a different (eta, pt) grid")
                return 1
        neta, npt = len(eta_e) - 1, len(pt_e) - 1

        def val(v, ie, ip):  # correctionlib multibinning content: first input slowest
            return tables[v]["content"][ie * npt + ip]

        max_tot_dev, max_ratio_dev = 0.0, 0.0
        rows = []
        for ie in range(neta):
            for ip in range(npt):
                r = {v: val(v, ie, ip) for v in VALTYPES}
                crack = str(eta_e[ie]) not in ("-inf", "inf") and str(eta_e[ie + 1]) not in ("-inf", "inf") \
                    and 1.4 < abs(0.5 * (float(eta_e[ie]) + float(eta_e[ie + 1]))) < 1.6
                if not crack:
                    tot = math.hypot(r["err_stat"], r["err_syst"])
                    max_tot_dev = max(max_tot_dev, abs((r["sfup"] - r["sf"]) - tot),
                                      abs((r["sf"] - r["sfdown"]) - tot))
                    if r["effMC"] > 0 and float(pt_e[ip + 1]) <= kRatioPtMax:
                        max_ratio_dev = max(max_ratio_dev, abs(r["effData"] / r["effMC"] - r["sf"]))
                rows.append([wp, eta_e[ie], eta_e[ie + 1], pt_e[ip], pt_e[ip + 1]] + [r[v] for v in VALTYPES])

        out = os.path.join(HERE, f"electron_sf_{YEAR}_{wp}.csv")
        with open(out, "w", newline="") as fo:
            fo.write(f"# EGM electron ID SF, working point {wp}, year {YEAR} -- flattened by skim/sf/extract_electron_sf.py\n")
            fo.write(f"# source: {rel_src} ({len(raw)} bytes, md5 {md5}); correction '{CORRECTION}' version {corr.get('version')}\n")
            fo.write(f"# extracted {time.strftime('%Y-%m-%d')}; eta = SUPERCLUSTER eta (signed), pt = electron pT; edges inclusive-low\n")
            fo.write(f"# sfup/sfdown = sf +- sqrt(err_stat^2 + err_syst^2) (verified per cell: max deviation {max_tot_dev:.2e}),\n")
            fo.write(f"#   although the file's description says 'sf +- syst'; max |effData/effMC - sf| = {max_ratio_dev:.2e} for pT < {kRatioPtMax:.0f}\n")
            fo.write(f"#   (above {kRatioPtMax:.0f} GeV sf is one merged value per eta row, NOT effData/effMC of the cell)\n")
            fo.write("# crack cells (1.444 < |eta| < 1.566) hold EGM placeholders (sf = 1)\n")
            w = csv.writer(fo)
            w.writerow(["wp", "eta_lo", "eta_hi", "pt_lo", "pt_hi"] + VALTYPES)
            for row in rows:
                w.writerow(row)
        print(f"[OK] {out}: {len(rows)} cells ({neta} eta x {npt} pT); "
              f"sfup-sf vs total error: max dev {max_tot_dev:.2e}; "
              f"effData/effMC vs sf (pT < {kRatioPtMax:.0f}): max dev {max_ratio_dev:.2e}")
        if max_tot_dev > 1e-5:
            print(f"[WARN] {wp}: sfup - sf is NOT the quadrature total error everywhere")
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))
