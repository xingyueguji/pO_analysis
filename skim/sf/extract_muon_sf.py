#!/usr/bin/env python3
"""skim/sf/extract_muon_sf.py -- reduce the Muon POG correctionlib file to the
corrections this analysis uses, so a small (~0.7 MB) provenance-stamped copy
can live in the repo instead of the 14 MB original.

    python3 extract_muon_sf.py <ScaleFactors_Muon_ID_ISO_2025_schemaV2.json> [out.json]

Keeps the schema-v2 layout unchanged (schema_version, description, corrections)
so correctionlib itself can still open the reduced file; only the `corrections`
list is filtered to KEEP below and a provenance paragraph is appended to the
description (source file name, size, md5, extraction date). The C++ reader
skim/muon_sf.h parses exactly this layout (binning eta -> binning pt ->
category scale_factors).

KEEP:
  NUM_TightID_DEN_TrackerMuons   the skim's muIDTight (CutBasedIdTight)      -> ID SF
  NUM_TightPFIso_DEN_TightID     PF relIso(dBeta, R=0.4) < 0.15 | TightID    -> ISO SF (the W cut, and the
                                 Z cut since skim_Zmm was harmonized from 0.20 to 0.15 on 2026-09-14)
"""
import hashlib
import json
import os
import sys
import time

KEEP = ["NUM_TightID_DEN_TrackerMuons", "NUM_TightPFIso_DEN_TightID"]


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    src = argv[1]
    out = argv[2] if len(argv) > 2 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "muon_sf_2025_TightID_PFIso_schemaV2.json")
    raw = open(src, "rb").read()
    md5 = hashlib.md5(raw).hexdigest()
    d = json.loads(raw)
    if d.get("schema_version") != 2:
        print(f"[ERR] schema_version {d.get('schema_version')} != 2")
        return 1
    names = [c["name"] for c in d["corrections"]]
    missing = [k for k in KEEP if k not in names]
    if missing:
        print(f"[ERR] corrections missing from {src}: {missing}")
        return 1
    red = {k: v for k, v in d.items() if k != "corrections"}
    red["corrections"] = [c for c in d["corrections"] if c["name"] in KEEP]
    prov = (f"\n\n[pO_analysis provenance] Reduced copy made by skim/sf/extract_muon_sf.py on "
            f"{time.strftime('%Y-%m-%d')} from {os.path.basename(src)} ({len(raw)} bytes, md5 {md5}); "
            f"kept only {', '.join(KEEP)} out of {len(names)} corrections. Nothing else changed.")
    red["description"] = red.get("description", "") + prov
    with open(out, "w") as fo:
        json.dump(red, fo, indent=1)
        fo.write("\n")
    print(f"[OK] {out}: {os.path.getsize(out)} bytes, {len(red['corrections'])} corrections; source md5 {md5}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
