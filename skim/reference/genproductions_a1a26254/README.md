# Reference copy: CMS genproductions POWHEG scripts (archive a1a26254)

Copied 2026-09-02 from
`/Users/zhenghuang/Downloads/genproductions_scripts-a1a2625485abc436ed10809327323d5d5daefc5e.zip`
(cms-sw/genproductions at commit a1a2625485abc436ed10809327323d5d5daefc5e), the
scripts that produced the pO POWHEG samples (`POWHEG_9p62TeV_2025Run3`, producer
anstahll). Only the files needed to explain the 217 per-event LHE weights stored
in `hiEvtAnalyzer/HiTree::ttbar_w` are kept. **Read-only reference** -- the
decoder `skim/lhe_weights.C` points at these files by `file:line`, so do not edit
them (re-copy from a newer archive instead and update the pointers).

| file | role in the weight list |
|---|---|
| `bin/Powheg/make_rwl.py` | writes `pwg-rwl.dat`, the `<initrwgt>` header POWHEG reads at generation. The `"EPPS21" in Period` branch (L481-540) = the 3x3 scale grid (L43-54) + a "hessian" PDF group + a "replica" PDF group = 110 weights = `ttbar_w[0..109]`. |
| `bin/Powheg/runcmsgrid_powheg.sh` | the gridpack run script. L252-266: if the card has `nPDFerrSet`, rerun `pwhg_main` 107 times with `rwl_add 1`, once per `nPDFerrSet = 1..107`, appending weight ids 9001-9107 (`EPPS21_variation`) = `ttbar_w[110..216]`. |
| `bin/Powheg/run_pwg_condor.py` | the driver. L324-326: for an `EPPS21_*` ion, `defaultPDF = 14600` (CT18ANLO, used as `lhapdf=` in the scale block) and `period = "Run3_" + ion`. L648-667: the card must have `lhans1 == lhans2 == 14600`. |
| `bin/Powheg/Templates/runGetSource_template.sh` | L55 calls `make_rwl.py`; L110-116 apply the EPPS21 patches and fetch the `EPPS21NLOR_16` R-factor grid for the `*_O` period. |
| `bin/Powheg/patches/EPPS21/*.patch`, `EPPS21.f` | the nuclear PDF itself: for the beam with `ia >= 16` the LHAPDF proton PDF is multiplied by the EPPS21 R-factors of error set `nPDFerrSet` (`lhapdf6if_nPDF.patch` L38-40, L80-95) and isospin-averaged (L63-67); `pdfcalls_nPDF.patch` passes `ia1`/`ia2` per beam. |
| `MetaData/npdflist_O_5f_run3.dat` | the O16 nPDF list of the `Run3_O` branch (LHAPDF-grid nPDFs incl. TUJU21/nNNPDF3.0). NOT what this production used -- kept to show the alternative. |

Consequence for the analysis (see `skim/lhe_weights.C` header): the O16 PDF is
`f_O = R_EPPS21(nPDFerrSet) x f_p^LHAPDF`, and the two weight kinds vary the two
factors separately -- `lhapdf=X` weights (idx 9-109) swap `f_p` on both beams
with R fixed at set 1, `nPDFerrSet` weights (idx 110-216) vary R with `f_p` fixed
at CT18ANLO central. Both idx 44 (`lhapdf=14600`) and idx 110 (`nPDFerrSet=1`)
are therefore exactly the nominal.
