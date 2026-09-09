#!/bin/bash
# skim/lhe_env.sh -- environment for the LHE-weight python step (lhe_updown.py).
#
# LHAPDF is built from source under ~/local/lhapdf (6.5.6) with python bindings
# for the SAME interpreter that carries PyROOT on this machine (Homebrew
# python 3.12 -- `root-config --python-version`), so one script can
# `import ROOT` and `import lhapdf` together. The brew tap formula binds
# python 3.10 and cannot do that. Build recipe (once):
#   cd ~/local/src && curl -L -o LHAPDF-6.5.6.tar.gz "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.6.tar.gz"
#   tar xf LHAPDF-6.5.6.tar.gz && cd LHAPDF-6.5.6
#   ./configure --prefix=$HOME/local/lhapdf PYTHON=/opt/homebrew/bin/python3.12 && make -j8 && make install
#   # the set (only its .info is read): the shipped pdfsets.index predates it,
#   # so fetch the tarball directly:
#   cd ~/local/lhapdf/share/LHAPDF && curl -sLO https://lhapdfsets.web.cern.ch/current/EPPS21nlo_CT18Anlo_O16.tar.gz \
#     && tar xzf EPPS21nlo_CT18Anlo_O16.tar.gz && rm EPPS21nlo_CT18Anlo_O16.tar.gz
# Usage: `source lhe_env.sh` (bash 3.2 safe). Override the prefix with LHAPDF_PREFIX.
LHAPDF_PREFIX="${LHAPDF_PREFIX:-$HOME/local/lhapdf}"
export PATH="$LHAPDF_PREFIX/bin:$PATH"
export PYTHONPATH="$LHAPDF_PREFIX/lib/python3.12/site-packages${PYTHONPATH:+:$PYTHONPATH}"
export LHAPDF_DATA_PATH="$LHAPDF_PREFIX/share/LHAPDF"
export LHE_PYTHON="${LHE_PYTHON:-/opt/homebrew/bin/python3.12}"
