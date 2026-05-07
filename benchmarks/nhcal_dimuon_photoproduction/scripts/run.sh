#!/bin/sh

# ROOT
export ROOTLIB=$(root-config --libdir)

# PYTHIA8
export PYTHIA8=$(pythia8-config --prefix)

# LHAPDF
export LHAPDFSYS=$(lhapdf-config --prefix)

# Paths
export PATH=$PATH:$LHAPDFSYS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LHAPDFSYS/lib:$PYTHIA8/lib:$ROOTLIB

# Pythia data
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc

SCRIPT_DIR="$(dirname "$0")"

INDEX=${1:?Usage: $0 <indeks> [OUTDIR]}
OUTDIR=${2:-"$(pwd)/sim_input"}

SCRIPT_DIR="$(dirname "$0")"

mkdir -p "$OUTDIR"

OUTBASE="${OUTDIR}/DiMuon_ep_18x275GeV.${INDEX}"
LOG="${OUTDIR}/DiMuon_ep_18x275GeV_${INDEX}.log"

"${SCRIPT_DIR}/pythiaDiMuon" "${SCRIPT_DIR}/ep_DiMuon.cmnd" "$OUTBASE" | tee "$LOG"

