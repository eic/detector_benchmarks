#!/bin/sh

#starpro
#setup 64b
#setup ROOT 5.34.38
#source bin/thisroot.csh

export ROOTLIB=/opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/root-6.32.08-zlapeub2v5vfuu3m6bqz6igubplfiios/lib/root

# PYTHIA8
export PYTHIA8=/opt/local

# LHAPDF
export LHAPDFSYS=/opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/lhapdf-6.5.4-eu6c7akvcf5zv3otrikiyvijfcx3wzvh

# Uzupełnienie ścieżek
export PATH=$PATH:$LHAPDFSYS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LHAPDFSYS/lib:$PYTHIA8/lib:$ROOTLIB

# Dane Pythii
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc

INDEX=${1:?Użycie: $0 <indeks> [OUTDIR]}
OUTDIR=${2:-"$(pwd)/sim_input"}

SCRIPT_DIR="$(dirname "$0")"

mkdir -p "$OUTDIR"

OUTBASE="${OUTDIR}/DiMuon_ep_18x275GeV.${INDEX}"
LOG="${OUTDIR}/DiMuon_ep_18x275GeV_${INDEX}.log"

"${SCRIPT_DIR}/pythiaDiMuon" "${SCRIPT_DIR}/ep_DiMuon.cmnd" "$OUTBASE" | tee "$LOG"

