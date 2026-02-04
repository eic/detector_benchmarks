#!/bin/sh

echo "" >&2
echo "=== DEBUG: Looking for ROOT ===" >&2
find /opt -name "libCore.so*" 2>/dev/null | head -5 || echo "libCore.so not found" >&2
find /opt -name "root-config" 2>/dev/null | head -5 || echo "root-config not found" >&2
which root 2>/dev/null || echo "root not in PATH" >&2

echo "" >&2
echo "=== DEBUG: Looking for PYTHIA8 ===" >&2
find /opt -name "libpythia8*" 2>/dev/null | head -5 || echo "libpythia8 not found" >&2

echo "" >&2
echo "=== DEBUG: Looking for LHAPDF ===" >&2
find /opt -name "liblhapdf*" 2>/dev/null | head -5 || echo "liblhapdf not found" >&2

# ROOT
export ROOTLIB=/opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/root-6.32.08-zlapeub2v5vfuu3m6bqz6igubplfiios/lib/root

# PYTHIA8
export PYTHIA8=/opt/local

# LHAPDF
export LHAPDFSYS=/opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/lhapdf-6.5.4-eu6c7akvcf5zv3otrikiyvijfcx3wzvh

# Paths
export PATH=$PATH:$LHAPDFSYS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LHAPDFSYS/lib:$PYTHIA8/lib:$ROOTLIB

# Pythia data
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc

echo "" >&2
echo "=== DEBUG: Environment variables ===" >&2
echo "ROOTLIB=$ROOTLIB" >&2
echo "PYTHIA8=$PYTHIA8" >&2
echo "LHAPDFSYS=$LHAPDFSYS" >&2
echo "PATH=$PATH" >&2
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >&2

echo "" >&2
echo "=== DEBUG: pythiaDiMuon dependencies ===" >&2
SCRIPT_DIR="$(dirname "$0")"
ldd "${SCRIPT_DIR}/pythiaDiMuon" 2>&1 >&2 || echo "ldd failed" >&2

INDEX=${1:?Usage: $0 <indeks> [OUTDIR]}
OUTDIR=${2:-"$(pwd)/sim_input"}

SCRIPT_DIR="$(dirname "$0")"

mkdir -p "$OUTDIR"

OUTBASE="${OUTDIR}/DiMuon_ep_18x275GeV.${INDEX}"
LOG="${OUTDIR}/DiMuon_ep_18x275GeV_${INDEX}.log"

"${SCRIPT_DIR}/pythiaDiMuon" "${SCRIPT_DIR}/ep_DiMuon.cmnd" "$OUTBASE" | tee "$LOG"

