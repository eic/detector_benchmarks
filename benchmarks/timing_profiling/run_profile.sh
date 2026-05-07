#!/bin/bash

# Detector Profiling Benchmark for npsim
# This script runs npsim with PerformanceProfileSteppingAction from the local build

# Default values
if [[ ! -n "${JUGGLER_N_EVENTS}" ]] ; then
  export JUGGLER_N_EVENTS=10
fi

if [[ ! -n "${DETECTOR_CONFIG}" ]] ; then
  export DETECTOR_CONFIG="epic_craterlake"
fi

# Generate a timestamp
timestamp=$(date '+%Y-%m-%d_%H-%M-%S')

# Set paths relative to detector_profiling root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DETECTOR_PROFILING_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
NPSIM_DIR="${DETECTOR_PROFILING_ROOT}/npsim"
SIM_CAMPAIGN_DIR="${DETECTOR_PROFILING_ROOT}/simulation_campaign_hepmc3"

# Output directories following benchmark conventions
mkdir -p sim_output
mkdir -p results/timing_profiling

# Check if PerformanceProfileSteppingAction plugin is built
PLUGIN_LIB="${NPSIM_DIR}/build/lib/libNPDetPlugins.so"
if [ ! -f "${PLUGIN_LIB}" ]; then
  echo "Error: PerformanceProfileSteppingAction plugin not built."
  echo "Expected: ${PLUGIN_LIB}"
  echo ""
  echo "Please build npsim first:"
  echo "  cd ${NPSIM_DIR}"
  echo "  mkdir -p build && cd build"
  echo "  cmake .. && make"
  exit 1
fi

echo "========================================"
echo "Detector Timing Profiling Benchmark"
echo "========================================"
echo "Events:            ${JUGGLER_N_EVENTS}"
echo "Detector:          ${DETECTOR_CONFIG}"
echo "Timestamp:         ${timestamp}"
echo "Plugin library:    ${PLUGIN_LIB}"
echo "========================================"
echo ""

# Run using the modified simulation_campaign_hepmc3 scripts
cd "${SIM_CAMPAIGN_DIR}"

# Set environment variables for profiling run
export EBEAM=${EBEAM:-18}
export PBEAM=${PBEAM:-275}
export DETECTOR_VERSION=${DETECTOR_VERSION:-"main"}
export TAG_PREFIX="PROFILE/${timestamp}"
export TAG_SUFFIX="PROFILE_${timestamp}"
export COPYRECO=false
export COPYLOG=false
export USERUCIO=false
export PERFORMANCE_PROFILE=true

echo "Running npsim with PerformanceProfileSteppingAction..."
echo ""

# Run the simulation using the modified scripts/run.sh
scripts/run.sh \
  EVGEN/DIS/NC/18x275/minQ2=1/pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1 \
  hepmc3.tree.root \
  ${JUGGLER_N_EVENTS}

echo ""
echo "========================================"
echo "Profiling Complete!"
echo "========================================"
echo ""

# Find the generated histos.root file
HISTOS_FILE=$(find /tmp -name "PROFILE_${timestamp}_*.histos.root" 2>/dev/null | head -1)

if [ -f "${HISTOS_FILE}" ]; then
  # Copy to benchmark output directory
  OUTPUT_HISTOS="${DETECTOR_PROFILING_ROOT}/detector_benchmarks/sim_output/timing_profile_${timestamp}.histos.root"
  cp "${HISTOS_FILE}" "${OUTPUT_HISTOS}"
  echo "Histogram file copied to: ${OUTPUT_HISTOS}"
  echo ""

  # Also find and copy the simulation output
  SIM_FILE=$(find /tmp -name "PROFILE_${timestamp}_*.edm4hep.root" 2>/dev/null | head -1)
  if [ -f "${SIM_FILE}" ]; then
    OUTPUT_SIM="${DETECTOR_PROFILING_ROOT}/detector_benchmarks/sim_output/timing_profile_${timestamp}.edm4hep.root"
    cp "${SIM_FILE}" "${OUTPUT_SIM}"
    echo "Simulation output copied to: ${OUTPUT_SIM}"
  fi

  # Copy log files
  LOG_FILE=$(find /tmp -name "PROFILE_${timestamp}_*.npsim.log" 2>/dev/null | head -1)
  if [ -f "${LOG_FILE}" ]; then
    OUTPUT_LOG="${DETECTOR_PROFILING_ROOT}/detector_benchmarks/sim_output/timing_profile_${timestamp}.npsim.log"
    cp "${LOG_FILE}" "${OUTPUT_LOG}"
    echo "Log file copied to: ${OUTPUT_LOG}"
  fi
else
  echo "Warning: Could not find histogram file in /tmp"
  exit 1
fi

echo ""
echo "========================================"
