#!/bin/bash

PARTICLE="neutron"

function print_the_help {
  echo "USAGE: ${0} [--particle=$PARTICLE] "
  echo "  OPTIONS: "
  echo "            --particle   particle to generate and simulate"
  exit 
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      shift # past argument
      print_the_help
      ;;
    --particle)
      PARTICLE="$2"
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $1"
      print_the_help
      shift # past argument
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ ! -n  "${JUGGLER_DETECTOR}" ]] ; then 
  export JUGGLER_DETECTOR="topside"
fi

if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=5000
fi

#a place where I can run my events when I want to. - AMJ
#export JUGGLER_N_EVENTS=1000

if [[ ! -n  "${E_start}" ]] ; then
  export E_start=125.0
fi

if [[ ! -n  "${E_end}" ]] ; then
  export E_end=145.0
fi

export JUGGLER_FILE_NAME_TAG="zdc_${PARTICLE}"
export JUGGLER_GEN_FILE="${JUGGLER_FILE_NAME_TAG}.hepmc"

export JUGGLER_SIM_FILE="sim_output/sim_${JUGGLER_FILE_NAME_TAG}.edm4hep.root"
export JUGGLER_REC_FILE="rec_${JUGGLER_FILE_NAME_TAG}.root"

export RESULTS_PATH="results/far_forward/zdc/"
mkdir -p "${RESULTS_PATH}"

echo "JUGGLER_N_EVENTS = ${JUGGLER_N_EVENTS}"
echo "JUGGLER_DETECTOR = ${JUGGLER_DETECTOR}"

# Generate the input events
root -b -q "benchmarks/zdc/scripts/gen_zdc_particles.cxx+(${JUGGLER_N_EVENTS}, \"${PARTICLE}\", ${E_start}, ${E_end}, \"${JUGGLER_FILE_NAME_TAG}.hepmc\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script: generating input events"
  exit 1
fi

# Run geant4 simulations
ddsim --runType batch \
      -v WARNING \
      --part.minimalKineticEnergy 0.5*GeV  \
      --filter.tracker edep0 \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --compactFile ${DETECTOR_PATH}/${JUGGLER_DETECTOR_CONFIG}.xml \
      --inputFiles ${JUGGLER_FILE_NAME_TAG}.hepmc \
      --outputFile ${JUGGLER_SIM_FILE}
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running npdet"
  exit 1
fi

echo "Running analysis root scripts"
root -b -q "benchmarks/zdc/scripts/analysis_zdc_particles.cxx+(\"${JUGGLER_SIM_FILE}\", \"${RESULTS_PATH}\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running root analysis script"
  exit 1
fi

