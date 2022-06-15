#!/bin/bash

function print_the_help {
  echo "USAGE: ${0} -nevents <nevents> -p <particle> -e <energy>"
  echo "  OPTIONS: "
  echo "    -n,--nevents     number of events"
  echo "    -p,--particle    particle type"
  echo "    -e,--energy      particle energy"
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
    -n|--nevents)
      nevents="$2"
      shift # past argument
      shift # past value
      ;;
    -p|--particle)
      particle="$2"
      shift # past argument
      shift # past value
      ;;
    -e|--energy)
      energy="$2"
      shift
      shift
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

if [[ ! -n  "${nevents}" ]] ; then
  nevents="1"
fi

if [[ ! -n  "${particle}" ]] ; then
  particle="e-"
fi

if [[ ! -n  "${energy}" ]] ; then
  energy="1*GeV"
fi

if [[ ! -n "${JUGGLER_DETECTOR}" ]] ; then
  export JUGGLER_DETECTOR="athena"
fi

if [[ ! -n "${DETECTOR_PATH}" ]] ; then
  export DETECTOR_PATH="/opt/detector/share/athena"
fi

if [[ ! -n "${JUGGLER_INSTALL_PREFIX}" ]] ; then
  export JUGGLER_INSTALL_PREFIX="/usr/local"
fi

compact_path=${DETECTOR_PATH}/${JUGGLER_DETECTOR_CONFIG}.xml

echo "DETECTOR_PATH = ${DETECTOR_PATH}"
echo "JUGGLER_DETECTOR = ${JUGGLER_DETECTOR}"

# Run geant4 simulations
output_dir="data/timing/${particle}/${energy/\*/}"
output_file="sim_${nevents}.edm4hep.root"
mkdir -p ${output_dir}
timing_dir="results/timing/${particle}/${energy/\*/}"
timing_file="time_${nevents}events.log"
ddsim_file="npsim_${nevents}events.log"
mkdir -p ${timing_dir}
/usr/bin/time -v -o ${timing_dir}/time_${nevents}events.log \
  ddsim --runType batch \
      --printLevel WARNING \
      --filter.tracker edep0 \
      --numberOfEvents ${nevents} \
      --enableGun \
      --gun.energy "${energy}" \
      --gun.particle "${particle}" \
      --gun.thetaMin "45*deg" \
      --gun.thetaMax "135*deg" \
      --gun.distribution "cos(theta)" \
      --part.minimalKineticEnergy "1*TeV" \
      --compactFile ${compact_path} \
      --outputFile ${output_dir}/${output_file} \
    2>&1 > ${timing_dir}/${ddsim_file}
echo "For ${nevents} events:"
cat ${timing_dir}/${timing_file}

if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running npsim"
  exit 1
fi

