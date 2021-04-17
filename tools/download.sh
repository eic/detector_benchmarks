#!/bin/bash

## =============================================================================
## Download generator & reconstruction artifacts for one or more physics
## processes.
## =============================================================================

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/..
pushd ${PROJECT_ROOT}

PROCS=()
BRANCH="master"

function print_the_help {
  echo "USAGE:    -p process [-p process2] [-b git_branch]"
  echo "OPTIONS:"
  echo "          -p,--process  Physics process name (can be defined multiple
  times)."
  echo "          -b,--branch   Git branch to download artifacts from (D:
  $BRANCH)"
  echo "          -h,--help     Print this message"
  echo ""
  echo "  This script will download the relevant generator artifacts needed"
  echo "  for local testing of the benchmarks."
  exit
}

while [ $# -gt 0 ]
do
  key="$1"
  case $key in
    -p|--process)
      PROCS+=("$2")
      shift # past argument
      shift # past value
      ;;
    -b|--branch)
      BRANCH="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      print_the_help
      shift
      ;;
    *)    # unknown option
      echo "unknown option: $1"
      exit 1
      ;;
  esac
done

echo "Downloading generator & reconstruction artifacts for one or more physics processes"

if [ ${#PROCS[@]} -eq 0 ]; then
  echo "ERROR: need one or more processes: -p <process name> "
  exit 1
fi

for proc in ${PROCS[@]}; do
  echo "Dowloading artifacts for $proc (branch: $BRANCH)"
  wget https://eicweb.phy.anl.gov/EIC/benchmarks/physics_benchmarks/-/jobs/artifacts/$BRANCH/download?job=${proc}:generate -O results_gen.zip
  ## FIXME this needs to be smarter, probably through more flags...
  wget https://eicweb.phy.anl.gov/EIC/benchmarks/physics_benchmarks/-/jobs/artifacts/$BRANCH/download?job=${proc}:process -O results_rec.zip
  echo "Unpacking artifacts..."
  unzip -u -o results_gen.zip
  unzip -u -o results_rec.zip
  echo "Cleaning up..."
  rm results_???.zip
done
popd
echo "All done"
