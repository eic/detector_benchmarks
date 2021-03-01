#!/bin/bash

## =============================================================================
## Generic utility script to parse command line arguments for the various
## bash scripts that control the CI. This script should be source'd with
## command line arguments from a bash-like (non-POSIX) shell such as
## bash or zsh.
##
## To control some of the functionality of the script, you can set the following
## environment variables prior to calling the script:
##   - REQUIRE_DECAY:     require the --decay flag to be set
## =============================================================================

## Commented out because this should be taken care of by the 
## calling script to not enforce a fixed directory structure.
## make sure we launch this script from the project root directory
#PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/..
#pushd ${PROJECT_ROOT}

## =============================================================================
## Step 1: Process the command line arguments

function print_the_help {
  echo "USAGE:    --ebeam E --pbeam E --config C1 --decay D2"
  echo "          [--config C2 --decay D2 --decay D3 ...]"
  echo "REQUIRED ARGUMENTS:"
  echo "          --ebeam       Electron beam energy"
  echo "          --pbeam       Ion beam energy"
  echo "          --config      Generator configuration identifiers (at least one)"
  if [ ! -z ${REQUIRE_DECAY} ]; then
    echo "        --decay       Specific decay particle (e.g. muon)."
  fi
  if [ ! -z ${REQUIRE_LEADING} ]; then
    echo "        --leading     Leading particle of interest (e.g. jpsi)."
  fi
  echo "          -h,--help     Print this message"
  echo ""
  echo "  Generate multiple monte carlo samples for a desired process." 
  exit
}

## Required variables
EBEAM=
PBEAM=
DECAYS=
CONFIG=

while [ $# -gt 0 ]
do
  key="$1"
  case $key in
    --config)
      CONFIG="$2"
      shift # past argument
      shift # past value
      ;;
    --ebeam)
      EBEAM="$2"
      shift # past argument
      shift # past value
      ;;
    --pbeam)
      PBEAM="$2"
      shift # past argument
      shift # past value
      ;;
    --leading)
      LEADING="$2"
      shift # past argument
      shift # past value
      ;;
    --decay)
      DECAY="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      print_the_help
      exit 0
      ;;
    *)    # unknown option
      echo "unknown option"
      exit 1
      ;;
  esac
done

if [ -z $CONFIG ]; then
  echo "ERROR: CONFIG not defined: --config <config>"
  print_the_help
  exit 1
elif [ -z $EBEAM ]; then
  echo "ERROR: EBEAM not defined: --ebeam <energy>"
  print_the_help
  exit 1
elif [ -z $PBEAM ]; then
  echo "ERROR: PBEAM not defined: --pbeam <energy>"
  print_the_help
  exit 1
elif [ -z $LEADING ] && [ ! -z $REQUIRE_LEADING ]; then
  echo "ERROR: LEADING not defined: --leading <channel>"
  print_the_help
  exit 1
elif [ ! -z $LEADING ] && [ -z $REQUIRE_LEADING ]; then
  echo "ERROR: LEADING flag specified but not required"
  print_the_help
  exit 1
elif [ -z $DECAY ] && [ ! -z $REQUIRE_DECAY ]; then
  echo "ERROR: DECAY not defined: --decay <channel>"
  print_the_help
  exit 1
elif [ ! -z $DECAY ] && [ -z $REQUIRE_DECAY ]; then
  echo "ERROR: DECAY flag specified but not required"
  print_the_help
  exit 1
fi

## Export the configured variables
export CONFIG
export EBEAM
export PBEAM
if [ ! -z $REQUIRE_LEADING ]; then
  export LEADING
fi
if [ ! -z $REQUIRE_DECAY ]; then
  export DECAY
fi
