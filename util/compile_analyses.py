#!/usr/bin/env python3

"""
Compile all root analysis scripts under
benchmarks/<BENCHMARK>/analysis/*.cxx

Doing this step here rather than during the main benchmark script has
multiple advantages:
    1. Get feedback on syntax errors early on, without wasting compute resources
    2. Avoid race conditions for large benchmarks run in parallel
    3. Make it easier to properly handle the root build directory, as
       this has to exist prior to our attempt to compile, else all will
       fail (this is probably an old bug in root...)

Analysis scripts are expected to have extension 'cxx' and be located in the analysis
subdirectory
"""

## Our analysis path and file extension for glob
ANALYSIS_PATH=r'benchmarks/{}/analysis'
ANALYSIS_EXT = r'cxx'

import argparse
import os
from pathlib import Path

## Exceptions for this module
class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass

class PathNotFoundError(Exception):
    '''Path does not exist.

    Attributes:
        path: the path name
        message: error message
    '''
    def __init__(self, path):
        self.file = file
        self.message = 'No such directory: {}'.format(file)
class NoAnalysesFoundError(Exception):
    '''Did not find any analysis scripts to complile

    Attributes:
        path: the analysis path
        message: error message
    '''
    def __init__(self, path):
        self.file = file
        self.message = 'No analysis found (extension \'{}\' in path: {}'.format(file,
                ANALYSIS_EXT)

class CompilationError(Exception):
    '''Raised when we failed to compile an analysis script

    Attributes:
        file: analysis file name
        path: analysis path
        message: error message
    '''
    def __init__(self, file):
        self.file = file
        self.message = "Analysis '{}' failed to compile".format(file)

parser = argparse.ArgumentParser()
parser.add_argument(
        'benchmark',
        help='A benchmarks for which to compile the analysis scripts.')

def compile_analyses(benchmark):
    '''Compile all analysis scripts for a benchmark.'''
    print("Compiling all analyis scripts for '{}'".format(benchmark))

    ## Ensure our build directory exists
    _init_build_dir(benchmark)

    ## Get a list of all analysis scripts
    _compile_all(benchmark)

    ## All done!
    print('All analyses for', benchmark, 'compiled successfully')

def _init_build_dir(benchmark):
    '''Initialize our ROOT build directory (if using one).'''
    print(' --> Initializing ROOT build directory ...')
    build_prefix = os.getenv('ROOT_BUILD_DIR')
    if build_prefix is None:
        print('    --> ROOT_BUILD_DIR not set, no action needed.')
        return
    ## deduce the root build directory
    pwd = os.getenv('PWD')
    build_dir = '{}/{}/{}'.format(build_prefix, pwd, ANALYSIS_PATH.format(benchmark))
    print("    --> Ensuring directory '{}' exists".format(build_dir))
    os.system('mkdir -p {}'.format(build_dir))

def _compile_all(benchmark):
    '''Compile all analysis for this benchmark.'''
    print(' --> Compiling analysis scripts')
    anadir = Path(ANALYSIS_PATH.format(benchmark))
    if not anadir.exists():
        raise PathNotFoundError(anadir)
    ana_list = []
    for file in anadir.glob('*.{}'.format(ANALYSIS_EXT)):
        ana_list.append(file)
        print('    --> Compiling:', file, flush=True)
        err = os.system(_compile_cmd(file))
        if err:
            raise CompilationError(file)
    if len(ana_list) == 0:
        raise NoAnalysesFoundError(anadir)

def _compile_cmd(file):
    '''Return a one-line shell command to compile an analysis script.'''
    return r'bash -c "root -q -b -e \".L {}+\""'.format(file)

if __name__ == "__main__":
    args = parser.parse_args()
    compile_analyses(args.benchmark)
