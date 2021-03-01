#!/usr/bin/env python3

"""
Combine the json files from the individual benchmark tests into
a final master json file combining all benchmarks.

Benchmark results are expected to be all json files in the results
directory.
"""

## Our master definition file, the benchmark project directory
MASTER_FILE=r'benchmarks/benchmarks.json'

## Our results directory
RESULTS_PATH=r'results'

## Output json file with all benchmark results
OUTPUT_FILE=r'results/summary.json'

import argparse
import json
from pathlib import Path

## Exceptions for this module
class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass
class FileNotFoundError(Error):
    '''File does not exist.

    Attributes:
        file: the file name
        message: error message
    '''
    def __init__(self, file):
        self.file = file
        self.message = 'No such file or directory: {}'.format(file)

class InvalidDefinitionError(Error):
    '''Raised for missing keys in the definitions.

    Attributes:
        key: the missing key
        file: the definition file
        message: error message
    '''
    def __init__(self, key, file):
        self.key = key
        self.file = file
        self.message = "key '{}' not found in '{}'".format(key, file)

class InvalidResultError(Error):
    '''Raised for invalid benchmark result value.

    Attributes:
        key: the missing key
        value: the invalid value
        file: the benchmark definition file
        message: error message
    '''
    def __init__(self, key, value, file):
        self.key = key
        self.value = value
        self.file = file
        self.message = "value '{}' for key '{}' invalid in benchmark file '{}'".format(
                value, key, file)
    
def collect_benchmarks():
    '''Collect all benchmark results and write results to a single file.'''
    print("Collecting all benchmark results")

    ## load the test definition for this benchmark
    results = _load_master()

    ## collect the test results
    results['benchmarks'] = _load_benchmarks()
    
    ## calculate aggregate test statistics
    results = _aggregate_results(results)

    ## save results to output file
    _save(results)

    ## Summarize results
    for bm in results['benchmarks']:
        _print_benchmark(bm)
    _print_summary(results)

def _load_master():
    '''Load master definition.'''
    master_file = Path(MASTER_FILE)
    if not master_file.exists():
        raise FileNotFoundError(master_file)
    print('  --> Loading master definition from:', master_file)
    results = None
    with master_file.open() as f:
        results = json.load(f)
    ## ensure this is a valid benchmark file
    for key in ('name', 'title', 'description'):
        if not key in results:
            raise InvalidDefinitionError('target', master_file)
    return results

def _load_benchmarks():
    '''Load all benchmark results from the results folder.'''
    print('  --> Collecting all benchmarks')
    rootdir = Path(RESULTS_PATH)
    results = []
    for file in rootdir.glob('*.json'):
        print('    --> Loading file:', file, '... ', end='')
        with open(file) as f:
            bm = json.load(f)
            ## skip files that don't include test results
            if not 'tests' in bm:
                print('skipped (does not contain benchmark results).')
                continue
            ## check if these are valid benchmark results,
            ## raise exception otherwise
            for key in ('name', 'title', 'description', 'target', 'n_tests',
                    'n_pass', 'n_fail', 'n_error', 'maximum', 'sum', 'value',
                    'result'):
                if not key in bm:
                    raise InvalidDefinitionError(key, file)
            if bm['result'] not in ('pass', 'fail', 'error'):
                raise InvalidResultError('result', bm['result'], file)
            ## Append to our test results
            results.append(bm)
            print('done')
    return results

def _aggregate_results(results):
    '''Aggregate benchmark results.'''
    print('  --> Aggregating benchmark statistics')
    results['n_benchmarks'] = len(results['benchmarks'])
    results['n_pass'] = len([1 for t in results['benchmarks'] if t['result'] == 'pass'])
    results['n_fail'] = len([1 for t in results['benchmarks'] if t['result'] == 'fail'])
    results['n_error'] = len([1 for t in results['benchmarks'] if t['result'] == 'error'])
    if results['n_error'] > 0:
        results['result'] = 'error'
    elif results['n_fail'] == 0:
        results['result'] = 'pass'
    else:
        results['result'] = 'fail'
    return results

def _save(results):
    '''Save aggregated benchmark results'''
    ofile = Path(OUTPUT_FILE)
    print('  --> Saving results to:', ofile)
    with ofile.open('w') as f:
        json.dump(results, f, indent=4)

def _print_benchmark(bm):
    '''Print benchmark summary to the terminal.'''
    print('====================================================================')
    print('  Summary for:', bm['title'])
    print('  Pass: {}, Fail: {}, Error: {} out of {} total tests'.format(
        bm['n_pass'], bm['n_fail'], bm['n_error'],
        bm['n_tests']))
    print('  Weighted sum: {} / {}'.format(bm['sum'], bm['maximum']))
    print('  kBenchmark value: {} (target: {})'.format(
        bm['value'], bm['target']))
    print('  ===> status:', bm['result'])

def _print_summary(results):
    '''Print master benchmark summary to the terminal.'''
    print('====================================================================')
    print('MASTER BENCHMARK SUMMARY FOR:', results['title'].upper())
    print('Pass: {}, Fail: {}, Error: {} out of {} total benchmarks'.format(
        results['n_pass'], results['n_fail'], results['n_error'],
        results['n_benchmarks']))
    print('===> status:', results['result'])
    print('====================================================================')


if __name__ == "__main__":
    try:
        collect_benchmarks()
    except Error as e:
        print()
        print('ERROR', e.message)
