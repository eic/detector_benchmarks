#!/usr/bin/env python3

"""
Collect the json files from individual benchmark tests into
a larger json file that combines all benchmark information,
and do additional accounting for the benchmark.

Tests results are expected to have the following file name and directory
structure:
   results/<BENCHMARK_NAME>/**/<SOME_NAME>.json
where ** implies we check recursively check all sub-directories of <BENCHMARK_NAME>

Internally, we will look for the "tests" keyword in each of these
files to identify them as benchmark components.
"""

## Our benchmark definition file, stored in the benchmark root directory
BENCHMARK_FILE=r'benchmarks/{}/benchmark.json'

## Our benchmark results directory
RESULTS_PATH=r'results/{}'

## Output json file with benchmark results
OUTPUT_FILE=r'results/{}.json'

import argparse
import json
from pathlib import Path

## Exceptions for this module
class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass
class FileNotFoundError(Exception):
    '''File does not exist.

    Attributes:
        file: the file name
        message: error message
    '''
    def __init__(self, file):
        self.file = file
        self.message = 'No such file or directory: {}'.format(file)

class InvalidBenchmarkDefinitionError(Exception):
    '''Raised for missing keys in the benchmark definition.

    Attributes:
        key: the missing key
        file: the benchmark definition file
        message: error message
    '''
    def __init__(self, key, file):
        self.key = key
        self.file = file
        self.message = "key '{}' not found in benchmark file '{}'".format(key, file)

class InvalidTestDefinitionError(Exception):
    '''Raised for missing keys in the test result.

    Attributes:
        key: the missing key
        file: the test result file
        message: error message
    '''
    def __init__(self, key, file):
        self.key = key
        self.file = file
        self.message = "key '{}' not found in test file '{}'".format(key, file)
class InvalidTestResultError(Exception):
    '''Raised for invalid test result value.

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
        self.message = "value '{}' for key '{}' invalid in test file '{}'".format(
                value, key, file)
    
    
parser = argparse.ArgumentParser()
parser.add_argument(
        'benchmark',
        action='append',
        help='One or more benchmarks for which to collect test results.')

def collect_results(benchmark):
    '''Collect benchmark tests and write results to file.'''
    print("Collecting results for benchmark '{}'".format(benchmark))

    ## load the test definition for this benchmark
    results = _load_benchmark(benchmark)

    ## collect the test results
    results['tests'] = _load_tests(benchmark)
    
    ## calculate aggregate test statistics
    results = _aggregate_results(results)

    ## save results to output file
    _save(benchmark, results)

    ## Summarize results
    _print_summary(results)

def _load_benchmark(benchmark):
    '''Load benchmark definition.'''
    benchfile = Path(BENCHMARK_FILE.format(benchmark))
    if not benchfile.exists():
        raise FileNotFoundError(benchfile)
    print('  --> Loading benchmark definition from:', benchfile)
    results = None
    with benchfile.open() as f:
        results = json.load(f)
    ## ensure this is a valid benchmark file
    for key in ('name', 'title', 'description', 'target'):
        if not key in results:
            raise InvalidBenchmarkDefinitionError('target', benchfile)
    return results

def _load_tests(benchmark):
    '''Loop over all test results in benchmark folder and return results.'''
    print('  --> Collecting all test results')
    rootdir = Path(RESULTS_PATH.format(benchmark))
    results = []
    for file in rootdir.glob('**/*.json'):
        print('    --> Loading file:', file, '... ', end='')
        with open(file) as f:
            new_results = json.load(f)
            ## skip files that don't include test results
            if not 'tests' in new_results:
                print('not a test result')
                continue
            ## check if these are valid test results,
            ## raise exception otherwise
            for test in new_results['tests']:
                for key in ('name', 'title', 'description', 'quantity', 'target',
                        'value', 'result'):
                    if not key in test:
                        raise InvalidTestDefinitionError(key, file)
                if test['result'] not in ('pass', 'fail', 'error'):
                    raise InvalidTestResultError('result', test['result'], file)
                ## ensure 'weight' key present, defaulting to 1 in needed
                if not 'weight' in test:
                    test['weight'] = 1.
                ## Append to our test results
                results.append(test)
            print('done')
    return results

def _aggregate_results(results):
    '''Aggregate test results for our benchmark.'''
    print('  --> Aggregating benchmark statistics')
    results['target'] = float(results['target'])
    results['n_tests'] = len(results['tests'])
    results['n_pass'] = len([1 for t in results['tests'] if t['result'] == 'pass'])
    results['n_fail'] = len([1 for t in results['tests'] if t['result'] == 'fail'])
    results['n_error'] = len([1 for t in results['tests'] if t['result'] == 'error'])
    results['maximum'] = sum([t['weight'] for t in results['tests']])
    results['sum'] = sum([t['weight'] for t in results['tests'] if t['result'] == 'pass'])
    if (results['n_tests'] > 0):
        results['value'] = results['sum'] / results['maximum']
        if results['n_error'] > 0:
            results['result'] = 'error'
        elif results['value'] >= results['target']:
            results['result'] = 'pass'
        else:
            results['result'] = 'fail'
    else:
        results['value'] = -1
        results['result'] = 'error'
    return results

def _save(benchmark, results):
    '''Save benchmark results'''
    ofile = Path(OUTPUT_FILE.format(benchmark))
    print('  --> Saving benchmark results to:', ofile)
    with ofile.open('w') as f:
        json.dump(results, f, indent=4)

def _print_summary(results):
    '''Print benchmark summary to the terminal.'''
    print('====================================================================')
    print('Summary for:', results['title'])
    print('Pass: {}, Fail: {}, Error: {} out of {} total tests'.format(
        results['n_pass'], results['n_fail'], results['n_error'],
        results['n_tests']))
    print('Weighted sum: {} / {}'.format(results['sum'], results['maximum']))
    print('Benchmark value: {} (target: {})'.format(
        results['value'], results['target']))
    print('===> status:', results['result'])
    print('====================================================================')


if __name__ == "__main__":
    args = parser.parse_args()
    for benchmark in args.benchmark:
        collect_results(benchmark)
