#!/usr/bin/env python3

""" 
This script will run a CI generator or processing script for multiple configurations.

Author: Sylvester Joosten <sjoosten@anl.gov>
"""

import os
import argparse
from multiprocessing import Pool, get_context
from tempfile import NamedTemporaryFile

class InvalidArgumentError(Exception):
    pass

parser = argparse.ArgumentParser()
parser.add_argument(
        'command',
        help="Script to be launched in parallel")
parser.add_argument(
        '--energy', '-e',
        dest='energies',
        action='append',
        help='One or more beam energy pairs (e.g. 10x100)',
        required=True)
parser.add_argument(
        '--config', '-c',
        dest='configs',
        action='append',
        help='One or more configurations',
        required=True)
parser.add_argument(
        '--leading',
        dest='leads',
        action='append',
        help='One or more leading particles(opt.)',
        required=False)
parser.add_argument(
        '--decay',
        dest='decays',
        action='append',
        help='One or more decay channels (opt.)',
        required=False)
parser.add_argument(
        '--nproc',
        dest='nproc',
        default=5,
        type=int,
        help='Number of processes to launch in parallel',
        required=False)

def worker(command):
    '''Execute the command in a system call, with the supplied argument string.'''
    ## use a temporary file to capture the terminal output, and then
    ## print the terminal output once the command finishes
    with NamedTemporaryFile() as f:
        cmd = [command, ' 2>&1 >', f.name]
        cmd = ' '.join(cmd)
        print("Executing '{}'".format(cmd))
        ret = os.system(cmd)
        with open(f.name) as log:
            print(log.read())
        return ret

if __name__ == '__main__':
    args = parser.parse_args()
    print('Launching CI script in parallel for multiple settings')
    for e in args.energies:
        beam_setting = e.split('x')
        if not beam_setting[0].isnumeric() or not beam_setting[1].isnumeric():
            print("Error: invalid beam energy setting:", e)
            raise InvalidArgumentError

    if not os.path.exists(args.command):
        print("Error: Script not found:", args.command)
        raise InvalidArgumentError

    if args.nproc < 1 or args.nproc > 50:
        print("Error: Invalid process limit (should be 1-50):", args.nproc)
        raise InvalidArgumentError

    print(' - command: {}'.format(args.command))
    print(' - energies: {}'.format(args.energies))
    print(' - config: {}'.format(args.configs))
    print(' - nproc: {}'.format(args.nproc))
    if (args.leads):
        print(' - leading: {}'.format(args.leads))
    if (args.decays):
        print(' - decay: {}'.format(args.decays))

    ## Expand our command and argument list for all combinatorics
    cmds = []
    decays = args.decays if args.decays else [None]
    leads = args.leads if args.leads else [None]
    for e in args.energies:
        for c in args.configs:
            for l in leads:
                for d in decays:
                    beam_setting = e.split('x')
                    cmd = [args.command,
                           '--ebeam', beam_setting[0],
                           '--pbeam', beam_setting[1],
                           '--config', c]
                    if l is not None:
                        cmd += ['--leading', l]
                    if d is not None:
                        cmd += ['--decay', d]
                    cmds.append(' '.join(cmd))

    ## create a process pool
    ## note that I'm using themultiprocessing.get_context function to setup
    ## a context where subprocesses are created using the new "spawn" process
    ## which avoids deadlocks that sometimes happen in the default dispatch
    with get_context('spawn').Pool(processes=args.nproc) as pool:
        return_values = pool.map(worker, cmds)
        ## check if we all exited nicely, else exit with status 1
        if not all(ret == 0 for ret in return_values):
            n_fail = sum([1 for ret in return_values if ret != 0])
            print('ERROR, {} of {} jobs failed'.format(n_fail, len(cmds)))
            print('Return values:', [ret for ret in return_values if ret != 0])
            exit(1)

    ## That's all!
