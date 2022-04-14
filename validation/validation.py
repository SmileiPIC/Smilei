#!/usr/bin/env python

"""
This script can do three things:
  (1) generate validation reference(s) for given benchmark(s)
  (2) compare benchmark(s) to their reference(s)
  (3) show visually differences between benchmark(s) and their reference(s)

Usage
#######
Type this command for details: 'python validation.py -h'

Files used in the validation process:
#####################################
The `validation` directory which contains:
  - the `references` directory with one file for each benchmark
  - the `analyses` directory with one validation file for each benchmark
  - a `workdirs` directory, created during the validation process
  - archived `workdirs` for previous versions of smilei

A `workdirs` contains:
 - the smilei binary : `smilei`
 - the compilation output : `compilation_out`
 - a directory `wd_<input_file>/<o>/<m>`, containing the output files
 or
 - the compilation errors file : `compilation_errors`

Steps of the script:
####################
1 - Compilation:
----------------
If the `workdirs` directory lacks a smilei binary (or it is too old),
then the `workdirs` is backed up, and a new compilation occurs.
If compiling errors occur, `compilation_errors` is created and the script exits with status 3.

2 - Execution:
--------------
If `wd_<input_file>/<o>/<m>` does not exist then:
- it is created
- smilei is executed in that directory for the requested benchmark
- if execution fails, the script exits with status 2

3 - Validation:
---------------
Loops through all requested benchmarks
    Runs the benchmark in the `workdir`
    If requested to generate references:
        Executes the `validate_*` script and stores the result as reference data
    If requested to compare to previous references
        Executes the `validate_*` script and compares the result to the reference data
    If requested to show differences to previous references
        Executes the `validate_*` script and plots the result vs. the reference data

Exit status:
============
0  validated
1  validation fails
2  execution fails
3  compilation fails
4  bad option

Remark:
=======
You may define the SMILEI_ROOT environment variable to use a different installation folder
"""

from sys import argv, exit
from glob import glob
from getopt import getopt, GetoptError
from easi import Validation

def usage():
    print( 'Usage: validation.py [-c] [-h] [-v] [-b <bench_case>] [-o <nb_OMPThreads>] [-m <nb_MPIProcs>] [-g | -s] [-r <nb_restarts>] [-t <max_time>] [-k <compile_mode>] [-p <partition name>] [-l <logs_folder>]' )
    print( '    Try `validation.py -h` for more details' )

# Get command-line options
try:
    external_options, remainder = getopt(
        argv[1:],
        'o:m:b:r:k:p:gshvcl:t:a:n:',
        ['OMP=', 'MPI=', 'BENCH=', 'RESTARTS=', 'PARTITION=', 'GENERATE', 'SHOW', 'HELP', 'VERBOSE', 'COMPILE_ONLY', 'COMPILE_MODE=', 'LOG=', 'time=', 'account=', 'nodes=', 'resource-file=']
    )
except GetoptError as err:
    usage()
    exit(4)

# Process options
options = {}
for opt, arg in external_options:
    if opt in ('-o', '--OMP'):
        options['omp'] = int(arg)
    elif opt in ('-m', '--MPI'):
        options['mpi'] = int(arg)
    elif opt in ('-n', '--nodes'):
        options['nodes'] = int(arg)
    elif opt in ('--resource-file',):
        options['resource-file'] = arg
    elif opt in ('-b', '--BENCH'):
        options['bench'] = arg
    elif opt in ('-p', '--PARTITION'):
        options['partition'] = arg
    elif opt in ('-c', '--COMPILE_ONLY'):
        options['compile_only'] = True
    elif opt in ('-k', '--COMPILE_MODE'):
        options['compile_mode'] = arg
    elif opt in ('-t', '--time'):
        options['max_time'] = arg
    elif opt in ('-a', '--account'):
        options['account'] = arg
    elif opt in ('-g', '--GENERATE'):
        options['generate'] = True
    elif opt in ('-s', '--SHOW'):
        options['showdiff'] = True
    elif opt in ('-v', '--VERBOSE'):
        options['verbose'] = True
    elif opt in ('-r', '--RESTARTS'):
        try:
            options['nb_restarts'] = int(arg)
            if options['nb_restarts'] < 0: raise
        except:
            print("Error: the number of restarts (option -r) must be a positive integer")
            exit(4)
    elif opt in ('-l', '--LOG'):
        options['log'] = arg
    elif opt in ('-h', '--HELP'):
        print( """
Options:
  -b <bench_case>
       <bench_case> : benchmark(s) to validate. Accepts wildcards.
       DEFAULT : All benchmarks are validated.
  -o <nb_OMPThreads>
       <nb_OMPThreads> : number of OpenMP threads used for the execution
       DEFAULT : 4
  -m <nb_MPIProcs>
       <nb_MPIProcs> : number of MPI processes used for the execution
       DEFAULT : 4
  -n <nb_nodes>
       <nb_nodes> : number of nodes to use for the execution
       DEFAULT : computed so that 1 mpi = 1 socket
  --resource-file <file>
       <file> : file containing benchmark-specific resources to override values in -o, -m and -n
  -g
       Generates the references
  -s
       Plot differences with references (python -i option required to keep figures on screen)
  -c
       Compilation only
  -k
       Compilation using config=... See make help for details
  -r <nb_restarts>
       <nb_restarts> : number of restarts to run, as long as the simulations provide them.
       DEFAULT : 0 (meaning no restarts, only one simulation)
  -v
       Verbose mode
  -t <max time>
       <max time>: format hh:mm:ss
  -p <partition name>
       <partition name>: partition name on super-computers
                           - ruche
                           - tornado
                           - jollyjumper
                           - irene_skylake
                           - irene_a64fx
  -a <account id>
  --account <account id>
       <account id>: account/project id given by some super-computer facilities
  -l
       Log some performance info in the directory `logs`
""")
        exit(0)


v = Validation( **options )
v.compile()
if 'compile_only' not in options:
    v.run_all()