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


# IMPORTS
import sys, os, re, json, pickle, numpy as np
from time import sleep, ctime, strftime
from math import ceil
from glob import glob
from shutil import rmtree, copy2
from getopt import getopt, GetoptError
from inspect import stack
from socket import gethostname
from subprocess import call, check_call, check_output, CalledProcessError
from os import path
s = os.sep

try:
    execfile
except: # python3
    def execfile(file):
        exec(compile(open(file).read(), file, 'exec'), globals())

try:
    raw_input
except: # python3
    raw_input = input


def usage():
    print( 'Usage: validation.py [-c] [-h] [-v] [-b <bench_case>] [-o <nb_OMPThreads>] [-m <nb_MPIProcs>] [-g | -s] [-r <nb_restarts>] [-k <compile_mode>] [-p <partition name>] [-l <logs_folder>]' )
    print( '    Try `validation.py -h` for more details' )

# --------------------
# Paths & variables
# --------------------

# Dictionnary containing all parameters
parameters = {}

if "SMILEI_ROOT" in os.environ :
    parameters['smilei_root'] = os.environ["SMILEI_ROOT"]+s
else:
    parameters['smilei_root'] = path.dirname(path.abspath(stack()[0][1]))+s+".."+s
parameters['smilei_root'] = path.abspath(parameters['smilei_root'])+s
parameters['smilei_benchmark_path'] = parameters['smilei_root']+"benchmarks"+s
parameters['smilei_script_path'] = parameters['smilei_root']+"scripts"+s
parameters['smilei_validation_path'] = parameters['smilei_root']+"validation"+s
parameters['smilei_reference_path'] = parameters['smilei_validation_path']+"references"+s
parameters['smilei_analyse_path'] = parameters['smilei_validation_path']+"analyses"+s

parameters['exec_script'] = 'exec_script.sh'
parameters['exec_script_output'] = 'exec_script.out'
parameters['output_file'] = 'smilei_exe.out'

WORKDIR_BASE = parameters['smilei_root']+"validation"+s+"workdirs"
SMILEI_W = WORKDIR_BASE+s+"smilei"
SMILEI_R = parameters['smilei_root']+s+"smilei"

SMILEI_TOOLS_W = WORKDIR_BASE+s+"smilei_tables"
SMILEI_TOOLS_R = parameters['smilei_root']+s+"smilei_tables"

COMPILE_ERRORS = WORKDIR_BASE+s+'compilation_errors'
COMPILE_OUT = WORKDIR_BASE+s+'compilation_out'
COMPILE_OUT_TMP = WORKDIR_BASE+s+'compilation_out_temp'

POINCARE = "poincare"
LLR = "llrlsi-gw"
HOSTNAME = gethostname()

# Get the current version of Smilei
INITIAL_DIRECTORY = os.getcwd()
os.chdir(parameters['smilei_root'])
gitversion = check_output( "echo `git log -n 1 --format=%h`-", shell=True ).decode()[:-1]
if 'CI_COMMIT_BRANCH' in os.environ:
    gitversion += os.environ['CI_COMMIT_BRANCH']
else:
    gitversion += check_output("echo `git rev-parse --abbrev-ref HEAD`", shell=True ).decode()[:-1]
os.chdir(INITIAL_DIRECTORY)

# Load the happi module
sys.path.insert(0, parameters['smilei_root'])
import happi
# import tools
sys.path.append(parameters['smilei_validation_path'] + '/lib/')
from log import *
# import specific machine modules
from ruche import *
from llr import *
from other import *
from poincare import *
from irene_a64fx import *
from irene_skylake import *

# ------------------------
# Get command-line options
# ------------------------

# General options
options = {}

# Default values
options['omp'] = 12
options['mpi'] = 4
options['ppn'] = 12
options['max_time'] = "00:10:00"
options['execution'] = False
options['verbose'] = False
options['bench']=""
options['compile_only'] = False
options['generate'] = False
options['showdiff'] = False
options['nb_restarts'] = 0
options['compile_mode']=""
options['log'] = False
options['partition'] = "jollyjumper"
options['accout'] = "" # Account for some super-computers

try:
    external_options, remainder = getopt(
        sys.argv[1:],
        'o:m:b:r:k:p:gshvcl:t:a:',
        ['OMP=', 'MPI=', 'BENCH=', 'RESTARTS=', 'PARTITION=', 'GENERATE', 'SHOW', 'HELP', 'VERBOSE', 'COMPILE_ONLY', 'COMPILE_MODE=', 'LOG=', 'time=', 'account=']
    )
except GetoptError as err:
    usage()
    sys.exit(4)

# Process options
for opt, arg in external_options:
    if opt in ('-o', '--OMP'):
        options['execution'] = True
        options['omp'] = int(arg)
    elif opt in ('-m', '--MPI'):
        options['execution'] = True
        options['mpi'] = int(arg)
    elif opt in ('-b', '--BENCH'):
        options['bench'] = arg
    elif opt in ('-p', '--PARTITION'):
        options['partition'] = arg
    elif opt in ('-c', '--COMPILE_ONLY'):
        options['compile_only']=True
    elif opt in ('-k', '--COMPILE_MODE'):
        options['compile_mode']=arg
    elif opt in ('-t', '--time'):
        options['max_time']=arg
    elif opt in ('-a', '--account'):
        options['accout']=arg
    elif opt in ('-h', '--HELP'):
        print( "-b")
        print( "     -b <bench_case>")
        print( "       <bench_case> : benchmark(s) to validate. Accepts wildcards.")
        print( "     DEFAULT : All benchmarks are validated.")
        print( "-o")
        print( "     -o <nb_OMPThreads>")
        print( "       <nb_OMPThreads> : number of OpenMP threads used for the execution")
        print( "     DEFAULT : 4")
        print( "-m")
        print( "     -m <nb_MPIProcs>")
        print( "       <nb_MPIProcs> : number of MPI processes used for the execution")
        print( "     DEFAULT : 4")
        print( "-p")
        print( "     -p <partition name>")
        print( "       <partition name>: partition name on super-computers")
        print( "                         - ruche")
        print( "                         - tornado")
        print( "                         - jollyjumper")
        print( "                         - irene_skylake")
        print( "                         - irene_a64fx")
        print( "-t")
        print( "     -t <max time>")
        print( "       <max time>: format hh:mm:ss")
        print( "-a")
        print( "     -a --account <account id>")
        print( "       <account id>: account/project id given by some super-computer facilities")
        print( "-g")
        print( "     Generates the references")
        print( "-s")
        print( "     Plot differences with references (python -i option required to keep figures on screen)")
        print( "-r")
        print( "     -r <nb_restarts>")
        print( "       <nb_restarts> : number of restarts to run, as long as the simulations provide them.")
        print( "     DEFAULT : 0 (meaning no restarts, only one simulation)")
        print( "-c")
        print( "     Compilation only")
        print( "-k")
        print( "     Compilation using config=... See make help for details")
        print( "-v")
        print( "     Verbose mode")
        print( "-l")
        print( "     Log some performance info in the directory `logs`")
        sys.exit(0)
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
            sys.exit(4)
    elif opt in ('-l', '--LOG'):
        options['log'] = True
        if path.isabs(arg):
            SMILEI_LOGS = arg + s
        else:
            SMILEI_LOGS = INITIAL_DIRECTORY + s + arg + s

# Manage some stuff according to options
MAKE = "make" + (" config=%s"%options['compile_mode'] if options['compile_mode'] else "")

options['max_time_seconds'] = np.sum(np.array(options['max_time'].split(":"),dtype=int)*np.array([3600,60,1]))

if options['generate'] and options['showdiff']:
    usage()
    sys.exit(4)

# Build the list of the requested input files
list_validation = [path.basename(b) for b in glob(parameters['smilei_analyse_path']+"validate_tst*py")]
if options['bench'] == "":
    parameters['benchmarks'] = [path.basename(b) for b in glob(parameters['smilei_benchmark_path']+"tst*py")]
else:
    parameters['benchmarks'] = glob( parameters['smilei_benchmark_path'] + options['bench'] )
    parameters['benchmarks'] = [b.replace(parameters['smilei_benchmark_path'],'') for b in parameters['benchmarks']]
parameters['benchmarks'] = [b for b in parameters['benchmarks'] if "validate_"+b in list_validation]
if not parameters['benchmarks']:
    raise Exception("Input file(s) "+options['bench']+" not found, or without validation file")

if options['verbose']:
    print( "")
    print( "The list of input files to be validated is:\n\t"+"\n\t".join(parameters['benchmarks']))
    print( "")


print(HOSTNAME)

# Define commands according to the host
if LLR in HOSTNAME:
    options['ppn'] = {"jollyjumper":12, "tornado":18}[options['partition']]
    if options['ppn'] % options['omp'] != 0:
        print("Smilei cannot be run with "+str(options['omp'])+" threads on "+HOSTNAME+" and partition "+options['partition'])
        sys.exit(4)
    NODES = int(ceil(options['mpi']/2.))
    NPERSOCKET = 1
    COMPILE_COMMAND = str(MAKE)+' -j '+str(options['ppn'])+' > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
    COMPILE_TOOLS_COMMAND = 'make tables > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
    CLEAN_COMMAND = 'make clean > /dev/null 2>&1'
    RUN_COMMAND = "mpirun --mca mpi_warn_on_fork 0 -mca orte_num_sockets 2 -mca orte_num_cores "+str(options['ppn']) + " -map-by ppr:"+str(NPERSOCKET)+":socket:"+"pe="+str(options['omp']) + " -n "+str(options['mpi'])+" -x OMP_NUM_THREADS -x OMP_SCHEDULE "+WORKDIR_BASE+s+"smilei %s >"+parameters['output_file']+" 2>&1"
    RUN = run_llr
elif "ruche" in HOSTNAME:
    options['ppn'] = 20
    if options['ppn'] < options['omp'] :
        print("Smilei cannot be run with "+str(options['omp'])+" threads on "+HOSTNAME+" and partition "+options['partition'])
        sys.exit(4)
    NODES=int(ceil(options['mpi']/2.))
    NPERSOCKET = 1
    COMPILE_COMMAND = str(MAKE)+' -j 20 machine=ruche > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
    COMPILE_TOOLS_COMMAND = 'make tables > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
    CLEAN_COMMAND = 'make clean > /dev/null 2>&1'
    RUN_COMMAND ="srun "+WORKDIR_BASE+s+"smilei %s >"+parameters['output_file']+" 2>&1"
    RUN = run_ruche
elif POINCARE in HOSTNAME:
    COMPILE_COMMAND = str(MAKE)+' -j 6 > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
    COMPILE_TOOLS_COMMAND = 'make tables > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
    CLEAN_COMMAND = 'module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu ; unset LD_PRELOAD ; export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0 > /dev/null 2>&1;make clean > /dev/null 2>&1'
    RUN_COMMAND = "mpirun -np "+str(options['mpi'])+" "+WORKDIR_BASE+s+"smilei %s >"+parameters['output_file']
    RUN = run_poincare
else: # Local computers
    # Determine the correct MPI command
    mpi_version = str(check_output("mpirun --version", shell=True))
    if re.search("Open MPI", mpi_version, re.I):
        v = re.search("\d\d?\.\d\d?\.\d\d?", mpi_version).group() # Full version number
        v = int(v.split(".")[0]) # Major version number
        if v > 1:
            MPIRUN = "mpirun --oversubscribe -np "
        else:
            MPIRUN = "mpirun -mca btl tcp,sm,self -np "
    else:
        MPIRUN = "mpirun -np "

    COMPILE_COMMAND = str(MAKE)+' -j4 > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
    COMPILE_TOOLS_COMMAND = 'make tables > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
    CLEAN_COMMAND = 'make clean > /dev/null 2>&1'
    RUN_COMMAND = "export OMP_NUM_THREADS="+str(options['omp'])+"; "+MPIRUN+str(options['mpi'])+" "+WORKDIR_BASE+s+"smilei %s >"+parameters['output_file']
    RUN = run_other

# ---------------
# Compilation
# ---------------

if options['verbose']:
    print("---------------------------")
    print("Compiling Smilei")
    print("---------------------------")

# Get state of smilei bin in root folder
os.chdir(parameters['smilei_root'])
STAT_SMILEI_R_OLD = os.stat(SMILEI_R) if path.exists(SMILEI_R) else ' '

# CLEAN
# If no smilei bin in the workdir, or it is older than the one in smilei directory, clean to force compilation
mkdir(WORKDIR_BASE)
if path.exists(SMILEI_R) and (not path.exists(SMILEI_W) or date(SMILEI_W)<date(SMILEI_R)):
    call(CLEAN_COMMAND , shell=True)

try:
    # Remove the compiling errors files
    if path.exists(WORKDIR_BASE+s+COMPILE_ERRORS) :
        os.remove(WORKDIR_BASE+s+COMPILE_ERRORS)
    # Compile
    RUN( COMPILE_COMMAND, parameters['smilei_root'], 'compilation', options, parameters )
    os.rename(COMPILE_OUT_TMP, COMPILE_OUT)
    if STAT_SMILEI_R_OLD!=os.stat(SMILEI_R) or date(SMILEI_W)<date(SMILEI_R): # or date(SMILEI_TOOLS_W)<date(SMILEI_TOOLS_R) :
        # if new bin, archive the workdir (if it contains a smilei bin)
        # and create a new one with new smilei and compilation_out inside
        if path.exists(SMILEI_W): # and path.exists(SMILEI_TOOLS_W):
            workdir_archiv(SMILEI_W)
        copy2(SMILEI_R,SMILEI_W)
        #copy2(SMILEI_TOOLS_R,SMILEI_TOOLS_W)
        if options['compile_only']:
            if options['verbose']:
                print("Smilei validation succeed.")
            exit(0)
    else:
        if options['compile_only'] :
            if options['verbose']:
                print("Smilei validation not needed.")
            exit(0)

except CalledProcessError as e:
    # if compiling errors, archive the workdir (if it contains a smilei bin),
    # create a new one with compilation_errors inside and exit with error code
    workdir_archiv(SMILEI_W)
    os.rename(COMPILE_ERRORS,WORKDIR_BASE+s+COMPILE_ERRORS)
    if options['verbose']:
        print("Smilei validation cannot be done : compilation failed. " + str(e.returncode))
    sys.exit(3)

if options['verbose']:
    print()

# -------------------------------------------
# Define functions and classes for validation
# -------------------------------------------

def loadReference(bench_name):
    try:
        try:
            with open(parameters['smilei_reference_path']+s+bench_name+".txt", 'rb') as f:
                return pickle.load(f, fix_imports=True, encoding='latin1')
        except:
            with open(parameters['smilei_reference_path']+s+bench_name+".txt", 'r') as f:
                return pickle.load(f)
    except:
        print("Unable to find the reference data for "+bench_name)
        sys.exit(1)

def matchesWithReference(data, expected_data, data_name, precision):
    # ok if exactly equal (including strings or lists of strings)
    try:
        if expected_data == data:
            return True
    except:
        pass
    # If numbers:
    try:
        double_data = np.array(np.double(data), ndmin=1)
        if precision is not None:
            error = np.abs( double_data-np.array(np.double(expected_data), ndmin=1) )
            max_error_location = np.unravel_index(np.argmax(error), error.shape)
            max_error = error[max_error_location]
            if max_error < precision:
                return True
            print("Reference quantity '"+data_name+"' does not match the data (required precision "+str(precision)+")")
            print("Max error = "+str(max_error)+" at index "+str(max_error_location))
        else:
            if np.all(double_data == np.double(expected_data)):
                return True
            print("Reference quantity '"+data_name+"' does not match the data")
    except Exception as e:
        print("Reference quantity '"+data_name+"': unable to compare to data")
        print( e )
    return False


# DEFINE A CLASS TO CREATE A REFERENCE
class CreateReference(object):
    def __init__(self, bench_name):
        self.reference_file = parameters['smilei_reference_path']+s+bench_name+".txt"
        self.data = {}

    def __call__(self, data_name, data, precision=None):
        self.data[data_name] = data

    def write(self):
        with open(self.reference_file, "wb") as f:
            pickle.dump(self.data, f, protocol=2)
        size = path.getsize(self.reference_file)
        if size > 1000000:
            print("Reference file is too large ("+str(size)+"B) - suppressing ...")
            os.remove(self.reference_file)
            sys.exit(2)
        if options['verbose']:
            print("Created reference file "+self.reference_file)

# DEFINE A CLASS TO COMPARE A SIMULATION TO A REFERENCE
class CompareToReference(object):
    def __init__(self, bench_name):
        self.data = loadReference(bench_name)

    def __call__(self, data_name, data, precision=None):
        # verify the name is in the reference
        if data_name not in self.data.keys():
            print("Reference quantity '"+data_name+"' not found")
            sys.exit(1)
        expected_data = self.data[data_name]
        if not matchesWithReference(data, expected_data, data_name, precision):
            print("Reference data:")
            print(expected_data)
            print("New data:")
            print(data)
            print()
            global _dataNotMatching
            _dataNotMatching = True

# DEFINE A CLASS TO VIEW DIFFERENCES BETWEEN A SIMULATION AND A REFERENCE
class ShowDiffWithReference(object):
    def __init__(self, bench_name):
        self.data = loadReference(bench_name)

    def __call__(self, data_name, data, precision=None):
        import matplotlib.pyplot as plt
        plt.ion()
        print("Showing differences about '"+data_name+"'")
        print("--------------------------")
        # verify the name is in the reference
        if data_name not in self.data.keys():
            print("\tReference quantity not found")
            expected_data = None
        else:
            expected_data = self.data[data_name]
        print_data = False
        # First, check whether the data matches
        if not matchesWithReference(data, expected_data, data_name, precision):
            global _dataNotMatching
            _dataNotMatching = True
        # try to convert to array
        try:
            data_float = np.array(data, dtype=float)
            expected_data_float = np.array(expected_data, dtype=float)
        # Otherwise, simply print the result
        except:
            print("\tQuantity cannot be plotted")
            print_data = True
            data_float = None
        # Manage array plotting
        if data_float is not None:
            if expected_data is not None and data_float.shape != expected_data_float.shape:
                print("\tReference and new data do not have the same shape: "+str(expected_data_float.shape)+" vs. "+str(data_float.shape))
            if expected_data is not None and data_float.ndim != expected_data_float.ndim:
                print("\tReference and new data do not have the same dimension: "+str(expected_data_float.ndim)+" vs. "+str(data_float.ndim))
                print_data = True
            elif data_float.size == 0:
                print("\t0D quantity cannot be plotted")
                print_data = True
            elif data_float.ndim == 1:
                nplots = 2
                if expected_data is None or data_float.shape != expected_data_float.shape:
                    nplots = 1
                fig = plt.figure()
                fig.suptitle(data_name)
                print("\tPlotting in figure "+str(fig.number))
                ax1 = fig.add_subplot(nplots,1,1)
                ax1.plot( data_float, label="new data" )
                ax1.plot( expected_data_float, label="reference data" )
                ax1.legend()
                if nplots == 2:
                    ax2 = fig.add_subplot(nplots,1,2)
                    ax2.plot( data_float-expected_data_float )
                    ax2.set_title("difference")
            elif data_float.ndim == 2:
                nplots = 3
                if expected_data is None:
                    nplots = 1
                elif data_float.shape != expected_data_float.shape:
                    nplots = 2
                fig = plt.figure()
                fig.suptitle(data_name)
                print("\tPlotting in figure "+str(fig.number))
                ax1 = fig.add_subplot(1,nplots,1)
                im = ax1.imshow( data_float )
                ax1.set_title("new data")
                plt.colorbar(im)
                if nplots > 1:
                    ax2 = fig.add_subplot(1,nplots,2)
                    im = ax2.imshow( expected_data_float )
                    ax2.set_title("reference data")
                    plt.colorbar( im )
                if nplots > 2:
                    ax3 = fig.add_subplot(1,nplots,nplots)
                    im = ax3.imshow( data_float-expected_data_float )
                    ax3.set_title("difference")
                    plt.colorbar( im )
                plt.draw()
                plt.show()
            else:
                print("\t"+str(data_float.ndim)+"D quantity cannot be plotted")
                print_data = True
        # Print data if necessary
        if print_data:
            if expected_data is not None:
                print("\tReference data:")
                print(expected_data)
            print("\tNew data:")
            print(data)


# ---------------
# Run benchmarks
# ---------------

_dataNotMatching = False
for BENCH in parameters['benchmarks'] :

    SMILEI_BENCH = parameters['smilei_benchmark_path'] + BENCH

    # Create the workdir path
    WORKDIR = WORKDIR_BASE+s+'wd_'+path.basename(path.splitext(BENCH)[0])
    mkdir(WORKDIR)
    WORKDIR += s+str(options['mpi'])
    mkdir(WORKDIR)
    WORKDIR += s+str(options['omp'])
    mkdir(WORKDIR)

    # If there are restarts, prepare a Checkpoints block in the namelist
    RESTART_INFO = ""
    if options['nb_restarts'] > 0:
        # Load the namelist
        namelist = happi.openNamelist(SMILEI_BENCH)
        niter = namelist.Main.simulation_time / namelist.Main.timestep
        # If the simulation does not have enough timesteps, change the number of restarts
        if options['nb_restarts'] > niter - 4:
            options['nb_restarts'] = max(0, niter - 4)
            if options['verbose'] :
                print("Not enough timesteps for restarts. Changed to "+str(options['nb_restarts'])+" restarts")
    if options['nb_restarts'] > 0:
        # Find out the optimal dump_step
        dump_step = int( (niter+3.) / (options['nb_restarts']+1) )
        # Prepare block
        if len(namelist.Checkpoints) > 0:
            RESTART_INFO = (" \""
                + "Checkpoints.keep_n_dumps="+str(options['nb_restarts'])+";"
                + "Checkpoints.dump_minutes=0.;"
                + "Checkpoints.dump_step="+str(dump_step)+";"
                + "Checkpoints.exit_after_dump=True;"
                + "Checkpoints.restart_dir=%s;"
                + "\""
            )
        else:
            RESTART_INFO = (" \"Checkpoints("
                + "keep_n_dumps="+str(options['nb_restarts'])+","
                + "dump_minutes=0.,"
                + "dump_step="+str(dump_step)+","
                + "exit_after_dump=True,"
                + "restart_dir=%s,"
                + ")\""
            )
        del namelist

    # Prepare logging
    if options['log']:
        log = Log(SMILEI_LOGS + BENCH + ".log")

    # Loop restarts
    for irestart in range(options['nb_restarts']+1):

        RESTART_WORKDIR = WORKDIR + s + "restart%03d"%irestart

        options['execution'] = True
        if not path.exists(RESTART_WORKDIR):
            os.mkdir(RESTART_WORKDIR)
        elif options['generate']:
            options['execution'] = False

        os.chdir(RESTART_WORKDIR)

        # Copy of the databases
        # For the cases that need a database
        # if BENCH in [
        #         "tst1d_09_rad_electron_laser_collision.py",
        #         "tst1d_10_pair_electron_laser_collision.py",
        #         "tst2d_08_synchrotron_chi1.py",
        #         "tst2d_09_synchrotron_chi0.1.py",
        #         "tst2d_v_09_synchrotron_chi0.1.py",
        #         "tst2d_v_10_multiphoton_Breit_Wheeler.py",
        #         "tst2d_10_multiphoton_Breit_Wheeler.py",
        #         "tst2d_15_qed_cascade_particle_merging.py",
        #         "tst3d_15_magnetic_shower_particle_merging.py"
        #     ]:
        #     try :
        #         # Copy the database
        #         check_call(['cp '+SMILEI_DATABASE+'/*.h5 '+RESTART_WORKDIR], shell=True)
        #     except CalledProcessError:
        #         if options['verbose'] :
        #             print(  "Execution failed to copy databases in ",RESTART_WORKDIR)
        #         sys.exit(2)

        # If there are restarts, adds the Checkpoints block
        arguments = SMILEI_BENCH
        if options['nb_restarts'] > 0:
            if irestart == 0:
                RESTART_DIR = "None"
            else:
                RESTART_DIR = "'"+WORKDIR+s+("restart%03d"%(irestart-1))+s+"'"
            arguments += RESTART_INFO % RESTART_DIR

        # Run smilei
        if options['execution']:
            if options['verbose']:
                print('Running '+BENCH+' on '+HOSTNAME+' with '+str(options['omp'])+'x'+str(options['mpi'])+' OMPxMPI' + ((", restart #"+str(irestart)) if irestart>0 else ""))
            RUN( RUN_COMMAND % arguments, RESTART_WORKDIR, 'execution', options, parameters )

        # Check the output for errors
        errors = []
        search_error = re.compile('error', re.IGNORECASE)
        with open(parameters['output_file'],"r") as fout:
            errors = [line for line in fout if search_error.search(line)]
        if errors:
            if options['verbose']:
                print("")
                print("Errors appeared while running the simulation:")
                print("---------------------------------------------")
                for error in errors:
                    print(error)
            sys.exit(2)

        # Scan some info for logging
        if options['log']:
            log.scan(parameters['output_file'])

    # Append info in log file
    if options['log']:
        log.append()

    os.chdir(WORKDIR)

    # Find the validation script for this bench
    validation_script = parameters['smilei_analyse_path'] + "validate_" + BENCH
    if options['verbose']: print("")
    if not path.exists(validation_script):
        print("Unable to find the validation script "+validation_script)
        sys.exit(1)

    # If required, options['generate'] the references
    if options['generate']:
        if options['verbose']:
            print( '----------------------------------------------------')
            print( 'Generating reference for '+BENCH)
            print( '----------------------------------------------------')
        Validate = CreateReference(BENCH)
        execfile(validation_script)
        Validate.write()

    # Or plot differences with respect to existing references
    elif options['showdiff']:
        if options['verbose']:
            print( '----------------------------------------------------')
            print( 'Viewing differences for '+BENCH)
            print( '----------------------------------------------------')
        Validate = ShowDiffWithReference(BENCH)
        execfile(validation_script)
        if _dataNotMatching:
            print("Benchmark "+BENCH+" did NOT pass")

    # Otherwise, compare to the existing references
    else:
        if options['verbose']:
            print( '----------------------------------------------------')
            print( 'Validating '+BENCH)
            print( '----------------------------------------------------')
        Validate = CompareToReference(BENCH)
        execfile(validation_script)
        if _dataNotMatching:
            sys.exit(1)

    # Clean workdirs, goes here only if succeeded
    os.chdir(WORKDIR_BASE)
    rmtree( WORKDIR_BASE+s+'wd_'+path.basename(path.splitext(BENCH)[0]), True )

    if options['verbose']: print( "")

if _dataNotMatching:
    print( "Errors detected")
else:
    print( "Everything passed")
os.chdir(INITIAL_DIRECTORY)
