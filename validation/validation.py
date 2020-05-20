#!/usr/bin/env python

"""
This script can do three things:
  (1) generate validation reference(s) for given benchmark(s)
  (2) compare benchmark(s) to their reference(s)
  (3) show visually differences between benchmark(s) and their reference(s)

Usage
#######
python validation.py [-c] [-h] [-v] [-b <bench_case>] [-o <nb_OMPThreads>] [-m <nb_MPIProcs>] [-g | -s] [-r <nb_restarts>]

For help on options, try 'python validation.py -h'


Here are listed the files used in the validation process:
#########################################################
The "validation" directory which contains:
  - the "references" directory with one file for each benchmark
  - the "analyses" directory with one validation file for each benchmark
  - a "workdirs" directory, created during the validation process
  - archived "workdirs" for previous versions of smilei

A "workdirs" contains:
 - the smilei binary : "smilei"
 - the compilation output : "compilation_output"
 - a directory wd_<input_file>/<o>/<m> directory, containing the output files
 or
 - the compilation errors file : "compilation_errors"

The different steps of the script are the folllowing :
######################################################
Compilation step :
+++++++++++++++
If the "workdirs" directory lacks a smilei binary, or it is too old),
then the "workdirs" is backed up, and a new compilation occurs.
If compiling errors occur, "compilation_errors" is created and the script exits with status 3.

Execution step :
+++++++++++++++
If wd_<input_file>/<o>/<m> does not exist then:
- it is created
- smilei is executed in that directory for the requested benchmark
- if execution fails, the script exits with status 2

Validation step :
+++++++++++++++
Loops through all requested benchmarks
	Runs the benchmark in the workdir
	If requested to generate references:
		Executes the "validate_*" script and stores the result as reference data
	If requested to compare to previous references
		Executes the "validate_*" script and compares the result to the reference data
	If requested to show differences to previous references
		Executes the "validate_*" script and plots the result vs. the reference data

Exit status of the script
+++++++++++++++++++++++++
0  validated
1  validation fails
2  execution fails
3  compilation fails
4  bad option

Remark
+++++++
This script may run anywhere: you can define a SMILEI_ROOT environment variable
"""


# IMPORTS
import sys, os, re, json
import pickle
import numpy as np
from os import path
from time import sleep, ctime, strftime
from math import ceil
from glob import glob
from shutil import rmtree, copy2
from getopt import getopt, GetoptError
from inspect import stack
from socket import gethostname
from subprocess import call, check_call, check_output, CalledProcessError
s = os.sep
INITIAL_DIRECTORY = os.getcwd()

def mkdir(dir):
	if not path.exists(dir):
		os.mkdir(dir)

# DEFINE THE execfile function for python3
try:
	execfile
except:
	def execfile(file):
		exec(compile(open(file).read(), file, 'exec'), globals())

# DEFINE THE raw_input function for python3
try:
	raw_input
except:
	raw_input = input

# SMILEI PATH VARIABLES
if "SMILEI_ROOT" in os.environ :
	SMILEI_ROOT=os.environ["SMILEI_ROOT"]+s
else :
	SMILEI_ROOT = path.dirname(path.abspath(stack()[0][1]))+s+".."+s
	#SMILEI_ROOT = os.getcwd()+s+".."+s
SMILEI_ROOT = path.abspath(SMILEI_ROOT)+s
SMILEI_SCRIPTS = SMILEI_ROOT+"scripts"+s
SMILEI_VALIDATION = SMILEI_ROOT+"validation"+s
SMILEI_REFERENCES = SMILEI_VALIDATION+"references"+s
SMILEI_BENCHS = SMILEI_ROOT+"benchmarks"+s

# SCRIPTS VARIABLES
EXEC_SCRIPT = 'exec_script.sh'
EXEC_SCRIPT_OUT = 'exec_script.out'
SMILEI_EXE_OUT = 'smilei_exe.out'

# Get the current version of Smilei
os.chdir(SMILEI_ROOT)
gitversion = check_output( "echo `git log -n 1 --format=%h`-", shell=True ).decode()[:-1]
if 'CI_COMMIT_BRANCH' in os.environ:
	gitversion += os.environ['CI_COMMIT_BRANCH']
else:
	gitversion += check_output("echo `git rev-parse --abbrev-ref HEAD`", shell=True ).decode()[:-1]
os.chdir(INITIAL_DIRECTORY)

# Load the happi module
sys.path.insert(0, SMILEI_ROOT)
import happi

# OTHER VARIABLES
POINCARE = "poincare"
LLR = "llrlsi-gw"
PARTITION = "jollyjumper"
HOSTNAME = gethostname()

# DIR VARIABLES
WORKDIR = ""

# DEFAULT VALUES FOR OPTIONS
OMP = 12
MPI = 4
PPN = 12
EXECUTION = False
VERBOSE = False
BENCH=""
COMPILE_ONLY = False
GENERATE = False
SHOWDIFF = False
nb_restarts = 0
COMPILE_MODE=""
LOG = False

# TO PRINT USAGE
def usage():
	print( 'Usage: validation.py [-c] [-h] [-v] [-b <bench_case>] [-o <nb_OMPThreads>] [-m <nb_MPIProcs>] [-g | -s] [-r <nb_restarts>] [-k <compile_mode>] [-p <partition name>] [-l <logs_folder>]' )
	print( '    Try `validation.py -h` for more details' )

# GET COMMAND-LINE OPTIONS
try:
	options, remainder = getopt(
		sys.argv[1:],
		'o:m:b:r:k:p:gshvcl:',
		['OMP=', 'MPI=', 'BENCH=', 'RESTARTS=', 'PARTITION=', 'GENERATE', 'SHOW', 'HELP', 'VERBOSE', 'COMPILE_ONLY', 'COMPILE_MODE=', 'LOG='])
except GetoptError as err:
	usage()
	sys.exit(4)

# PROCESS THE OPTIONS
for opt, arg in options:
	if opt in ('-o', '--OMP'):
		EXECUTION = True
		OMP = int(arg)
	elif opt in ('-m', '--MPI'):
		EXECUTION = True
		MPI = int(arg)
	elif opt in ('-b', '--BENCH'):
		BENCH = arg
	elif opt in ('-p', '--PARTITION'):
		PARTITION = arg
	elif opt in ('-c', '--COMPILE_ONLY'):
		COMPILE_ONLY=True
	elif opt in ('-k', '--COMPILE_MODE'):
		COMPILE_MODE=arg
	elif opt in ('-h', '--HELP'):
		print( "-b")
		print( "     -b <bench_case>")
		print( "       <bench_case> : benchmark(s) to validate. Accepts wildcards.")
		print( "       <bench_case>=? : ask input for a benchmark")
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
		GENERATE = True
	elif opt in ('-s', '--SHOW'):
		SHOWDIFF = True
	elif opt in ('-v', '--VERBOSE'):
		VERBOSE = True
	elif opt in ('-r', '--RESTARTS'):
		try:
			nb_restarts = int(arg)
			if nb_restarts < 0: raise
		except:
			print("Error: the number of restarts (option -r) must be a positive integer")
			sys.exit(4)
	elif opt in ('-l', '--LOG'):
		LOG = True
		if path.isabs(arg):
			SMILEI_LOGS = arg + s
		else:
			SMILEI_LOGS = INITIAL_DIRECTORY + s + arg + s

if GENERATE and SHOWDIFF:
	usage()
	sys.exit(4)

# Build the list of the requested input files
list_bench = [path.basename(b) for b in glob(SMILEI_BENCHS+"tst*py")]
list_validation = [path.basename(b) for b in glob(SMILEI_VALIDATION+"analyses"+s+"validate_tst*py")]
list_bench = [b for b in list_bench if "validate_"+b in list_validation]
if BENCH == "":
	SMILEI_BENCH_LIST = list_bench
elif BENCH == "?":
	VERBOSE = True
	os.chdir(SMILEI_SCRIPTS)
	#- Propose the list of all the input files
	print( '\n'.join(list_bench))
	#- Choose an input file name in the list
	print( 'Enter an input file from the list above:')
	BENCH = raw_input()
	SMILEI_BENCH_LIST = [ BENCH ]
	while BENCH not in list_bench:
		print( "Input file "+BENCH+" invalid. Try again.")
		BENCH = raw_input()
		SMILEI_BENCH_LIST = [ BENCH ]
elif BENCH in list_bench:
	SMILEI_BENCH_LIST = [ BENCH ]
elif glob( SMILEI_BENCHS+BENCH ):
	BENCH = glob( SMILEI_BENCHS+BENCH )
	list_all = glob(SMILEI_BENCHS+"tst*py")
	for b in BENCH:
		if b not in list_all:
			if VERBOSE:
				print( "Input file "+b+" invalid.")
			sys.exit(4)
		SMILEI_BENCH_LIST= []
		for b in BENCH:
				if b.replace(SMILEI_BENCHS,'') in list_bench:
					SMILEI_BENCH_LIST.append( b.replace(SMILEI_BENCHS,'') )
		BENCH = SMILEI_BENCH_LIST
else:
	if VERBOSE:
		print( "Input file "+BENCH+" invalid.")
	sys.exit(4)

if VERBOSE :
	print( "")
	print( "The list of input files to be validated is:\n\t"+"\n\t".join(SMILEI_BENCH_LIST))
	print( "")

# GENERIC FUNCTION FOR WORKDIR ORGANIZATION

def date(BIN_NAME):
	statbin = os.stat(BIN_NAME)
	return statbin.st_ctime
def date_string(BIN_NAME):
	date_integer = date(BIN_NAME)
	date_time = ctime(date_integer)
	return date_time.replace(" ","-")
def workdir_archiv(BIN_NAME) :
	if path.exists(SMILEI_W):
		ARCH_WORKDIR = WORKDIR_BASE+'_'+date_string(SMILEI_W)
		os.rename(WORKDIR_BASE, ARCH_WORKDIR)
		mkdir(WORKDIR_BASE)

# PLATFORM-DEPENDENT INSTRUCTIONS FOR RUNNING PARALLEL COMMAND
def RUN_POINCARE(command, dir):
	# Create script
	with open(EXEC_SCRIPT, 'w') as exec_script_desc:
		print( "ON POINCARE NOW")
		exec_script_desc.write(
			"# environnement \n"
			+"module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu 2>&1 > /dev/null\n"
			+"unset LD_PRELOAD\n"
			+"export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0\n"
			+"# \n"
			+"# execution \n"
			+"export OMP_NUM_THREADS="+str(OMP)+"\n"
			+command+" \n"
			+"exit $?  "
		)
	# Run command
	COMMAND = "/bin/bash "+EXEC_SCRIPT+" > "+EXEC_SCRIPT_OUT+" 2>&1"
	try:
		check_call(COMMAND, shell=True)
	except CalledProcessError:
		# if execution fails, exit with exit status 2
		if VERBOSE :
			print(  "Execution failed for command `"+command+"`")
			COMMAND = "/bin/bash cat "+SMILEI_EXE_OUT
			try :
				check_call(COMMAND, shell=True)
			except CalledProcessError:
				print(  "cat command failed")
				sys.exit(2)
		if dir==WORKDIR:
			os.chdir(WORKDIR_BASE)
			rmtree(WORKDIR)
		sys.exit(2)

def RUN_LLR(command, dir):
	"""
	Run the command `command` on the LLR system.

	Inputs:
	- command: command to run
	- dir: working directory
	"""
	EXIT_STATUS="100"
	exit_status_fd = open(dir+s+"exit_status_file", "w")
	exit_status_fd.write(str(EXIT_STATUS))
	exit_status_fd.close()
	# Create script
	with open(EXEC_SCRIPT, 'w') as exec_script_desc:
		#NODES=((int(MPI)*int(OMP)-1)/24)+1
		NODES=int(ceil(MPI/2.))
		if (PARTITION == "jollyjumper"):
			PPN = 24
		elif (PARTITION == "tornado"):
			PPN = 36
		exec_script_desc.write(
			"#PBS -l nodes="+str(NODES)+":ppn="+str(PPN)+" \n"
			+"#PBS -q default \n"
			+"#PBS -j oe\n"
			+"module purge\n"
			+"unset MODULEPATH;\n"
			+"module use /opt/exp_soft/vo.llr.in2p3.fr/modulefiles_el7\n"
			+"module load compilers/icc/17.4.196\n"
			+"module load python/2.7.9\n"
			+"module load mpi/openmpi/1.6.5-ib-icc\n"
			+"module load hdf5/1.8.19-icc-omp1.6.5\n"
			+"module load h5py/hdf5_1.8.19-icc-omp1.6.5-py2.7.9\n"
			+"module load compilers/gcc/4.9.2\n"
			+"export OMP_NUM_THREADS="+str(OMP)+" \n"
			+"export OMP_SCHEDULE=DYNAMIC \n"
			+"export KMP_AFFINITY=verbose \n"
			+"export PATH=$PATH:/opt/exp_soft/vo.llr.in2p3.fr/GALOP/beck \n"
			+"module load fftw/3.3.7-opm-1.6.5-icc-17 \n"
			+"export LIBPXR=/home/llr/galop/derouil/applications/picsar/lib \n"
			+"export LD_LIBRARY_PATH=$LIBPXR:$LD_LIBRARY_PATH \n"
			+"ulimit -s unlimited \n"
			+"#Specify the number of sockets per node in -mca orte_num_sockets \n"
			+"#Specify the number of cores per sockets in -mca orte_num_cores \n"
			+"cd "+dir+" \n"
			+"module list 2> module.log\n"
			+command+" \n"
			+"echo $? > exit_status_file \n"
		)
	# Run command
	if (PARTITION=="jollyjumper"):
		COMMAND = "PBS_DEFAULT=llrlsi-jj.in2p3.fr qsub  "+EXEC_SCRIPT
	elif (PARTITION=="tornado"):
		COMMAND = "PBS_DEFAULT=poltrnd.in2p3.fr qsub  "+EXEC_SCRIPT
	try:
		check_call(COMMAND, shell=True)
	except CalledProcessError:
		# if command qsub fails, exit with exit status 2
		#Retry once in case the server was rebooting
		if VERBOSE :
			print(  "qsub command failed once: `"+COMMAND+"`")
			print(  "Wait and retry")
		sleep(10)
		try:
			check_call(COMMAND, shell=True)
		except CalledProcessError:
			if dir==WORKDIR:
				os.chdir(WORKDIR_BASE)
				rmtree(WORKDIR)
			if VERBOSE :
				print(  "qsub command failed twice: `"+COMMAND+"`")
				print(  "Exit")
			sys.exit(2)
	if VERBOSE:
		print( "Submitted job with command `"+command+"`")
	while ( EXIT_STATUS == "100" ) :
		sleep(5)
		exit_status_fd = open(dir+s+"exit_status_file", "r+")
		EXIT_STATUS = exit_status_fd.readline()
		exit_status_fd.close()
	if ( int(EXIT_STATUS) != 0 )  :
		if VERBOSE :
			print(  "Execution failed for command `"+command+"`")
			COMMAND = "cat "+SMILEI_EXE_OUT
			try :
				check_call(COMMAND, shell=True)
			except CalledProcessError:
				print(  "cat command failed")
				sys.exit(2)
		sys.exit(2)

def RUN_OTHER(command, dir):
	"""
	Run the command `command` on an arbitrary system.

	Inputs:
	- command: command to run
	- dir: working directory
	"""
	try :
		if VERBOSE:
			print( "Trying command `"+command+"`")
		check_call(command, shell=True)
	except CalledProcessError:
		if VERBOSE :
			print(  "Execution failed for command `"+command+"`")
		sys.exit(2)


# SET DIRECTORIES
if VERBOSE :
  print( "Compiling Smilei")

os.chdir(SMILEI_ROOT)
WORKDIR_BASE = SMILEI_ROOT+"validation"+s+"workdirs"
SMILEI_W = WORKDIR_BASE+s+"smilei"
SMILEI_R = SMILEI_ROOT+s+"smilei"
if path.exists(SMILEI_R):
	STAT_SMILEI_R_OLD = os.stat(SMILEI_R)
else :
	STAT_SMILEI_R_OLD = ' '
SMILEI_TOOLS_W = WORKDIR_BASE+s+"smilei_tables"
SMILEI_TOOLS_R = SMILEI_ROOT+s+"smilei_tables"
if path.exists(SMILEI_TOOLS_R):
	STAT_SMILEI_TOOLS_R_OLD = os.stat(SMILEI_TOOLS_R)
else :
	STAT_SMILEI_TOOLS_R_OLD = ' '
COMPILE_ERRORS=WORKDIR_BASE+s+'compilation_errors'
COMPILE_OUT=WORKDIR_BASE+s+'compilation_out'
COMPILE_OUT_TMP=WORKDIR_BASE+s+'compilation_out_temp'

MAKE='make'
if COMPILE_MODE:
	MAKE += " config="+COMPILE_MODE

# Find commands according to the host
if LLR in HOSTNAME :
	if (PARTITION=="jollyjumper"):
		PPN = 12
	elif (PARTITION=="tornado"):
		PPN = 18
	if PPN % OMP != 0:
		print(  "Smilei cannot be run with "+str(OMP)+" threads on "+HOSTNAME+" and partition "+PARTITION)
		sys.exit(4)
	NODES=int(ceil(MPI/2.))
	NPERSOCKET = 1
	COMPILE_COMMAND = str(MAKE)+' -j 12 > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
	COMPILE_COMMAND = str(MAKE)+' -j '+str(PPN)+' > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
	COMPILE_TOOLS_COMMAND = 'make tables > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
	CLEAN_COMMAND = 'make clean > /dev/null 2>&1'
	RUN_COMMAND = "mpirun -mca orte_num_sockets 2 -mca orte_num_cores "+str(PPN)+" -cpus-per-proc "+str(OMP)+" --npersocket "+str(NPERSOCKET)+" -n "+str(MPI)+" -x OMP_NUM_THREADS -x OMP_SCHEDULE "+WORKDIR_BASE+s+"smilei %s >"+SMILEI_EXE_OUT+" 2>&1"
	RUN = RUN_LLR
elif POINCARE in HOSTNAME :
	#COMPILE_COMMAND = 'module load intel/15.0.0 openmpi hdf5/1.8.10_intel_openmpi python gnu > /dev/null 2>&1;make -j 6 > compilation_out_temp 2>'+COMPILE_ERRORS
	#CLEAN_COMMAND = 'module load intel/15.0.0 openmpi hdf5/1.8.10_intel_openmpi python gnu > /dev/null 2>&1;make clean > /dev/null 2>&1'
	COMPILE_COMMAND = str(MAKE)+' -j 6 > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
	COMPILE_COMMAND = str(MAKE)+' -j 6 > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
	COMPILE_TOOLS_COMMAND = 'make tables > '+COMPILE_OUT_TMP+' 2>'+COMPILE_ERRORS
	CLEAN_COMMAND = 'module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu ; unset LD_PRELOAD ; export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0 > /dev/null 2>&1;make clean > /dev/null 2>&1'
	RUN_COMMAND = "mpirun -np "+str(MPI)+" "+WORKDIR_BASE+s+"smilei %s >"+SMILEI_EXE_OUT
	RUN = RUN_POINCARE
# Local computers
else:
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
	RUN_COMMAND = "export OMP_NUM_THREADS="+str(OMP)+"; "+MPIRUN+str(MPI)+" "+WORKDIR_BASE+s+"smilei %s >"+SMILEI_EXE_OUT
	RUN = RUN_OTHER

# CLEAN
# If the workdir does not contains a smilei bin, or it contains one older than the the smilei bin in directory smilei, force the compilation in order to generate the compilation_output
mkdir(WORKDIR_BASE)
if path.exists(SMILEI_R) and (not path.exists(SMILEI_W) or date(SMILEI_W)<date(SMILEI_R)):
	call(CLEAN_COMMAND , shell=True)

# COMPILE
try :
	# Remove the compiling errors files
	if path.exists(WORKDIR_BASE+s+COMPILE_ERRORS) :
		os.remove(WORKDIR_BASE+s+COMPILE_ERRORS)
	# Compile
	RUN( COMPILE_COMMAND, SMILEI_ROOT )
	#RUN( COMPILE_TOOLS_COMMAND, SMILEI_ROOT )
	os.rename(COMPILE_OUT_TMP, COMPILE_OUT)
	if STAT_SMILEI_R_OLD!=os.stat(SMILEI_R) or date(SMILEI_W)<date(SMILEI_R): # or date(SMILEI_TOOLS_W)<date(SMILEI_TOOLS_R) :
		# if new bin, archive the workdir (if it contains a smilei bin)
		# and create a new one with new smilei and compilation_out inside
		if path.exists(SMILEI_W): # and path.exists(SMILEI_TOOLS_W):
			workdir_archiv(SMILEI_W)
		copy2(SMILEI_R,SMILEI_W)
		#copy2(SMILEI_TOOLS_R,SMILEI_TOOLS_W)
		if COMPILE_ONLY:
			if VERBOSE:
				print( "Smilei validation succeed.")
			exit(0)
	else:
		if COMPILE_ONLY :
			if VERBOSE:
				print(  "Smilei validation not needed.")
			exit(0)

except CalledProcessError as e:
	# if compiling errors, archive the workdir (if it contains a smilei bin),
	# create a new one with compilation_errors inside and exit with error code
	workdir_archiv(SMILEI_W)
	os.rename(COMPILE_ERRORS,WORKDIR_BASE+s+COMPILE_ERRORS)
	if VERBOSE:
		print( "Smilei validation cannot be done : compilation failed. " + str(e.returncode))
	sys.exit(3)
if VERBOSE: print( "")


def findReference(bench_name):
	try:
		try:
			with open(SMILEI_REFERENCES+s+bench_name+".txt", 'rb') as f:
				return pickle.load(f, fix_imports=True, encoding='latin1')
		except:
			with open(SMILEI_REFERENCES+s+bench_name+".txt", 'r') as f:
				return pickle.load(f)
	except:
		print( "Unable to find the reference data for "+bench_name)
		sys.exit(1)

def matchesWithReference(data, expected_data, data_name, precision):
	# ok if exactly equal (including strings or lists of strings)
	try:
		if expected_data == data:
			return True
	except: pass
	# If numbers:
	try:
		double_data = np.array(np.double(data), ndmin=1)
		if precision is not None:
			error = np.abs( double_data-np.array(np.double(expected_data), ndmin=1) )
			max_error_location = np.unravel_index(np.argmax(error), error.shape)
			max_error = error[max_error_location]
			if max_error < precision:
				return True
			print( "Reference quantity '"+data_name+"' does not match the data (required precision "+str(precision)+")")
			print( "Max error = "+str(max_error)+" at index "+str(max_error_location))
		else:
			if np.all(double_data == np.double(expected_data)):
				return True
			print( "Reference quantity '"+data_name+"' does not match the data")
	except Exception as e:
		print( "Reference quantity '"+data_name+"': unable to compare to data")
		print( e )
	return False


# DEFINE A CLASS TO CREATE A REFERENCE
class CreateReference(object):
	def __init__(self, bench_name):
		self.reference_file = SMILEI_REFERENCES+s+bench_name+".txt"
		self.data = {}

	def __call__(self, data_name, data, precision=None):
		self.data[data_name] = data

	def write(self):
		with open(self.reference_file, "wb") as f:
			pickle.dump(self.data, f, protocol=2)
		size = path.getsize(self.reference_file)
		if size > 1000000:
			print( "Reference file is too large ("+str(size)+"B) - suppressing ...")
			os.remove(self.reference_file)
			sys.exit(2)
		if VERBOSE:
			print( "Created reference file "+self.reference_file)

# DEFINE A CLASS TO COMPARE A SIMULATION TO A REFERENCE
class CompareToReference(object):
	def __init__(self, bench_name):
		self.data = findReference(bench_name)
	
	def __call__(self, data_name, data, precision=None):
		# verify the name is in the reference
		if data_name not in self.data.keys():
			print( "Reference quantity '"+data_name+"' not found")
			sys.exit(1)
		expected_data = self.data[data_name]
		if not matchesWithReference(data, expected_data, data_name, precision):
			print( "Reference data:")
			print( expected_data)
			print( "New data:")
			print( data)
			print( "" )
			global _dataNotMatching
			_dataNotMatching = True
			
# DEFINE A CLASS TO VIEW DIFFERENCES BETWEEN A SIMULATION AND A REFERENCE
class ShowDiffWithReference(object):
	def __init__(self, bench_name):
		self.data = findReference(bench_name)
	
	def __call__(self, data_name, data, precision=None):
		import matplotlib.pyplot as plt
		plt.ion()
		print( "Showing differences about '"+data_name+"'")
		print( "--------------------------")
		# verify the name is in the reference
		if data_name not in self.data.keys():
			print( "\tReference quantity not found")
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
			print( "\tQuantity cannot be plotted")
			print_data = True
			data_float = None
		# Manage array plotting
		if data_float is not None:
			if expected_data is not None and data_float.shape != expected_data_float.shape:
				print( "\tReference and new data do not have the same shape: "+str(expected_data_float.shape)+" vs. "+str(data_float.shape))
			if expected_data is not None and data_float.ndim != expected_data_float.ndim:
				print( "\tReference and new data do not have the same dimension: "+str(expected_data_float.ndim)+" vs. "+str(data_float.ndim))
				print_data = True
			elif data_float.size == 0:
				print( "\t0D quantity cannot be plotted")
				print_data = True
			elif data_float.ndim == 1:
				nplots = 2
				if expected_data is None or data_float.shape != expected_data_float.shape:
					nplots = 1
				fig = plt.figure()
				fig.suptitle(data_name)
				print( "\tPlotting in figure "+str(fig.number))
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
				print( "\tPlotting in figure "+str(fig.number))
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
				print( "\t"+str(data_float.ndim)+"D quantity cannot be plotted")
				print_data = True
		# Print data if necessary
		if print_data:
			if expected_data is not None:
				print( "\tReference data:")
				print( expected_data)
			print( "\tNew data:")
			print( data)


# DEFINE A CLASS FOR LOGGING DATA
class Log:
	pattern1 = re.compile(""
		+"[\n\t\s]+(Time[ _]in[ _]time[ _]loop)\s+([.0-9]+)\s+([<.0-9]+)\% coverage"
		+"([\n\t\s]+([\w ]+)\s+([.0-9]+)\s+([<.0-9]+)\%){2,15}"
	)
	pattern2 = re.compile(""
		+"[\t\s]+([\w ]+):?\s+([.0-9]+)\s+([<.0-9]+)\%"
	)
	
	def __init__(self, log_file):
		mkdir(SMILEI_LOGS)
		self.log_file = log_file
		self.data = {}
	
	def scan(self, output):
		# Open file and find the pattern
		with open(output, 'r') as f:
			text = f.read()
		found = re.search(self.pattern1, text)
		if not found:
			print( "WARNING: Unable to log data from "+output)
			return
		lines = found.string[found.start():found.end()].split("\n")[1:]
		matches = [re.search(self.pattern2, l).groups() for l in lines]
		# Get timers values and add to current timers
		for m in matches:
			key = m[0].replace(" ", "")
			if key == "Time_in_time_loop": key = "Timeintimeloop"
			value = float(m[1])
			if key in self.data:
				self.data[key] += value
			else:
				self.data[key] = value
	
	def append(self):
		if self.data:
			# Append commit and date to current data
			self.data["commit"] = gitversion
			self.data["date"] = strftime("%Y_%m_%d_%H:%M:%S")
			# Open previous database
			try:
				with open(self.log_file, 'r') as f:
					db = json.load(f)
			except:
				db = {}
			maxlen = max([len(v) for v in db.values()] or [0])
			# Update the database
			for k,v in self.data.items():
				if k in db:
					db[k] += [None]*(maxlen-len(db[k])) + [v]
				else:
					db[k] = maxlen*[None] + [v]
			# Overwrite the file
			with open(self.log_file, 'w+') as f:
				json.dump(db, f)
	

# RUN THE BENCHMARKS
_dataNotMatching = False
for BENCH in SMILEI_BENCH_LIST :

	SMILEI_BENCH = SMILEI_BENCHS + BENCH

	# CREATE THE WORKDIR CORRESPONDING TO THE INPUT FILE
	WORKDIR = WORKDIR_BASE+s+'wd_'+path.basename(path.splitext(BENCH)[0])
	mkdir(WORKDIR)

	WORKDIR += s+str(MPI)
	mkdir(WORKDIR)

	WORKDIR += s+str(OMP)
	mkdir(WORKDIR)

	# If there are restarts, prepare a Checkpoints block to the namelist
	RESTART_INFO = ""
	if nb_restarts > 0:
		# Load the namelist
		namelist = happi.openNamelist(SMILEI_BENCH)
		niter = namelist.Main.simulation_time / namelist.Main.timestep
		# If the simulation does not have enough timesteps, change the number of restarts
		if nb_restarts > niter - 4:
			nb_restarts = max(0, niter - 4)
			if VERBOSE :
				print("Not enough timesteps for restarts. Changed to "+str(nb_restarts)+" restarts")
	if nb_restarts > 0:
		# Find out the optimal dump_step
		dump_step = int( (niter+3.) / (nb_restarts+1) )
		# Prepare block
		if len(namelist.Checkpoints) > 0:
			RESTART_INFO = (" \""
				+ "Checkpoints.keep_n_dumps="+str(nb_restarts)+";"
				+ "Checkpoints.dump_minutes=0.;"
				+ "Checkpoints.dump_step="+str(dump_step)+";"
				+ "Checkpoints.exit_after_dump=True;"
				+ "Checkpoints.restart_dir=%s;"
				+ "\""
			)
		else:
			RESTART_INFO = (" \"Checkpoints("
				+ "keep_n_dumps="+str(nb_restarts)+","
				+ "dump_minutes=0.,"
				+ "dump_step="+str(dump_step)+","
				+ "exit_after_dump=True,"
				+ "restart_dir=%s,"
				+ ")\""
			)
		del namelist

	# Prepare logging
	if LOG:
		log = Log(SMILEI_LOGS + BENCH + ".log")
	
	# Loop restarts
	for irestart in range(nb_restarts+1):

		RESTART_WORKDIR = WORKDIR + s + "restart%03d"%irestart

		EXECUTION = True
		if not path.exists(RESTART_WORKDIR):
			os.mkdir(RESTART_WORKDIR)
		elif GENERATE:
			EXECUTION = False

		os.chdir(RESTART_WORKDIR)

		# Copy of the databases
		# For the cases that need a database
		# if BENCH in [
		# 		"tst1d_09_rad_electron_laser_collision.py",
		# 		"tst1d_10_pair_electron_laser_collision.py",
		# 		"tst2d_08_synchrotron_chi1.py",
		# 		"tst2d_09_synchrotron_chi0.1.py",
		# 		"tst2d_v_09_synchrotron_chi0.1.py",
		# 		"tst2d_v_10_multiphoton_Breit_Wheeler.py",
		# 		"tst2d_10_multiphoton_Breit_Wheeler.py",
		# 		"tst2d_15_qed_cascade_particle_merging.py",
		# 		"tst3d_15_magnetic_shower_particle_merging.py"
		# 	]:
		# 	try :
		# 		# Copy the database
		# 		check_call(['cp '+SMILEI_DATABASE+'/*.h5 '+RESTART_WORKDIR], shell=True)
		# 	except CalledProcessError:
		# 		if VERBOSE :
		# 			print(  "Execution failed to copy databases in ",RESTART_WORKDIR)
		# 		sys.exit(2)

		# If there are restarts, adds the Checkpoints block
		SMILEI_NAMELISTS = SMILEI_BENCH
		if nb_restarts > 0:
			if irestart == 0:
				RESTART_DIR = "None"
			else:
				RESTART_DIR = "'"+WORKDIR+s+("restart%03d"%(irestart-1))+s+"'"
			SMILEI_NAMELISTS += RESTART_INFO % RESTART_DIR

		# RUN smilei IF EXECUTION IS TRUE
		if EXECUTION :
			if VERBOSE:
				print( 'Running '+BENCH+' on '+HOSTNAME+' with '+str(OMP)+'x'+str(MPI)+' OMPxMPI' + ((", restart #"+str(irestart)) if irestart>0 else ""))
			RUN( RUN_COMMAND % SMILEI_NAMELISTS, RESTART_WORKDIR)

		# CHECK THE OUTPUT FOR ERRORS
		errors = []
		search_error = re.compile('error', re.IGNORECASE)
		with open(SMILEI_EXE_OUT,"r") as fout:
			for line in fout:
				if search_error.search(line):
					errors += [line]
		if errors:
			if VERBOSE :
				print( "")
				print("Errors appeared while running the simulation:")
				print("---------------------------------------------")
				for error in errors:
					print(error)
			sys.exit(2)

		# Scan some info for logging
		if LOG:
			log.scan(SMILEI_EXE_OUT)
	
	# Append info in log file
	if LOG:
		log.append()
	
	os.chdir(WORKDIR)

	# FIND THE VALIDATION SCRIPT FOR THIS BENCH
	validation_script = SMILEI_VALIDATION + "analyses" + s + "validate_"+BENCH
	if VERBOSE: print( "")
	if not path.exists(validation_script):
		print( "Unable to find the validation script "+validation_script)
		sys.exit(1)
	
	# IF REQUIRED, GENERATE THE REFERENCES
	if GENERATE:
		if VERBOSE:
			print( '----------------------------------------------------')
			print( 'Generating reference for '+BENCH)
			print( '----------------------------------------------------')
		Validate = CreateReference(BENCH)
		execfile(validation_script)
		Validate.write()

	# OR PLOT DIFFERENCES WITH RESPECT TO EXISTING REFERENCES
	elif SHOWDIFF:
		if VERBOSE:
			print( '----------------------------------------------------')
			print( 'Viewing differences for '+BENCH)
			print( '----------------------------------------------------')
		Validate = ShowDiffWithReference(BENCH)
		execfile(validation_script)
		if _dataNotMatching:
			print("Benchmark "+BENCH+" did NOT pass")
	
	# OTHERWISE, COMPARE TO THE EXISTING REFERENCES
	else:
		if VERBOSE:
			print( '----------------------------------------------------')
			print( 'Validating '+BENCH)
			print( '----------------------------------------------------')
		Validate = CompareToReference(BENCH)
		execfile(validation_script)
		if _dataNotMatching:
			sys.exit(1)
		# DATA DID NOT MATCH REFERENCES
	
	# CLEAN WORKDIRS, GOES HERE ONLY IF SUCCEED
	os.chdir(WORKDIR_BASE)
	rmtree( WORKDIR_BASE+s+'wd_'+path.basename(path.splitext(BENCH)[0]), True )

	if VERBOSE: print( "")

if _dataNotMatching:
	print( "Errors detected")
else:
	print( "Everything passed")
os.chdir(INITIAL_DIRECTORY)
