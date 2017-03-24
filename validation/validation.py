##!/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/bin/python 

"""
This script can do two things:
  (1) generate validation reference(s) for given benchmark(s)
  (2) compare benchmark(s) to their reference(s)

Usage
#######
python validation.py [-c] [-h] [-b <bench_case> [-o <nb_OMPThreads>] [-m <nb_MPIProcs>] [-g] [-v]]

For help on options, try 'python validation.py -h'


Here are listed the files used in the validation process:
#########################################################
The "validation" directory which contains:
  - the "references" directory with one file for each benchmark
  - validation files (prefixed with "validate_") for each benchmark
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
import sys, os, re, glob
import shutil, getopt, inspect, socket, pickle
from subprocess import check_call,CalledProcessError,call
s = os.sep

# SMILEI PATH VARIABLES
if "SMILEI_ROOT" in os.environ :
	SMILEI_ROOT=os.environ["SMILEI_ROOT"]+s
else :
	SMILEI_ROOT = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))+s+".."+s
	#SMILEI_ROOT = os.getcwd()+s+".."+s
SMILEI_ROOT = os.path.abspath(SMILEI_ROOT)+s
SMILEI_SCRIPTS = SMILEI_ROOT+"scripts"+s
SMILEI_VALIDATION = SMILEI_ROOT+"validation"+s
SMILEI_REFERENCES = SMILEI_VALIDATION+"references"+s
SMILEI_BENCHS = SMILEI_ROOT+"benchmarks"+s

# Load the Smilei module
execfile(SMILEI_SCRIPTS+"Diagnostics.py")

# OTHER VARIABLES
POINCARE = "poincare"
JOLLYJUMPER = "llrlsi-gw"
HOSTNAME = socket.gethostname()

# DEFAULT VALUES FOR OPTIONS
OMP = 2
MPI = 2
EXECUTION = False
VERBOSE = False
BENCH=""
COMPILE_ONLY = False
GENERATE = False

# TO PRINT USAGE
def usage():
	print 'Usage: validation.py [-c] [-h] [-b <bench_case> [-o <nb_OMPThreads>] [-m <nb_MPIProcs>] [-g] [-v]]'

# GET COMMAND-LINE OPTIONS
try:
	options, remainder = getopt.getopt(
		sys.argv[1:],
		'o:m:b:ghvc',
		['OMP=', 'MPI=', 'BENCH=', 'COMPILE_ONLY=', 'GENERATE=', 'HELP=', 'VERBOSE='])
except getopt.GetoptError as err:
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
	elif opt in ('-c', '--COMPILEONLY'):
		COMPILE_ONLY=True
	elif opt in ('-h', '--HELP'):
		print "-b"
		print "     -b input_file"
		print "       Chooses the benchmark(s) to validate."
		print "       input_file : input file(s) to validate. Accepts wildcards."
		print "       input_file=? : prompts a list of all possible input files"
		print "     DEFAULT : All benchmarks are validated."  
		print "-o"
		print "     -o omp_threads"
		print "       omp_threads : number of OpenMP threads used for the execution of smilei (option -e must be present)"
		print "     DEFAULT : 2"  
		print "-m"
		print "     -m mpi_procs"
		print "       mpi_procs : number of MPI processus used for the execution of smilei (option -e must be present)"
		print "     DEFAULT : 2"
		print "-g"
		print "     Generation of references only"
		print "-c"
		print "     Compilation only"
		print "-v"
		print "     Verbose"
		exit()
	elif opt in ('-g', '--GENERATE'):
		GENERATE = True
	elif opt in ('-v', '--VERBOSE'):
		VERBOSE = True

# TEST IF THE NUMBER OF THREADS IS COMPATIBLE WITH THE HOST
if JOLLYJUMPER in HOSTNAME :
	if 12 % OMP != 0:
		print  "Smilei cannot be run with "+str(OMP)+" threads on "+HOSTNAME
		sys.exit(4)  
	NPERSOCKET = 12/OMP

# Build the list of the requested input files
list_bench = [os.path.basename(b) for b in glob.glob(SMILEI_BENCHS+"tst*py")]
if BENCH == "":
	SMILEI_BENCH_LIST = list_bench
elif BENCH == "?":
	VERBOSE = True
	os.chdir(SMILEI_SCRIPTS)
	#- Propose the list of all the input files
	print '\n'.join(list_bench)
	#- Choose an input file name in the list
	print 'Enter an input file from the above list:'
	BENCH = raw_input()
	SMILEI_BENCH_LIST = [ BENCH ]
	while BENCH not in list_bench:
		print "Input file "+BENCH+" invalid. Try again."
		BENCH = raw_input()
		SMILEI_BENCH_LIST = [ BENCH ]
elif BENCH in list_bench:
	SMILEI_BENCH_LIST = [ BENCH ]
elif glob.glob( BENCH ):
	BENCH = glob.glob( BENCH )    
	for b in BENCH:
		if b not in list_bench:
			if VERBOSE:
				print "Input file "+b+" invalid."
			sys.exit(4)
	SMILEI_BENCH_LIST = BENCH
else:
	if VERBOSE:
		print "Input file "+BENCH+" invalid."
	sys.exit(4)

if VERBOSE :
	print ""
	print "The list of input files to be validated is:\n\t"+"\n\t".join(SMILEI_BENCH_LIST)
	print ""

# COMPILE SMILEI

import time
def date(BIN_NAME):
	statbin = os.stat(BIN_NAME)
	return statbin.st_ctime
def date_string(BIN_NAME):
	date_integer = date(BIN_NAME)
	date_time = time.ctime(date_integer)
	return date_time.replace(" ","-")
def workdir_archiv(BIN_NAME) :
	if os.path.exists(SMILEI_W):
		ARCH_WORKDIR = WORKDIR_BASE+'_'+date_string(SMILEI_W)
		os.rename(WORKDIR_BASE,ARCH_WORKDIR)
		os.mkdir(WORKDIR_BASE)

if VERBOSE :
  print "Compiling Smilei"

os.chdir(SMILEI_ROOT)
WORKDIR_BASE = SMILEI_ROOT+"validation"+s+"workdirs"
SMILEI_W = WORKDIR_BASE+s+"smilei"
SMILEI_R = SMILEI_ROOT+s+"smilei"
if os.path.exists(SMILEI_R):
	STAT_SMILEI_R_OLD = os.stat(SMILEI_R)
else :
	STAT_SMILEI_R_OLD = ' '
COMPILE_ERRORS=WORKDIR_BASE+s+'compilation_errors'
COMPILE_OUT=WORKDIR_BASE+s+'compilation_out'

# Find compile command according to the host
if JOLLYJUMPER in HOSTNAME :
	COMPILE_COMMAND = 'unset MODULEPATH;module use /opt/exp_soft/vo.llr.in2p3.fr/modulefiles; module load compilers/icc/16.0.109 mpi/openmpi/1.6.5-ib-icc python/2.7.10 hdf5 compilers/gcc/4.8.2  > /dev/null 2>&1;make -j 2 > compilation_out_temp 2>'+COMPILE_ERRORS  
	CLEAN_COMMAND = 'unset MODULEPATH;module use /opt/exp_soft/vo.llr.in2p3.fr/modulefiles; module load compilers/icc/16.0.109 mpi/openmpi/1.6.5-ib-icc python/2.7.10 hdf5 compilers/gcc/4.8.2 > /dev/null 2>&1;make clean > /dev/null 2>&1'
elif POINCARE in HOSTNAME :
	#COMPILE_COMMAND = 'module load intel/15.0.0 openmpi hdf5/1.8.10_intel_openmpi python gnu > /dev/null 2>&1;make -j 6 > compilation_out_temp 2>'+COMPILE_ERRORS     
	#CLEAN_COMMAND = 'module load intel/15.0.0 openmpi hdf5/1.8.10_intel_openmpi python gnu > /dev/null 2>&1;make clean > /dev/null 2>&1'
	COMPILE_COMMAND = 'module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu ; unset LD_PRELOAD ; export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0 > /dev/null 2>&1;make -j 6 > compilation_out_temp 2>'+COMPILE_ERRORS
	CLEAN_COMMAND = 'module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu ; unset LD_PRELOAD ; export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0 > /dev/null 2>&1;make clean > /dev/null 2>&1'
else:
	COMPILE_COMMAND = 'make -j4 > compilation_out_temp 2>'+COMPILE_ERRORS
	CLEAN_COMMAND = 'make clean > /dev/null 2>&1'


# If the workdir does not contains a smilei bin, or it contains one older than the the smilei bin in directory smilei, force the compilation in order to generate the compilation_output
if not os.path.exists(WORKDIR_BASE):
	os.mkdir(WORKDIR_BASE)
if os.path.exists(SMILEI_R) and (not os.path.exists(SMILEI_W) or date(SMILEI_W)<date(SMILEI_R)):
	call(CLEAN_COMMAND , shell=True)
try :
	if os.path.exists(WORKDIR_BASE+s+COMPILE_ERRORS) :
		os.remove(WORKDIR_BASE+s+COMPILE_ERRORS)
	check_call(COMPILE_COMMAND, shell=True)
	os.rename('compilation_out_temp',COMPILE_OUT)
	if STAT_SMILEI_R_OLD!=os.stat(SMILEI_R) or date(SMILEI_W)<date(SMILEI_R):
		# if new bin, archive the workdir (if it contains a smilei bin)  and create a new one with new smilei and compilation_out inside
		if os.path.exists(SMILEI_W):
			workdir_archiv(SMILEI_W)
		shutil.copy2(SMILEI_R,SMILEI_W)
		if COMPILE_ONLY:
			if VERBOSE:
				print  "Smilei validation succeed."
			exit(0)
	else: 
		if COMPILE_ONLY :
			if VERBOSE:
				print  "Smilei validation not needed."
			exit(0)
except CalledProcessError,e:
	# if compiling errors, archive the workdir (if it contains a smilei bin), create a new one with compilation_errors inside and exit with error code
	workdir_archiv(SMILEI_W)
	os.rename(COMPILE_ERRORS,WORKDIR_BASE+s+COMPILE_ERRORS)
	if VERBOSE:
		print "Smilei validation cannot be done : compilation failed." ,e.returncode
	sys.exit(3)


# DEFINE A CLASS TO CREATE A REFERENCE
class CreateReference(object):
	def __init__(self, bench_name):
		self.reference_file = SMILEI_REFERENCES+s+bench_name+".txt"
		self.data = {}
	
	def __call__(self, data_name, data, precision=None):
		self.data[data_name] = data
	
	def write(self):
		with open(self.reference_file, "w") as f:
			pickle.dump(self.data, f)
		size = os.path.getsize(self.reference_file)
		if size > 1000000:
			print "Reference file is too large ("+str(size)+"B) - suppressing ..."
			os.remove(self.reference_file)
			sys.exit(2)
		if VERBOSE:
			print "Created reference file "+self.reference_file

# DEFINE A CLASS TO COMPARE A SIMULATION TO A REFERENCE
class CompareToReference(object):
	def __init__(self, bench_name):
		try:
			with open(SMILEI_REFERENCES+s+bench_name+".txt", 'r') as f:
				self.data = pickle.load(f)
		except:
			print "Unable to find the reference data for "+bench_name
			sys.exit(1)
	
	def __call__(self, data_name, data, precision=None):
		# verify the name is in the reference
		if data_name not in self.data.keys():
			print "Reference quantity '"+data_name+"' not found"
			sys.exit(1)
		expected_data = self.data[data_name]
		# ok if exactly equal (including strings or lists of strings)
		try   :
			if expected_data == data: return
		except: pass
		# If numbers:
		try:
			double_data = np.double(data)
			if precision is not None:
				error = np.abs( double_data-np.double(expected_data) )
				max_error_index = np.argmax(error)
				max_error = error[max_error_index]
				if max_error < precision: return
				print "Reference quantity '"+data_name+"' does not match the data (required precision "+str(precision)
				print "Max error = "+str(max_error)+" at index "+str(max_error_index)
			else:
				if np.all(double_data == np.double(expected_data)): return
				print "Reference quantity '"+data_name+"' does not match the data"
		except:
			print "Reference quantity '"+data_name+"': unable to compare to data"
		print "Reference data:"
		print expected_data
		print "New data:"
		print data
		sys.exit(1)

# RUN THE BENCHMARKS

for BENCH in SMILEI_BENCH_LIST :
	
	SMILEI_BENCH = SMILEI_BENCHS + BENCH
	
	# CREATE THE WORKDIR CORRESPONDING TO THE INPUT FILE AND GO INTO                
	WORKDIR = WORKDIR_BASE+s+'wd_'+os.path.basename(os.path.splitext(BENCH)[0])
	if not os.path.exists(WORKDIR):
		os.mkdir(WORKDIR)
	os.chdir(WORKDIR)
	
	WORKDIR += s+str(MPI)
	if not os.path.exists(WORKDIR):
		os.mkdir(WORKDIR)
	
	WORKDIR += s+str(OMP)
	if not os.path.exists(WORKDIR):
		os.mkdir(WORKDIR)
		EXECUTION = True
	else:
		EXECUTION = False
	
	os.chdir(WORKDIR)
	
	# define the name of the execution script
	EXEC_SCRIPT = 'exec_script.sh'
	EXEC_SCRIPT_OUT = 'exec_script.out'
	SMILEI_EXE_OUT = 'smilei_exe.out'
	
	# RUN smilei IF EXECUTION IS TRUE
	if EXECUTION :
		if VERBOSE:
			print 'Running '+BENCH+' on '+HOSTNAME+' with '+str(OMP)+'x'+str(MPI)+' OMPxMPI'
		
		#  depending on the host, build the script and run it
		if POINCARE in HOSTNAME:
			with open(EXEC_SCRIPT, 'w') as exec_script_desc:
				print "ON POINCARE NOW"
				exec_script_desc.write(
					"# environnement \n"
					+"module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu 2>&1 > /dev/null\n"
					+"unset LD_PRELOAD\n"
					+"export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0\n"
					+"# \n"
					+"# execution \n"
					+"export OMP_NUM_THREADS="+str(OMP)+"\n"
					+"mpirun -np "+str(MPI)+" "+WORKDIR_BASE+s+"smilei "+SMILEI_BENCH+" >"+SMILEI_EXE_OUT+" \n"
					+"exit $?  "
				)
			
			# RUN SMILEI
			COMMAND = "/bin/bash "+EXEC_SCRIPT+" > "+EXEC_SCRIPT_OUT+" 2>&1"
			try :
			  check_call(COMMAND, shell=True)
			except CalledProcessError,e:
			# if execution fails, exit with exit status 2
				if VERBOSE :
					print  "Smilei validation cannot be done : execution failed."
					os.chdir(WORKDIR)
					COMMAND = "/bin/bash cat "+SMILEI_EXE_OUT
					try :
						check_call(COMMAND, shell=True)
					except CalledProcessError,e:
						print  "cat command failed"
						sys.exit(2)
				os.chdir(WORKDIR_BASE)
				shutil.rmtree(WORKDIR)
				sys.exit(2)
		
		elif JOLLYJUMPER in HOSTNAME :
			with open(EXEC_SCRIPT, 'w') as exec_script_desc:
				NODES=((int(MPI)*int(OMP)-1)/24)+1
				exec_script_desc.write(
					"#PBS -l nodes="+str(NODES)+":ppn=24 \n"
					+"#PBS -q default \n"
					+"#PBS -j oe\n"
					+"#Set the correct modules available \n"
					+"unset MODULEPATH; \n"
					+"module use /opt/exp_soft/vo.llr.in2p3.fr/modulefiles \n"
					+"#Load compilers and mpi \n"
					+"module load compilers/icc/16.0.109 \n"
					+"module load mpi/openmpi/1.6.5-ib-icc \n"
					+"module load python/2.7.10 \n"
					+"module load hdf5 \n"
					+"module load compilers/gcc/4.8.2 \n"
					+" \n"
					+"# -loadbalance to spread the MPI processes among the different nodes. \n"
					+"# -bind-to-core to fix a given MPI process to a fixed set of cores. \n"
					+"# -cpus-per-proc 6 to set said set of cores to a size of 6 (half socket of JJ) which is also the number of omp threads. \n"
					+"export OMP_NUM_THREADS="+str(OMP)+" \n"
					+"export OMP_SCHEDULE=DYNAMIC \n"
					+"export KMP_AFFINITY=verbose \n"
					+"export PATH=$PATH:/opt/exp_soft/vo.llr.in2p3.fr/GALOP/beck \n"
					+"#Specify the number of sockets per node in -mca orte_num_sockets \n"
					+"#Specify the number of cores per sockets in -mca orte_num_cores \n"
					+"cd "+WORKDIR+" \n"
					+"mpirun -mca orte_num_sockets 2 -mca orte_num_cores 12 -cpus-per-proc "+str(OMP)+" --npersocket "+str(NPERSOCKET)+" -n "+str(MPI)
					+" -x $OMP_NUM_THREADS -x $OMP_SCHEDULE "+WORKDIR_BASE+s+"smilei "+SMILEI_BENCH+" >"+SMILEI_EXE_OUT+" 2>&1 \n"
					+"echo $? > exit_status_file \n"
				)
			
			EXIT_STATUS="100"
			exit_status_fd = open(WORKDIR+"/exit_status_file", "w+")
			exit_status_fd.write(str(EXIT_STATUS))
			exit_status_fd.seek(0)
			COMMAND = "PBS_DEFAULT=llrlsi-jj.in2p3.fr qsub  "+EXEC_SCRIPT
			try :
				check_call(COMMAND, shell=True)
			except CalledProcessError,e:
				# if commande qsub fails, exit with exit status 2
				exit_status_fd.close()  
				os.chdir(WORKDIR_BASE)
				shutil.rmtree(WORKDIR)           
				if VERBOSE :
					print  "Smilei validation cannot be done : qsub command failed."
					sys.exit(2)
			
			while ( EXIT_STATUS == "100" ) :
				time.sleep(10)
				EXIT_STATUS = exit_status_fd.readline()
				exit_status_fd.seek(0)
			if ( int(EXIT_STATUS) != 0 )  :
				if VERBOSE :
					print  "Smilei validation cannot be done : execution failed."
					os.chdir(WORKDIR)
					COMMAND = "cat "+SMILEI_EXE_OUT
					try :
						check_call(COMMAND, shell=True)
					except CalledProcessError,e:
						print  "cat command failed"
						sys.exit(2)
				exit_status_fd.close()  
				os.chdir(WORKDIR_BASE)
				#shutil.rmtree(WORKDIR)
				sys.exit(2)
		
		else:
			COMMAND = (
				"export OMP_NUM_THREADS="+str(OMP)+"; "
				+"mpirun -mca btl tcp,sm,self -np "+str(MPI)+" "+WORKDIR_BASE+s+"smilei "+SMILEI_BENCH+" >"+SMILEI_EXE_OUT
			)
			try :
				check_call(COMMAND, shell=True)
			except CalledProcessError,e:
				if VERBOSE :
					print  "Smilei could not run"
				sys.exit(2)
	
	# CHECK THE OUTPUT FOR ERRORS
	errors = []
	search_error = re.compile('error', re.IGNORECASE)
	with open(SMILEI_EXE_OUT,"r") as fout:
		for line in fout:
			if search_error.search(line):
				errors += [line]
	if errors:
		print ""
		print("Errors appeared while running the simulation:")
		print("---------------------------------------------")
		for error in errors:
			print(error)
		sys.exit(2)
	
	# FIND THE VALIDATION SCRIPT FOR THIS BENCH
	validation_script = SMILEI_VALIDATION + "validate_"+BENCH
	print ""
	if not os.path.exists(validation_script):
		print "Unable to find the validation script "+validation_script
		sys.exit(1)
	
	# IF REQUIRED, GENERATE THE REFERENCES
	if GENERATE:
		if VERBOSE:
			print 'Generating reference for '+BENCH
		Validate = CreateReference(BENCH)
		execfile(validation_script)
		Validate.write()
	
	# OTHERWISE, COMPARE TO THE EXISTING REFERENCES
	else:
		if VERBOSE:
			print 'Validating '+BENCH
		Validate = CompareToReference(BENCH)
		execfile(validation_script)
	

print "Everything passed"
