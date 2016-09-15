#! /usr/bin/python
##!/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/bin/python 
#
# This script checks the validity of the results of the execution of smilei defined by :
# - the bench file (default tst1d_0_em_propagation.py) given by the -b option
# - the number of OpenMP tasks (default 2) given by the -o option
# - the number of MPI processus (default 2) given by the -m option
# 
# Here are listed the files used in the validation process:
###########################################################
#- validation directory which contains:
#  - a directory named "references" which contains one file scalars.txt for each file in benchmarks and a list of precision values for each scalar.
#  - a workdirs directory for the current version of smilei.
#  - archived workdirs for previous versions of smilei
# Before the first execution of the script, workdirs is empty.
# Then, workdirs will contain :
#  - smilei binary : "smilei"
#  - compilation output : "compilation_output"
#  - a directory wd_<input_file> which contains one <o>_<m> directory 
#    The <o>_<m> directory contains the output files of the execution of smilei with o OpenMP threads and m MPI processus :
#    - exec_script.out 
#    - exec_script.sh  
#    - Fields.h5  
#    - scalars.txt  
#    - smilei_exe.out  
#    - smilei.py
#  or
#  - the compilation errors file : "compilation_errors"
#
# The different steps of the script are the folllowing :
########################################################
# Compilation step :
# +++++++++++++++
# If the "workdirs" directory does not contains a smilei binary, or it contains one older than the the smilei bin in directory smilei, 
# then the smilei binary in directory smilei is removed in order to force the compilation and so generate the compilation_output. 
# Then the compilation occurs.
# If a new smilei binary is created, then :
# if the "workdirs" directory contains a smilei bin, then it is archived and a new one is created with the new smilei binary 
# and the compilation output inside.
# If compiling errors occur, the workdir (if it contains a smilei bin) is archived and a new one is created with compilation_errors inside 
# and the script exits with status 3.
#
# Execution step :
# +++++++++++++++
# For the given execution defined by -b, -m, -o options, if the corresponding directory named wd_<bench_file>_<o>_<m> does not exist then :
# - it is created.
# - smilei is executed in this working directory. So, the different output files such as scalars.txt and the standard ouput of this execution are kept 
#   in this working directory.
# - if the execution fails then the script exits with status 2.
#
# Validation step :
# +++++++++++++++
# The validity of smilei is obtained by comparing the scalars contained in scalars.txt contained in the working directory "wd_<bench_file>/<o>_<m>" 
# with those found in "references" with a given precision or a precision value found in the list.
#
# Using the -a option, smilei is executed for all the bench files. For every execution all the existing scalars in the scalars.txt file are checked.
#
# Using the -s option, you have the choice to check several scalars with their corresponding precision:
#  - no -s option : scalar Utot is checked with the given precision (-p option, default found in the list of precision values))
#  - the -s option is a scalar name: this scalar is  checked with the given precision (-p option, default found in the list of precision values)
#  - the -s option is ? : several scalars that you can choose inside a provided list 
#  - the -s option is all : all the existing scalars inside the file scalars.txt 
#  In these 2 last cases the list of precision values for all the possible scalars is examined. 
#
# In the same way, -t option allows to validate smilei for one or more timesteps :         
#  - no -t option : the last time step is considered                                                                                   
#  - the -t option is a number : smilei is checked only for this time step                                                                               
#  - the -t option is ? : time steps that you can choose inside a provided list are considered
#  - the -t option is all : all the time steps are considered
#  In these 2 last cases the list of precision values for all the possible scalars is examined. 
#
# Smilei is validated (exit status 0) only in the case where all the scalars are considered (-s all option) 
# and validated at last time step (no -t option). Otherwise, exit status is one.
#
# Output of the script
#+++++++++++++++++++++
# If the -v option is present, the script gives information about what it does and lists the validity of the scalars
#
# Exit status of the script
#++++++++++++++++++++++++++
# 0  validation : all the scalars are validated for the last time step (this may happen only with the -s all option)                       
# 1  no validation : not all the scalars are checked or one is wrong
# 2  execution fails
# 3  compilation fails
# 4  bad option
#
# Remark
#++++++++
# This script may run anaywhere : you can define a SMILEI_ROOT variable which is a path containing a smilei version, 
# otherwise, default SMILEI_ROOT variable is the current directory (which must contain a smilei version).
#
# IMPORTS
import Diagnostics
import getopt
import shutil
import numpy as np
import glob
from subprocess import check_call,CalledProcessError,call
import os
import sys
import socket
#
# SMILEI PATH VARIABLES
if "SMILEI_ROOT" in os.environ :
  SMILEI_ROOT=os.environ["SMILEI_ROOT"]
else :
  SMILEI_ROOT = os.getcwd()+'/..'
SMILEI_ROOT = SMILEI_ROOT+'/'
SMILEI_SCRIPTS = SMILEI_ROOT+"scripts/"
SMILEI_REFERENCES = SMILEI_ROOT+"validation/references/"
SMILEI_BENCHS = SMILEI_ROOT+"benchmarks/"
sys.path.insert(0, SMILEI_SCRIPTS)
#
# OTHER VARIABLES
POINCARE="poincare"
JOLLYJUMPER="llrlsi-gw"
HOSTNAME=socket.gethostname()
# DEFAULT VALUES FOR OPTIONS
OMP = 2
MPI = 2
SCALAR_LIST_ARG = [ "Utot" ]
SCALAR_NAME = "Utot"
OPT_BENCH = False
OPT_TIMESTEP = False
EXECUTION = False
OPT_PRECISION = False
VERBOSE = False
VALID_ALL = False
BENCH=""
#
# FUNCTION FOR PARSING OPTIONS
def usage():
    print 'Usage: validation.py [-o <nb_OMPThreads> -m <nb_MPIProcs> -b <bench_case> -s <scalar> -p <precision> '
try:
  options, remainder = getopt.getopt(sys.argv[1:], 'o:m:b:s:p:t:hva', ['OMP=', 
                                                          'MPI='
                                                          'BENCH='
                                                          'SCALAR='
                                                          'PRECISION='
                                                          'TIMESTEP='
                                                          'HELP='
                                                          'VERBOSE='
                                                          'ALL='
                                                         ])
except getopt.GetoptError as err:
        usage()
        sys.exit(4)
#
def scalarValidation(SCALAR_NAME,t):
  # COMPARES SCALAR_NAME VALUE AGAINST THE REFERENCE VALUE AT TIME t :
  #
  # FIND THE VALUE OF SCALAR_NAME in  scalars.txt
  S = Diagnostics.Smilei(".")
  SCALAR = S.Scalar(SCALAR_NAME,timestep=t)
  VALEUR = SCALAR._getDataAtTime(t)
  #
  # FIND THE REFERENCE VALUE OF SCALAR_NAME                          
  shutil.copyfile('scalars.txt','scalars.txt_save')
  shutil.copyfile(SMILEI_REFERENCES+'/sca_'+BENCH, 'scalars.txt') 
  shutil.copyfile('smilei.py', 'smilei.py_save') 
  Sref = Diagnostics.Smilei(".")
  SCALAR_REFERENCE = Sref.Scalar(SCALAR_NAME,timestep=t)
  VALEUR_REFERENCE = SCALAR_REFERENCE._getDataAtTime(t)
  os.rename('scalars.txt_save','scalars.txt')
  os.rename('smilei.py_save','smilei.py')
  # 
  # COMPARE VALEUR AGAINST VALEUR_REFERENCE WITH PRECISION PRECISION
  if VERBOSE :
    print 'At time',t,',',SCALAR_NAME,'=',VALEUR,', Reference =',VALEUR_REFERENCE
  if float(abs(VALEUR - VALEUR_REFERENCE)) < float(PRECISION) :
    if VERBOSE :
      print 'Scalar',SCALAR_NAME,'is OK with precision',PRECISION,'\n'
    SCALAR_OK = True
  else :
    if VERBOSE :
      print 'Scalar',SCALAR_NAME,'is wrong according to precision',PRECISION,'\n'
    SCALAR_OK = False
  #
  return(SCALAR_OK) 
#
def scalarListValidation(SCALAR_NAME,t) :
  global PRECISION
  global VERBOSE
  # Test if either you want to choose one scalar in a list, or all the scalars, or a given scalar
  # Check these scalars with the function scalarValidation()
  #
  it = np.int64(t)
  L=Diagnostics.Scalar(".",scalar="Utot")
  LISTE_SCALARS = L.getScalars()
  if SCALAR_NAME == "?":
    VERBOSE = True
    # Propose the list of all the scalars
    L = Diagnostics.Scalar(".")
    print '\nEnter a scalar name from the above list (press <Enter> if no more scalar to check ):'
    SCALAR_NAME = raw_input()
    while SCALAR_NAME != "" :
      if SCALAR_NAME in LISTE_SCALARS :                           
        PRECISION = precision_d[SCALAR_NAME]
        VALIDATION = scalarValidation(SCALAR_NAME,it)
        print 'Enter a scalar name from the above list (press <Enter> if no more scalar to check ):'
        SCALAR_NAME = raw_input()
      else :
        print "Scalar", SCALAR_NAME,"is not valid. Enter a scalar name again."
        SCALAR_NAME = raw_input()
        PRECISION = precision_d[SCALAR_NAME]
  elif SCALAR_NAME == "all" or VALID_ALL :
    for SCALAR_NAME in LISTE_SCALARS:
      PRECISION = precision_d[SCALAR_NAME]
      VALIDATION = scalarValidation(SCALAR_NAME,it)
      if not VALIDATION :
        return(False)
    return(True)           
  elif len(SCALAR_LIST_ARG) > 1 :
    for SCALAR_NAME in SCALAR_LIST_ARG :
      if SCALAR_NAME in LISTE_SCALARS :                    
        PRECISION = precision_d[SCALAR_NAME]
        VALIDATION = scalarValidation(SCALAR_NAME,it)
      else :
        print "Scalar", SCALAR_NAME,"is not valid."
    return(False)
  else:
    if SCALAR_NAME in LISTE_SCALARS :                    
      if not OPT_PRECISION :
        PRECISION = precision_d[SCALAR_NAME]
      VALIDATION = scalarValidation(SCALAR_NAME,it)
    else :
      print "Scalar", SCALAR_NAME,"is not valid."
  return(False)
#
# PROCESS THE OPTIONS
for opt, arg in options:
    if opt in ('-o', '--OMP'):
        EXECUTION = True
        OMP = arg
    elif opt in ('-m', '--MPI'):
        EXECUTION = True
        MPI = arg
    elif opt in ('-b', '--BENCH'):
        BENCH=arg
        OPT_BENCH = True
    elif opt in ('-s', '--SCALAR'):
        SCALAR_NAME = arg
        SCALAR_LIST_ARG = arg.split()
        if len(SCALAR_LIST_ARG) == 1 : 
          SCALAR_NAME = arg
    elif opt in ('-p', '--PRECISION'):
        PRECISION = arg
        OPT_PRECISION = True
    elif opt in ('-t', '--TIMESTEP'):
        TIMESTEP = arg
        OPT_TIMESTEP=True
    elif opt in ('-h', '--HELP'):
        print "-s"
        print "     -s scalar_name"
        print "        scalar_name=? : a list of all possible scalars is provided : choose a name in this list and press enter when no more scalar to check"
        print "        scalar_name=all : all the existing scalars inside the file scalars.txt are checked"
        print "        scalar_name=<\"scalar_name1, scalar_name2, ...\"> : scalars scalar_name1, scalar_name2, ... are checked"
        print "        scalar_name=<scalar_name> : only <scalar_name> is checked"
        print "     DEFAULT : Utot\n"
        print "-t"
        print "     -t time_step"
        print "       time_step : time step at which scalars are checked"
        print "        time_step=? : a list of all possible time_steps is provided : choose a name in this list and press enter when no more time step to consider"
        print "        time_step=all : all the existing time steps inside the file scalars.txt are considered"
        print "        time_step=<time_step> : only <time_step> is considered"
        print "     DEFAULT : last time step\n"
        print "-p"
        print "     -p precision"
        print "       precision : precision with which <scalar_name> is checked"
        print "       This option is used only when scalar_name is a single scalar name"
        print "     DEFAULT : the precision of this scalar name find in the file references/precision_values\n"
        print "-b"
        print "     -b input_file"
        print "       input_file : input file corresponding to the execution being validated"
        print "       input_file=? : a list of all possible input file names is provided : choose a name in this list and press enter when no more input file to process"
        print "     DEFAULT : tst1d_0_em_propagation.py"  
        print "-o"
        print "     -o omp_threads"
        print "       omp_threads : number of OpenMP threads used for the execution of smilei (option -e must be present)"
        print "     DEFAULT : 2"  
        print "-m"
        print "     -m mpi_procs"
        print "       mpi_procs : number of MPI processus used for the execution of smilei (option -e must be present)"
        print "     DEFAULT : 2"  
        print "-v"
        print "     This option allows to print the messages on standard output"
        exit()
    elif opt in ('-v', '--VERBOSE'):
        VERBOSE = True
    elif opt in ('-a', '--ALL'):
        VALID_ALL = True
#
# TEST IF THE NUMBER OF THREADS IS COMPATIBLE WITH THE HOST
if JOLLYJUMPER in HOSTNAME :
  if (12 % OMP != 0) :
    print  "Smilei cannot be run with " ,OMP ," threads on ", HOSTNAME
    sys.exit(4)  
  NPERSOCKET=12/OMP
#
# CASE PRECISION NOT DEFINED OR MULTIPLE SCALARS REQUIRED : 
# CREATING A DICTIONNARY WITH THE LIST OF PRECISION VALUES FOUND IN ./references/precision_values
if OPT_PRECISION and ( SCALAR_NAME == "all"  or  SCALAR_NAME == "?" or VALID_ALL ) :
  print "\n WARNING : Precision option ignored since a list of scalars is required."
if not OPT_PRECISION or  SCALAR_NAME == "all"  or  SCALAR_NAME == "?"  or VALID_ALL :
  precision_d = {}
  with open(SMILEI_REFERENCES+"/precision_values") as f:
      for line in f:
         (key, val) = line.split()
         precision_d[key] = val
#
# PROCESS THE INPUT FILE
#
  #-  Build the list of the input files
os.chdir(SMILEI_BENCHS)
list_bench = glob.glob("tst*py")
if VALID_ALL or BENCH == "all":
  SMILEI_BENCH_LIST = list_bench 
elif OPT_BENCH == False :
  #-  IF BENCH NOT IN OPTIONS, SET THE DEFAULT
  BENCH = "tst1d_0_em_propagation.py"
  SMILEI_BENCH_LIST = [ BENCH ]
elif BENCH == "?":
  VERBOSE = True
  os.chdir(SMILEI_SCRIPTS)
  #- Propose the list of all the input files
  print '\n'.join(list_bench)
  #- Choose an input file name in the list
  print 'Enter an input file from the above list:'
  BENCH = raw_input()
  SMILEI_BENCH_LIST = [ BENCH ]
  while not  BENCH in list_bench:
    print "Input file",BENCH,"not valid. Enter an input file name again."
    BENCH = raw_input()
    SMILEI_BENCH_LIST = [ BENCH ]
elif BENCH in list_bench:
  SMILEI_BENCH_LIST = [ BENCH ]
else :
  if VERBOSE :
    print "Input file",BENCH,"not valid."
  sys.exit(4)  
#
# COMPILE  SMILEI
import time
def date(BIN_NAME):
  statbin=os.stat(BIN_NAME)
  return statbin.st_ctime
def date_string(BIN_NAME):
  date_integer=date(BIN_NAME)
  date_time=time.ctime(date_integer)
  return date_time.replace(" ","-")
def workdir_archiv(BIN_NAME) :
  if os.path.exists(SMILEI_W):
    ARCH_WORKDIR = WORKDIRS+'_'+date_string(SMILEI_W)
    os.rename(WORKDIRS,ARCH_WORKDIR)
    os.mkdir(WORKDIRS) 
  return 
if VERBOSE :
  print "\nRunning make command."
os.chdir(SMILEI_ROOT)
WORKDIRS = SMILEI_ROOT+"validation/workdirs"
SMILEI_W=WORKDIRS+"/smilei"
SMILEI_R=SMILEI_ROOT+"smilei"
if os.path.exists(SMILEI_R):
  STAT_SMILEI_R_OLD=os.stat(SMILEI_R)
else :
  STAT_SMILEI_R_OLD = ' '
COMPILE_ERRORS='compilation_errors'
COMPILE_OUT='compilation_out'
COMPILE_COMMAND = 'module load intel/15.0.0 openmpi hdf5/1.8.10_intel_openmpi python > /dev/null 2>&1;make -j 6 > compilation_out_temp 2>'+COMPILE_ERRORS     
CLEAN_COMMAND = 'module load intel/15.0.0 openmpi hdf5/1.8.10_intel_openmpi python > /dev/null 2>&1;make clean > /dev/null 2>&1'
# If the workdir does not contains a smilei bin, or it contains one older than the the smilei bin in directory smilei, force the compilation in order to generate the compilation_output
#import pdb;pdb.set_trace()
if not os.path.exists(WORKDIRS) :
  os.mkdir(WORKDIRS)
if ( os.path.exists(SMILEI_R) and not os.path.exists(SMILEI_W)) or (os.path.exists(SMILEI_R) and date(SMILEI_W) < date(SMILEI_R)) :
  call(CLEAN_COMMAND , shell=True)
try :
  if os.path.exists(WORKDIRS+'/'+COMPILE_ERRORS) :
    os.remove(WORKDIRS+'/'+COMPILE_ERRORS)
  check_call(COMPILE_COMMAND, shell=True)
  if STAT_SMILEI_R_OLD != os.stat(SMILEI_R) or date(SMILEI_W) < date(SMILEI_R):
# if new bin, archive the workdir (if it contains a smilei bin)  and create a new one with new smilei and compilation_out inside
    if os.path.exists(SMILEI_W):
      workdir_archiv(SMILEI_W)
    shutil.copy2(SMILEI_R,SMILEI_W)
    if STAT_SMILEI_R_OLD != os.stat(SMILEI_R):
      os.rename('compilation_out_temp',WORKDIRS+'/'+COMPILE_OUT)
except CalledProcessError,e:
# if compiling errors, archive the workdir (if it contains a smilei bin), create a new one with compilation_errors inside and exit with error code
  workdir_archiv(SMILEI_W)
  os.rename(COMPILE_ERRORS,WORKDIRS+'/'+COMPILE_ERRORS)
  if VERBOSE :
    print  "Smilei validation cannot be done : compilation failed." ,e.returncode
  sys.exit(3)
#
#import pdb;pdb.set_trace()
for BENCH in SMILEI_BENCH_LIST :
    # PRINT INFO           
    if VERBOSE :
     print '\nTesting smilei on', HOSTNAME,'with', OMP, 'OMP tasks,',MPI, 'MPI processes, case',BENCH,':\n'
    #
    SMILEI_BENCH = SMILEI_BENCHS+BENCH
    VALID_ALL_OK = False
    # CREATE THE WORKDIR CORRESPONDING TO THE INPUT FILE AND GO INTO                
    WORKDIR = WORKDIRS+'/wd_'+BENCH
    if not os.path.exists(WORKDIR) :
      os.mkdir(WORKDIR)
    os.chdir(WORKDIR)
    WORKDIR = WORKDIR+"/"+str(MPI)+"_"+str(OMP)
    if not os.path.exists(WORKDIR) :
      os.mkdir(WORKDIR)
      os.chdir(WORKDIR)
      EXECUTION = True
    else:
      os.chdir(WORKDIR)
      EXECUTION = False
    #
    #  RUN smilei IF EXECUTION IS TRUE
    if EXECUTION :
      # RUN SMILEI                   
        #  define the name of the execution script
      EXEC_SCRIPT = 'exec_script.sh'
      EXEC_SCRIPT_OUT = 'exec_script.out'
      SMILEI_EXE_OUT = 'smilei_exe.out'
      exec_script_desc = open(EXEC_SCRIPT, 'w')
      if VERBOSE :
        print "Executing smilei, please wait ...\n"
        #  depending on the host, build the script and run it
      if POINCARE in HOSTNAME :
        print "ON EST SUR POINCARE"
#        exec_script_desc.write( "\n")
        exec_script_desc.write( "\
    # environnement \n \
module load intel/15.0.0 openmpi  hdf5/1.8.10_intel_openmpi python\n  \
# \n \
# execution \n \
export OMP_NUM_THREADS="+str(OMP)+"\n \
mpirun -bind-to-socket -np "+str(MPI)+" "+WORKDIRS+"/smilei "+SMILEI_BENCH+" >"+SMILEI_EXE_OUT+" \n exit $?  ")
        exec_script_desc.close()
    #
        # RUN  SMILEI
        COMMANDE = "/bin/bash "+EXEC_SCRIPT+" > "+EXEC_SCRIPT_OUT+" 2>&1"
        try :
          check_call(COMMANDE, shell=True)
        except CalledProcessError,e:
        # if execution fails, exit with exit status 2
          os.chdir(WORKDIRS)
          shutil.rmtree(WORKDIR)
          if VERBOSE :
            print  "Smilei validation cannot be done : execution failed."
          sys.exit(2)
      elif JOLLYJUMPER in HOSTNAME :
        NODES=((int(MPI)*int(OMP)-1)/12)+1
        exec_script_desc.write( "\
#PBS -l nodes="+str(NODES)+":ppn=24 \n \
#PBS -q default \n \
#PBS -j oe\n \
#Set the correct modules available \n \
unset MODULEPATH; \n \
module use /opt/exp_soft/vo.llr.in2p3.fr/modulefiles \n \
#Load compilers and mpi \n \
module load compilers/icc/16.0.109 \n \
module load mpi/openmpi/1.6.5-ib-icc \n \
module load python/2.7.10 \n \
module load hdf5 \n \
 \n \
# -loadbalance to spread the MPI processes among the different nodes. \n \
# -bind-to-core to fix a given MPI process to a fixed set of cores. \n \
# -cpus-per-proc 6 to set said set of cores to a size of 6 (half socket of JJ) which is also the number of omp threads. \n \
export OMP_NUM_THREADS="+str(OMP)+" \n \
export OMP_SCHEDULE=DYNAMIC \n \
export KMP_AFFINITY=verbose \n \
export PATH=$PATH:/opt/exp_soft/vo.llr.in2p3.fr/GALOP/beck \n \
#Specify the number of sockets per node in -mca orte_num_sockets \n \
#Specify the number of cores per sockets in -mca orte_num_cores \n \
cd "+SMILEI_SCRIPTS+" \n \
mpirun -mca orte_num_sockets 2 -mca orte_num_cores 12 -cpus-per-proc "+str(OMP)+" --npersocket "+str(NPERSOCKET)+" -n "+str(MPI)+"\
 -x $OMP_NUM_THREADS -x $OMP_SCHEDULE "+WORKDIRS+"/smilei "+SMILEI_BENCH+" >"+SMILEI_EXE_OUT+"2>&1 \n \
exit $?  \n  ")
        exec_script_desc.close()
        import pdb;pdb.set_trace()
#    #
#    # RUN  SMILEI
#      COMMANDE = "/bin/qsub -sync y  "+EXEC_SCRIPT
#      try :
#        check_call(COMMANDE, shell=True)
#      except CalledProcessError,e:
#        # if execution fails, exit with exit status 2
#        os.chdir(WORKDIRS)
#        shutil.rmtree(WORKDIR)
#        if VERBOSE :
#          print  "Smilei validation cannot be done : execution failed."
#        sys.exit(2)
#
    # READ TIME STEPS AND VALIDATE 
    if VERBOSE :
      print "Testing scalars :\n"
    VALID_OK = False
#    import pdb;pdb.set_trace()
    L=Diagnostics.Scalar(".",scalar="Utot") # scalar argument must be anything except none in order times is defined
    LISTE_TIMESTEPS = L.getAvailableTimesteps()
    if OPT_TIMESTEP == False or VALID_ALL:
      # IF TIMESTEP NOT IN OPTIONS, COMPUTE THE LAST TIME STEP               
      TIMESTEP=LISTE_TIMESTEPS[-1]
      VALID_OK = scalarListValidation(SCALAR_NAME,TIMESTEP)
    elif TIMESTEP == "?":
      VERBOSE = True
      # Test if either you want to choose one time step in a list, or all the time steps, or a given time step
      # Validate smilei for a list of scalars (scalarListValidation function)
      #
      # Propose the list of all the timesteps
      print LISTE_TIMESTEPS
      print 'Enter a time step from the above list (press <Enter> if no more time step to check ):'
      # Test if either you want to choose one scalar in a list, or all the scalars, or a given scalar
      TIMESTEP_INPUT = raw_input()
      while TIMESTEP_INPUT != "" :
        TIMESTEP=TIMESTEP_INPUT     
        if int(TIMESTEP) in LISTE_TIMESTEPS :
          VALID = scalarListValidation(SCALAR_NAME,TIMESTEP) 
        else :
          print "Time step :", TIMESTEP,"is not valid. Enter a time step name again."
          TIMESTEP_INPUT = raw_input()
        print "Enter a time step name again (press <Enter> if no more time step to check ):"
        TIMESTEP_INPUT = raw_input()
    elif TIMESTEP == "all":
      for TIMESTEP in LISTE_TIMESTEPS:
        VALID = scalarListValidation(SCALAR_NAME,TIMESTEP)
    else:
        VALID = scalarListValidation(SCALAR_NAME,TIMESTEP)
    if VALID_OK :
      if VERBOSE :
        print "Smilei is valid for the list of scalars at last timestep, input file "+BENCH+", "+str(OMP)+" OpenMP tasks and "+str(MPI)+" MPI processus"
      if VALID_ALL :
        VALID_ALL_OK = True
    else :
      exit(1)
if VALID_ALL and VALID_ALL_OK  :
    exit(0)

