#!/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/bin/python 
#
# This script checks the validity of the results of smilei contained in the file scalars.txt, inside the directory scripts.
#
# One directory named "references" contains one file scalars.txt for each file in benchmarks and a list of precision values for each scalar.
#
# The validity of smilei is obtained by comparing the scalars contained in scalars.txt with those found in "references" with a given precision 
# or a precsion value found in the list.
# Either the existing scalars.txt file is examined or it is produced by a new execution of smilei in the following conditions :
# - -e, or -o, or -m option is present
# - -b option has changed since the last execution of smilei .
# In these 2 cases, smilei is executed with a given number of OpenMP tasks (default 2), 
# a given number of MPI processus (default 2), and the input file given by the -b option 
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
#  - the -s option is all : all the time steps are considered
#  In these 2 last cases the list of precision values for all the possible scalars is examined. 
#
# IMPORTS
import Diagnostics
import getopt
import sys
import os
import shutil
import numpy as np
#
# DEFAULT VALUES FOR OPTIONS
OMP = 2
MPI = 2
SCALAR_LIST_ARG = [ "Utot" ]
SCALAR_NAME = "Utot"
OPT_BENCH = False
OPT_TIMESTEP = False
EXECUTION = False
OPT_PRECISION = False
#
# FUNCTION FOR PARSING OPTIONS
def usage():
    print 'Usage: validation.py [-o <nb_OMPThreads> -m <nb_MPIProcs> -b <bench_case> -s <scalar> -p <precision> -e'
try:
  options, remainder = getopt.getopt(sys.argv[1:], 'o:m:b:s:p:t:he', ['OMP=', 
                                                          'MPI='
                                                          'BENCH='
                                                          'SCALAR='
                                                          'PRECISION='
                                                          'TIMESTEP='
                                                          'HELP='
                                                          'EXECUTION='
                                                         ])
except getopt.GetoptError as err:
        usage()
        sys.exit(2)
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
  shutil.copyfile('../references/sca_'+BENCH, 'scalars.txt') 
  shutil.copyfile('smilei.py', 'smilei.py_save') 
  Sref = Diagnostics.Smilei(".")
  SCALAR_REFERENCE = Sref.Scalar(SCALAR_NAME,timestep=t)
  VALEUR_REFERENCE = SCALAR_REFERENCE._getDataAtTime(t)
  shutil.copyfile('scalars.txt_save','scalars.txt')
  shutil.copyfile('smilei.py_save','smilei.py')
  # 
  # COMPARE VALEUR AGAINST VALEUR_REFERENCE WITH PRECISION PRECISION
  print 'At time',t,',',SCALAR_NAME,'=',VALEUR,', Reference =',VALEUR_REFERENCE
  if float(abs(VALEUR - VALEUR_REFERENCE)) < float(PRECISION) :
    print 'Scalar',SCALAR_NAME,'is OK with precision',PRECISION,'\n'
  else :
    print 'Scalar',SCALAR_NAME,'is wrong according to precision',PRECISION,'\n'
  #
  return 
#
def scalarListValidation(SCALAR_NAME,t) :
  global PRECISION
  # Test if either you want to choose one scalar in a list, or all the scalars, or a given scalar
  # Check these scalars with the function scalarValidation()
  #
  it = np.int64(t)
  L=Diagnostics.Scalar(".",scalar="Utot")
  LISTE_SCALARS = L.getScalars()
  if SCALAR_NAME == "?":
    # Propose the list of all the scalars
    L = Diagnostics.Scalar(".")
    print 'Enter a scalar name from the above list (press <Enter> if no more scalar to check ):'
    SCALAR_NAME = raw_input()
    while SCALAR_NAME != "" :
      if SCALAR_NAME in LISTE_SCALARS :                           
        PRECISION = precision_d[SCALAR_NAME]
        scalarValidation(SCALAR_NAME,it)
        print 'Enter a scalar name from the above list (press <Enter> if no more scalar to check ):'
        SCALAR_NAME = raw_input()
      else :
        print "Scalar", SCALAR_NAME,"is not valid. Enter a scalar name again."
        SCALAR_NAME = raw_input()
        PRECISION = precision_d[SCALAR_NAME]
  elif SCALAR_NAME == "all":
    for SCALAR_NAME in LISTE_SCALARS:
      PRECISION = precision_d[SCALAR_NAME]
      scalarValidation(SCALAR_NAME,it)
  elif len(SCALAR_LIST_ARG) > 1 :
    for SCALAR_NAME in SCALAR_LIST_ARG :
      if SCALAR_NAME in LISTE_SCALARS :                    
        PRECISION = precision_d[SCALAR_NAME]
        scalarValidation(SCALAR_NAME,it)
      else :
        print "Scalar", SCALAR_NAME,"is not valid."
  else:
    if SCALAR_NAME in LISTE_SCALARS :                    
      if not OPT_PRECISION :
        PRECISION = precision_d[SCALAR_NAME]
      scalarValidation(SCALAR_NAME,it)
    else :
      print "Scalar", SCALAR_NAME,"is not valid."
      sys.exit(2)
  return
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
        print "-e"
        print "     If -e is present, smilei will be executed with input file corresponding to -b option"
        print "     DEFAULT : tst1d_0_em_propagation.py"  
        print "-o"
        print "     -o omp_threads"
        print "       omp_threads : number of OpenMP threads used for the execution of smilei (option -e must be present)"
        print "     DEFAULT : 2"  
        print "-m"
        print "     -m mpi_procs"
        print "       mpi_procs : number of MPI processus used for the execution of smilei (option -e must be present)"
        print "     DEFAULT : 2"  

        exit()
    elif opt in ('-e', '--EXECUTION'):
        EXECUTION=True
#
# CASE PRECISION NOT DEFINED OR MULTIPLE SCALARS REQUIRED : 
# CREATING A DICTIONNARY WITH THE LIST OF PRECISION VALUES FOUND IN ./references/precision_values
if OPT_PRECISION and ( SCALAR_NAME == "all"  or  SCALAR_NAME == "?" ) :
  print "\n WARNING : Precision option ignored since a list of scalars is required."
if not OPT_PRECISION or  SCALAR_NAME == "all"  or  SCALAR_NAME == "?"  :
  precision_d = {}
  with open("../references/precision_values") as f:
      for line in f:
         (key, val) = line.split()
         precision_d[key] = val
#
# READ THE INPUT FILE
if OPT_BENCH == False :
  # IF BENCH NOT IN OPTIONS, SET THE DEFAULT
  BENCH = "tst1d_0_em_propagation.py"
elif BENCH == "?":
  # Build the list of the input files
  with open('list_bench.sh', 'w') as list_bench_desc :
    list_bench_desc.write( "cd ../benchmarks; ls tst*py \n")
  os.system("/bin/sh /gpfshome/mds/staff/mfle/GITLAB/smilei/scripts/list_bench.sh > /gpfshome/mds/staff/mfle/GITLAB/smilei/scripts/list_bench")
  with open('list_bench', 'r') as list_bench_desc :
    list_bench = list_bench_desc.read()
  # Propose the list of all the input files
  print list_bench           
  # Choose an input file name in the list
  print 'Enter an input file from the above list:'
  BENCH = raw_input()
  while not  BENCH in list_bench:
    print "Input file",BENCH,"not valid. Enter an input file name again."
    BENCH = raw_input()
#
# PRINT INFO           
print '\nTesting smilei with', OMP, 'OMP tasks,',MPI, 'MPI processes, case',BENCH,':'
#
# EXECUTE smilei IF -e OPTON IS SET OR INPUT FILE HAS CHANGED SINCE LAST EXECUTION : 
#    
if not EXECUTION :
#   OPEN smilei EXECUTION SCRIPT                   
  job_validation = open('validation.ll', 'rw')
#
#   TEST IF smilei MUST BE EXECUTED (-b option argument and input file name in validation.ll are different)
  line = ' '
  EXECUTION=True
  while EXECUTION and line != '' :
    line=job_validation.read()
    if BENCH in line:
      EXECUTION=False
  job_validation.close()
#
#  EXECUTE smilei IF EXECUTION IS TRUE
if EXECUTION : 
  # BUILD THE SMILEI EXECUTION SCRIPT                   
  job_validation = open('validation.ll', 'w')
  print "Executing smilei, please wait ..."
  job_validation.write( "\
# environnement \n \
module load intel/15.0.0 openmpi ddt hdf5/1.8.10_intel_openmpi python\n  \
# \n \
# compilation\n  \
cd ..\n  \
#make  openmpgnu \n  \
make  \n \
# execution \n \
cd  scripts \n \
export OMP_NUM_THREADS="+str(OMP)+"\n \
mpirun -bind-to-socket -np "+str(MPI)+" ../smilei ../benchmarks/"+BENCH+" >out_validation.ll\n")
#
  job_validation.close()
#
# RUN  SMILEI
  os.system("/bin/sh /gpfshome/mds/staff/mfle/GITLAB/smilei/scripts/validation.ll > /dev/null 2>&1")
#
# READ TIME STEPS AND VALIDATE 
L=Diagnostics.Scalar(".",scalar="Utot") # scalar argument must be anything except none in order times is defined
LISTE_TIMESTEPS = L.getAvailableTimesteps()
if OPT_TIMESTEP == False :
  # IF TIMESTEP NOT IN OPTIONS, COMPUTE THE LAST TIME STEP               
  TIMESTEP=LISTE_TIMESTEPS[-1]
  scalarListValidation(SCALAR_NAME,TIMESTEP)
elif TIMESTEP == "?":
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
      scalarListValidation(SCALAR_NAME,TIMESTEP) 
    else :
      print "Time step :", TIMESTEP,"is not valid. Enter a time step name again."
      TIMESTEP_INPUT = raw_input()
    print "Enter a time step name again (press <Enter> if no more time step to check ):"
    TIMESTEP_INPUT = raw_input()
elif TIMESTEP == "all":
  for TIMESTEP in LISTE_TIMESTEPS:
    scalarListValidation(SCALAR_NAME,TIMESTEP)
else:
    scalarListValidation(SCALAR_NAME,TIMESTEP)
