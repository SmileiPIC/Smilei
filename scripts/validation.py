#!/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/bin/python
# This script checks the validity of the results of smilei.
# wirh the given precision (default 0.5).
# A directory named reference has been created and contains one file scalars.txt for each file in benchmarks.
# Smilei is executed with a given number of OpenMP tasks (default 2), a given number of MPI processus (default 2), 
# and a given test case (default : tst1d_0_em_propagation.py) in the benchmark directory.
# Then the field  Utot is compared with the reference value, at a given time step (default : the last), 
# wirh the given precision (default 0.5).
#
# IMPORTS
import Diagnostics
import getopt
import sys
import os
import shutil
#
# COMPUTE THE LAST TIME STEP               
import smilei
T=smilei.timestep
ST=smilei.sim_time
t=int(ST/T)

# ARGUMENTS - OPTIONS
OMP = 2
MPI = 2
BENCH = "tst1d_0_em_propagation.py"
NOM_FIELD="Utot"
PRECISION = 0.5

def usage():
    print 'Usage: validation.py [-o <nb_OMPThreads> -m <nb_MPIProcs> -b <bench_case> -f <Field> -p <precision>'
try:
  options, remainder = getopt.getopt(sys.argv[1:], 'o:m:b:f:p:', ['OMP=', 
                                                          'MPI='
                                                          'BENCH='
                                                          'FIELD='
                                                          'PRECISION='
                                                         ])
except getopt.GetoptError as err:
        # print help information and exit:
#        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

for opt, arg in options:
    if opt in ('-o', '--OMP'):
        OMP = arg
    elif opt in ('-m', '--MPI'):
        MPI = arg
    elif opt in ('-b', '--BENCH'):
        BENCH = arg
    elif opt in ('-f', '--FIELD'):
        NOM_FIELD = arg
    elif opt in ('-p', '--PRECISION'):
        PRECISION = arg
   

# BUILD THE SMILEI EXECUTION SCRIPT                   
job_validation = open('validation.ll', 'w')
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
mpirun -bind-to-socket -np "+str(MPI)+" ../smilei ../benchmarks/"+BENCH+" \n")
#
job_validation.close()
#
# RUN  SMILEI
os.system("validation.ll > /dev/null 2>&1")
#

# FIND THE VALUE OF THE FIELD UTOT in  scalars.txt
S = Diagnostics.Smilei(".")
Diag = S.Scalar(NOM_FIELD)
FIELD = S.Scalar(NOM_FIELD,timestep=t)
VALEUR = FIELD.getData()[t]
#
#
# FIND THE REFERENCE VALUE OF UTOT                          
shutil.copyfile('../references/sca_'+BENCH, '../references/scalars.txt') 
shutil.copyfile('smilei.py', '../references/smilei.py') 
Sref = Diagnostics.Smilei("../references")
Diag = Sref.Scalar(NOM_FIELD)
FIELD = Sref.Scalar(NOM_FIELD,timestep=t)
VALEUR_REFERENCE = FIELD.getData()[t]
# 
# COMPARER VALEUR et VALEUR_REFERENCE A precision pres
if abs(VALEUR - VALEUR_REFERENCE) < PRECISION :
  print 'Testing smilei with', OMP, 'OMP tasks,',MPI, 'MPI processes,',\
  'case ',BENCH,' :\n',NOM_FIELD,' OK with precision :', PRECISION
else :
  print 'Testing smilei with : \n', OMP, 'OMP tasks\n',MPI, 'MPI processes\n',\
  'Case ',BENCH,'\nNOM_FIELD bad value :', VALEUR, 'instead of :',VALEUR_REFERENCE
print 'VALEUR',VALEUR

