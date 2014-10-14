# Usage : ./validation.sh -i input -e error -p nprocs -t nthreads
# input : input file, located in tests, tst3_thermal_expansion.in by default
# error : error tolerance, 0.001 by default
# reference : reference file, located in scripts 
# nprocs : number of MPI processus , 2 by default
# nthreads : number of OpenMP thrteads , 2 by default
#
# Object : this script tests local modifications by comparing energy computed at the last iteration with the one contained in a reference file; 
# this file is either provided by the -r option 
# or the file scripts/reference.<input>.<branch>,
# computed by version v1.2rc-149-g9d0d56b of smilei with 1 MPI proc, 1 OpenMP task and -O0 option ,
# and where branch is the current git branch
#
# Steps of the script :
# - compilation and execution of the code to validate  with  the input file, the  number of MPI procs and the number of MPI tasks given as parameters
# - if the difference between the energy computed at the last iteration and the reference exceeds the error value 
# then the script ends with the exit status 1.
#
# Marie Fle IDRIS
#
# environnement
module load intel openmpi ddt hdf5/1.8.10_intel_openmpi
#ROOT="/gpfshome/mds/staff/mfle/SMILEI/smilei"
ROOT=`cd ..; pwd `
while getopts "i:e:r:p:t:" OPTION
do
        case $OPTION in
        i )    INPUT=$OPTARG;;
        e )    ERROR=$OPTARG;;
        r )    REFERENCEFILE=${ROOT}/scripts/$OPTARG;;
        p )    NPROCS=$OPTARG;;
        t )    OMP_NUM_THREADS=$OPTARG;;
        esac
done
INPUT=${INPUT:-tst3_thermal_expansion.in}
INPUT_PATH=${ROOT}/tests/${INPUT}
ERROR=${ERROR:-0.01}
OUTPUT_VALIDATION=/tmp/output_validation
BRANCH=`git status | head -1 |grep -i branch|awk '{print $4}'`
REFERENCEFILE=${REFERENCEFILE:-${ROOT}/scripts/reference_${INPUT}.${BRANCH}}
#
# Search number of iterations in input file
NB_IT=`grep "res_time =" ${INPUT_PATH}| awk '{print $3}'`
# 
# Inside the reference file, find the energy which corresponds to the number of iterations
REFERENCE=`grep "it = ${NB_IT}" ${REFERENCEFILE}| awk '{ print $12 }'`
# 
# Compilation
cd ${ROOT}/src
make openmpintel
# 
# Execution de smilei avec 2 MPI, 2 OpenMP
export NPROCS=${NPROCS:-2}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
export OMP_SCHEDULE="static"
mpirun  -np ${NPROCS}  ${ROOT}/src/smilei ${INPUT_PATH} > ${OUTPUT_VALIDATION}
#
# Inside the output file, find the energy which corresponds to the number of iterations and compare with the reference
echo
grep "it = ${NB_IT}" ${OUTPUT_VALIDATION}| awk -v REFERENCE=$REFERENCE -v ERROR=$ERROR '{ print "Energy=",$12,", Reference=",REFERENCE,", Error=", ERROR,", Reference-energy=",REFERENCE-$12}' 
RESU=`grep "it = ${NB_IT}" ${OUTPUT_VALIDATION}| awk -v REFERENCE=$REFERENCE -v ERROR=$ERROR '{ if( abs(REFERENCE-$12)< ERROR) print "OK" } \
function abs(a) { \
   if(a<0) a=-a; \
   return a; \
}' `
if [ "$RESU" != "OK" ] 
then
echo Energy is false
exit 1
else
echo Energy is right 
fi


