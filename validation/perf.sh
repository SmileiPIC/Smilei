# OpenMP performance analysis of smilei depending of cluster width
#
# Marie Fle IDRIS 
#
# Usage : ./perf.sh -i input -c clrw
# input : input file, must be located in tests directory, default: tst3_thermal_expansion.in
# clrw : cluster width (default 1)
#
# * Compute n_space with value of res_space and grid_length in the input file
# * Read CLRW (cluster width) from stdout while CLRW does not divide n_space
# * Put CLRW into the input file
# * Submit the job perf.ll:
#   - executes smilei sequentially
#   - output the cpu_time of the main loop,  particules, maxwell, diagnostics, densities for 1,2, 4, 8, 16 OpenMP threads 
#   - compute the efficiency of the main loop,  particules, maxwell, diagnostics, densities for 2, 4, 8, 16 OpenMP threads 
#
# The results are printed in the file: PERFS_<input_file>_CLRW_<clrw>
#
# Note : the value of the OpenMP schedule can be changed in the job perf.ll
#
# environnement
ROOT="/gpfshome/mds/staff/mfle/SMILEI/smilei"
while getopts "i:c:" OPTION
do
        case $OPTION in
        i )    INPUT=$OPTARG;;
        c )    CLRW=${OPTARG};;
        esac
done
INPUT=${INPUT:-tst3_thermal_expansion.in}
CLRW=${CLRW:-1}
echo $INPUT
INPUT_PATH=${ROOT}/tests/${INPUT}
INPUT_CLRW=${INPUT}_CLRW
INPUT_PATH_CLRW=${ROOT}/tests/${INPUT_CLRW}
JOB_perf=${ROOT}/scripts/perf.ll
JOB_perf_temp=${ROOT}/scripts/perf_temp.ll
#
#
# Compute n_space with value of res_space
set -x
N_SPACE=`awk -v n_space=1 -F "=" '$1 == "res_space" {n_space=n_space*$2} $1 == "grid_length" {n_space=n_space*$2} END {print n_space}'   ${INPUT_PATH} `
#
# Read CLRW from stdout while CLRW does not divide n_space
cp ${INPUT_PATH} ${INPUT_PATH_CLRW}
if [ "$CLRW" != "1" ]
  then
  CONTINUE=`echo $N_SPACE $CLRW |  awk ' $1%$2 == "0" {print "1" }'`
  while [ "$CONTINUE" != "1" ]
    do
    echo Entrer un diviseur de $N_SPACE
    read CLRW 
    CONTINUE=`echo $N_SPACE $CLRW |  awk ' $1%$2 == "0" {print "1" }'`
  done
fi
#
# Put CLRW into the input file  
sed " 25a cluster_width=${CLRW} " ${INPUT_PATH} > ${INPUT_PATH_CLRW} 
#
# Build the job perf.ll                        
cat > ${JOB_perf} << 'eof'
#@ job_name = smilei
#@ output = $(job_name).$(jobid).out
#@ error  = $(job_name).$(jobid).err
#@ job_type = mpich
#@ wall_clock_limit = 02:30:00
#@ class = clallmds
#@ total_tasks      = 1
#@ environment = COPY_ALL
#@ queue
#
#
 function exe {
     mpirun  -hostfile ${LOADL_HOSTFILE} -np ${LOADL_TOTAL_TASKS} ${ROOT}/src/smilei ${INPUT_PATH} > ${OUTPUT_PERF}
 }
 function cpu_time {
     TIME_LOOP=`grep "Time in time loop" ${OUTPUT_PERF}| awk -v TIME_LOOP_SEQ=$TIME_LOOP_SEQ -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{print $6 }' `
     PARTICULES=`grep particles ${OUTPUT_PERF}| grep  "%" | awk -v PARTICULES_SEQ=$PARTICULES_SEQ -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{print $2 }'`
     MAXWELL=`grep maxwell ${OUTPUT_PERF} | awk  -v MAXWELL_SEQ=$MAXWELL_SEQ -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{print $2 }' `
     DIAGNOSTICS=`grep diagnostics ${OUTPUT_PERF} | awk -v DIAGNOSTICS_SEQ=$DIAGNOSTICS_SEQ -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{print $2 }' `
     DENSITIES=`grep densities ${OUTPUT_PERF} | awk -v DENSITIES_SEQ=$DENSITIES_SEQ  -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{print $2 }' `
#     echo ${CLRW}"          "${LOADL_TOTAL_TASKS}"          "${OMP_NUM_THREADS}"         "${OMP_SCHEDULE}"      "${TIME_LOOP}"      " ${PARTICULES}"      " ${MAXWELL}"     " ${DIAGNOSTICS}"     " ${DENSITIES}   >> ${PERFS}
     printf " %-11i %-11i %-11i %-11s %-11.2f %-11.2f %-11.2f %-11.2f %-11.2f\n" ${CLRW} ${LOADL_TOTAL_TASKS} ${OMP_NUM_THREADS} ${OMP_SCHEDULE} ${TIME_LOOP} ${PARTICULES}  ${MAXWELL} ${DIAGNOSTICS} ${DENSITIES}   >> ${PERFS}
 }
 function efficiency {
     echo "EFFICIENCY :"
     TIME_LOOP=`grep "Time in time loop" ${OUTPUT_PERF}| awk -v TIME_LOOP_SEQ=$TIME_LOOP_SEQ -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{print TIME_LOOP_SEQ/($6*NBMPI*NBOMP) }' `
     PARTICULES=`grep particles ${OUTPUT_PERF}| grep  "%" | awk -v PARTICULES_SEQ=$PARTICULES_SEQ -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{printf PARTICULES_SEQ/($2*NBMPI*NBOMP) }'`
     MAXWELL=`grep maxwell ${OUTPUT_PERF}|  awk  -v MAXWELL_SEQ=$MAXWELL_SEQ -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{print MAXWELL_SEQ/($2*NBMPI*NBOMP) }' `
     DIAGNOSTICS=`grep diagnostics ${OUTPUT_PERF} | awk -v DIAGNOSTICS_SEQ=$DIAGNOSTICS_SEQ -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{print DIAGNOSTICS_SEQ/($2*NBMPI*NBOMP) }' `
     DENSITIES=`grep densities ${OUTPUT_PERF} | awk -v DENSITIES_SEQ=$DENSITIES_SEQ  -v NBMPI=${LOADL_TOTAL_TASKS} -v NBOMP=$OMP_NUM_THREADS '{printf DENSITIES_SEQ/($2*NBMPI*NBOMP) }' `
     printf " %-11i %-11i %-11i %-11s %-11.2f %-11.2f %-11.2f %-11.2f %-11.2f\n" ${CLRW} ${LOADL_TOTAL_TASKS} ${OMP_NUM_THREADS} ${OMP_SCHEDULE} ${TIME_LOOP} ${PARTICULES}  ${MAXWELL} ${DIAGNOSTICS} ${DENSITIES}   >> ${PERFS}
 }
#
# Environment
set -x
module load intel openmpi ddt hdf5/1.8.10_intel_openmpi
ROOT="/gpfshome/mds/staff/mfle/SMILEI/smilei"
INPUT=input                           
INPUT_PATH=${ROOT}/tests/${INPUT}
OUTPUT_SEQ=/tmp/output_seq                         
OUTPUT_PERF=/tmp/output_efficiency                         
CLRW=`grep ^cluster_width ${INPUT_PATH} | awk -F "=" '{print $2}'`
PERFS="./PERFS_"${INPUT}_${CLRW}
echo ${PERFS}
# 
# Titles
ssh poincareint01 'cd /gpfshome/mds/staff/mfle/SMILEI/smilei/src ;git status | head -1 |grep -i branch' > /tmp/temp2
awk '{print $4}' /tmp/temp2 > /tmp/temp3
read BRANCH < /tmp/temp3
ssh poincareint01 'cd /gpfshome/mds/staff/mfle/SMILEI/smilei/src ;git log | head -1 | grep -i commit' > /tmp/temp2
read COMMIT < /tmp/temp2
echo "Branch : "${BRANCH}", "${COMMIT} >> ${PERFS}
echo "Namelist : " $INPUT >> ${PERFS}
echo "TEMPS CPU :" >> ${PERFS}
echo " CLRW        NB_PROC     NB_THREADS  SCHEDULE    Time_loop   Particules  Maxwell     Diagnostics Densities " >> ${PERFS}               
# 
# Sequential execution
export OMP_NUM_THREADS=1
export OMP_SCHEDULE="static"
mpirun  -hostfile ${LOADL_HOSTFILE} -np 1  ${ROOT}/src/smilei ${INPUT_PATH} > ${OUTPUT_SEQ}
#
# Register sequential perfs
TIME_LOOP_SEQ=`grep "Time in time loop" ${OUTPUT_SEQ} | awk '{print $6}' `
PARTICULES_SEQ=`grep particles ${OUTPUT_SEQ}| grep "%"| awk '{print $2}' `
MAXWELL_SEQ=`grep maxwell ${OUTPUT_SEQ} | awk '{print $2}' `
DIAGNOSTICS_SEQ=`grep diagnostics ${OUTPUT_SEQ} | awk '{print $2}' `
DENSITIES_SEQ=`grep densities ${OUTPUT_SEQ} | awk '{print $2}' `
#
# output sequentiel cpu time               
printf  " %-11i %-11i %-11i %-11s %-11.2f %-11.2f %-11.2f %-11.2f %-11.2f\n"   ${CLRW} ${LOADL_TOTAL_TASKS} ${OMP_NUM_THREADS} ${OMP_SCHEDULE} ${TIME_LOOP_SEQ} ${PARTICULES_SEQ}  ${MAXWELL_SEQ} ${DIAGNOSTICS_SEQ} ${DENSITIES_SEQ}   >> ${PERFS}
#
# Read CLRW
# cpu time with 2 OpenMP         
export OMP_NUM_THREADS=2
export OMP_SCHEDULE="static"
exe
cpu_time
#
# cpu time with 4 OpenMP        
export OMP_NUM_THREADS=4
export OMP_SCHEDULE="static,15"
exe
cpu_time
#
# cpu time with 8 OpenMP        
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="static,10"
exe
cpu_time
#
# cpu time with 16 OpenMP       
export OMP_NUM_THREADS=16
export OMP_SCHEDULE="static,5"
exe
cpu_time
#
echo "EFFICIENCY :" >> ${PERFS}
# Efficiency with 2 OpenMP         
export OMP_NUM_THREADS=2
export OMP_SCHEDULE="static"
exe
efficiency
#
# Efficiency with 4 OpenMP        
export OMP_NUM_THREADS=4
export OMP_SCHEDULE="static,15"
exe
efficiency
#
# Efficiency with 8 OpenMP        
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="static,10"
exe
efficiency
#
# Efficiency with 16 OpenMP       
export OMP_NUM_THREADS=16
export OMP_SCHEDULE="static,5"
exe
efficiency
#
rm ${OUTPUT_SEQ} ${OUTPUT_PERF}
eof
#
# Put the name of the input file into the job perf.ll
sed -e "1,\$s/input/${INPUT_CLRW}/g" ${JOB_perf} > ${JOB_perf_temp}
#
# Submit the job perf.ll
llsubmit -s ${JOB_perf_temp}
rm ${JOB_perf} ${JOB_perf_temp} ${INPUT_PATH_CLRW=} 
