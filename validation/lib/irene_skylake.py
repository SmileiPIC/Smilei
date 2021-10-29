from .tools import *
from math import ceil
import os
s = os.sep

def run_irene_skylake(command, dir, mode, options, parameters):
    """
    Run the job on the Irene Skylake system.

    Inputs:
    - command: command to run
    - dir: working directory
    - mode: compilation or execution
    """

    env = "module purge\n" \
            +"module load intel/20.0.4\n" \
            +"module load mpi/intelmpi/20.0.4\n" \
            +"module load flavor/hdf5/parallel hdf5/1.8.20\n" \
            +"module load python3\n" \
            +"export PATH=${HDF5_ROOT}/bin:${PATH}\n" \
            +"export LD_LIBRARY_PATH=${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}\n" \
            +"export HDF5_ROOT_DIR=${HDF5_ROOT}\n"

    if (mode == "compilation"):

        # Create script
        with open(parameters['exec_script'], 'w') as exec_script_desc:
            exec_script_desc.write(
                "#!/bin/bash\n"
                +"#MSUB --job-name=smilei\n"
                +"#MSUB -N 1\n"
                +"#MSUB -n 48\n"
                +"#MSUB --output=output\n"
                +"#MSUB --error=error\n"
                +"#MSUB -q skylake\n"
                +"#MSUB --time="+options['max_time']+"\n"
                +"#MSUB -A {}\n".format(options['accout'])
                +"#MSUB -m work,scratch\n"
                + env
                +"set -x\n"
                +"cd "+dir+" \n"
                +"cd ${BRIDGE_MSUB_PWD} \n"
                +"module list 2> module.log\n"
                +command+" \n"
                +"echo $? > exit_status_file \n"
            )

    else:

        # Create script
        with open(parameters['exec_script'], 'w') as exec_script_desc:
            #NODES=((int(options['mpi'])*int(options['omp'])-1)/24)+1
            NODES=int(ceil(options['mpi']/2.))
            options['ppn'] = 24
            exec_script_desc.write(
                "#!/bin/bash\n"
                +"#MSUB --job-name=smilei\n"
                +"#MSUB -N "+str(NODES)+"\n"
                +"#MSUB -n "+str(options['mpi'])+"\n"
                +"#MSUB -c "+str(options['omp'])+"\n"
                +"#MSUB --output=output\n"
                +"#MSUB --error=error\n"
                +"#MSUB -q skylake\n"
                +"#MSUB --time="+options['max_time']+"\n"
                +"#MSUB -A {}\n".format(options['accout'])
                +"#MSUB -m work,scratch\n"
                +env
                +"export OMP_NUM_THREADS="+str(options['omp'])+" \n"
                +"export OMP_SCHEDULE=DYNAMIC \n"
                +"export OMP_PROC_BIND=true \n"
                +"export KMP_AFFINITY=verbose \n"
                +"set -x\n"
                +"ulimit -s unlimited \n"
                +"cd ${BRIDGE_MSUB_PWD} \n"
                +"cd "+dir+" \n"
                +"module list 2> module.log\n"
                +command+" \n"
                +"echo $? > exit_status_file \n"
            )
    # Run command
    JOB = "ccc_msub  "+parameters['exec_script']
    launch_job(command, JOB, dir, options['max_time_seconds'], parameters['output_file'], repeat=2, verbose=options['verbose'])
