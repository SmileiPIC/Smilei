from tools import *
from math import ceil
from os import path
s = os.sep

def run_irene_a64fx(command, dir, mode, options, parameters):
    """
    Run the job on the Irene A64FX partition.

    Inputs:
    - command: command to run
    - dir: working directory
    - mode: compilation or execution
    """

    env = "module purge\n" \
            +"module load fujitsu mpi\n" \
            +"module load python3/3.8.10\n" \
            +"export HDF5_ROOT=/ccc/work/cont003/mds/lobetmat/Libraries/hdf5-1.12_fujitsu_a64fx/install/\n" \
            +"export PATH=${HDF5_ROOT}/bin:${PATH}\n" \
            +"export LD_LIBRARY_PATH=${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}\n" \
            +"export HDF5_ROOT_DIR=${HDF5_ROOT}\n" \
            +"export LD_LIBRARY_PATH=/ccc/products/ucx-1.10.1/system/default/lib:$LD_LIBRARY_PATH\n" \
            +"export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/ccc/products2/python3-3.8.10/Rhel_8__aarch64-a64fx/system/default/install_tree/python/3.8.10/s2azw3pgbfzhfcf44tvnh652pju2vtyj/lib:/ccc/products2/python3-3.8.10/Rhel_8__aarch64-a64fx/system/default/install_tree/gettext/0.21/i4aacmkl6pqxumiqfa36455yfhoojidl/lib\n"

    if (mode == "compilation"):

        # Create script
        with open(parameters['exec_script'], 'w') as exec_script_desc:
            #NODES=((int(options['mpi'])*int(options['omp'])-1)/24)+1
            NODES=int(ceil(options['mpi']/4.))
            options['ppn'] = 48
            exec_script_desc.write(
                "#!/bin/bash\n"
                +"#MSUB --job-name=smilei\n"
                +"#MSUB -N 1\n"
                +"#MSUB -n 48\n"
                +"#MSUB --output=output\n"
                +"#MSUB --error=error\n"
                +"#MSUB -q a64fx\n"
                +"#MSUB --time="+max_time+"\n"
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
            NODES=int(ceil(options['mpi']/4.))
            options['ppn'] = 48
            exec_script_desc.write(
                "#!/bin/bash\n"
                +"#MSUB --job-name=smilei\n"
                +"#MSUB -N "+str(NODES)+"\n"
                +"#MSUB -n "+str(options['mpi'])+"\n"
                +"#MSUB -c "+str(options['omp'])+"\n"
                +"#MSUB --output=output\n"
                +"#MSUB --error=error\n"
                +"#MSUB -q a64fx\n"
                +"#MSUB --time="+options['max_time']+"\n"
                +"#MSUB -A {}\n".format(options['accout'])
                +"#MSUB -m work,scratch\n"
                + env
                +"export OMP_NUM_THREADS="+str(options['omp'])+" \n"
                +"export OMP_SCHEDULE=DYNAMIC \n"
                +"export OMP_PROC_BIND=true \n"
                +"set -x\n"
                +"cd "+dir+" \n"
                +"cd ${BRIDGE_MSUB_PWD} \n"
                +"module list 2> module.log\n"
                +command+" \n"
                +"echo $? > exit_status_file \n"
            )

    # Run command
    JOB = "ccc_msub  "+parameters['exec_script']
    launch_job(command, JOB, dir, options['max_time_seconds'], parameters['output_file'], repeat=2)
