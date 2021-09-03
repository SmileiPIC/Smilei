from tools import *
from math import ceil
from os import path
s = os.sep

def run_poincare(command, dir, mode, options, parameters):
    # Create script
    with open(parameters['exec_script'], 'w') as exec_script_desc:
        exec_script_desc.write(
            "# environnement \n"
            +"module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu 2>&1 > /dev/null\n"
            +"unset LD_PRELOAD\n"
            +"export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0\n"
            +"# \n"
            +"# execution \n"
            +"export OMP_NUM_THREADS="+str(options['omp'])+"\n"
            +command+" \n"
            +"echo $? > exit_status_file \n"
            +"exit $? "
        )
    JOB = "/bin/bash "+parameters['exec_script']+" > "+EXEC_SCRIPT_OUT+" 2>&1"
    launch_job(command, JOB, dir, options['max_time_seconds'], parameters['output_file'], repeat=2)
