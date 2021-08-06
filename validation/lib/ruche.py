from launch_job import *
from math import ceil

def run_ruche(command, dir, mode, options):
    """
    Run the command `command` on the RUCHE system.

    Inputs:
    - command: command to run
    - dir: working directory
    """
    EXIT_STATUS="100"
    exit_status_fd = open(dir+s+"exit_status_file", "w")
    exit_status_fd.write(str(EXIT_STATUS))
    exit_status_fd.close()
    # Create script
    with open(parameters['exec_script'], 'w') as exec_script_desc:
        #NODES=((int(options['mpi'])*int(options['omp'])-1)/24)+1
        NODES=int(ceil(options['mpi']/2.))
        options['ppn'] = 20
        exec_script_desc.write(
            "#!/bin/bash\n"
            +"#SBATCH --job-name=smilei\n"
            +"#SBATCH --nodes="+str(NODES)+"\n"
            +"#SBATCH --ntasks="+str(options['mpi'])+"\n"
            +"#SBATCH --cpus-per-task="+str(options['omp'])+"\n"
            +"#SBATCH --output=output\n"
            +"#SBATCH --error=error\n"
            +"#SBATCH --time="+options['max_time']+"\n"
            +"#SBATCH --partition=cpu_short\n"
            +"module purge\n"
            +"module load zlib/1.2.9/gcc-9.2.0\n"
            +"module load anaconda3/2020.02/gcc-9.2.0\n"
            +"module load intel/19.0.3/gcc-4.8.5\n"
            +"module load intel-mpi/2019.3.199/intel-19.0.3.199\n"
            +"module load hdf5/1.10.6/intel-19.0.3.199-intel-mpi\n"
            +"export OMP_NUM_THREADS="+str(options['omp'])+" \n"
            +"export OMP_SCHEDULE=DYNAMIC \n"
            +"export OMP_PLACES=cores \n"
            +"export KMP_AFFINITY=verbose \n"
            +"export I_MPI_CXX=icpc \n"
            +"export HDF5_ROOT_DIR=/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/intel-19.0.3.199/hdf5-1.10.6-na3ilncuwbx2pdim2xaqwf23sgqza6qo \n"
            +"ulimit -s unlimited \n"
            +"cd ${SLURM_SUBMIT_DIR} \n"
            +"#Specify the number of sockets per node in -mca orte_num_sockets \n"
            +"#Specify the number of cores per sockets in -mca orte_num_cores \n"
            +"cd "+dir+" \n"
            +"module list 2> module.log\n"
            +command+" \n"
            +"echo $? > exit_status_file \n"
        )
    # Run command
    COMMAND = "sbatch  "+parameters['exec_script']
    try:
        check_call(COMMAND, shell=True)
    except CalledProcessError:
        # if command qsub fails, exit with exit status 2
        #Retry once in case the server was rebooting
        if VERBOSE :
            print(  "sbatch command failed once: `"+COMMAND+"`")
            print(  "Wait and retry")
        sleep(10)
        try:
            check_call(COMMAND, shell=True)
        except CalledProcessError:
            if dir==WORKDIR:
                os.chdir(WORKDIR_BASE)
                rmtree(WORKDIR)
            if VERBOSE :
                print(  "sbatch command failed twice: `"+COMMAND+"`")
                print(  "Exit")
            sys.exit(2)
    if VERBOSE:
        print( "Submitted job with command `"+command+"`")
        print( " -> max duration: {} s".format(max_time_seconds))
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
