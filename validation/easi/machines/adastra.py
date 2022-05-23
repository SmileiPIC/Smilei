from . import Machine


class MachineAdastra(Machine):
    """
    As of 17 of march 22, this class, while named Adastra, instead targets the 
    Adastra porting machines. The real Adastra environment should not change much.
    """

    # If you are editing this file, be carful with the python format brackets '{' '}'.
    # You may need to escape some.
    the_slurm_script = """#!/bin/bash
#SBATCH --job-name=smilei_validation
#SBATCH --nodes={the_node_count}               # Number of nodes
#SBATCH --ntasks={the_mpi_process_count}       # Number of MPI ranks
# #SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task={the_omp_thread_count} # Number of cores per MPI rank
#SBATCH --gpus-per-task={the_gpu_count}        # Number of gpu per MPI rank
#SBATCH --output=output
#SBATCH --error=output                         # stderr and stdout in the same file
#SBATCH --time={the_maximum_task_duration}

# Dump all executed commands (very, VERY verbose)
# set +x

# TODO(Etienne M): dunno the partition/feature/constraints for adastra yet
# # --feature=MI200
# #SBATCH --gpus=mi100:1 or mi200:1

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

# We should need only that to run the rest is loaded by default
module load rocm/4.5.0

# Info on the node
rocm-smi
rocminfo

# MPICH Gpu support
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_IPC_ENABLED=1
# export MPICH_OPTIMIZED_MEMCPY=2

# MPI to GPU binding
# #!/bin/bash
# # The following script could be useful to bind a mpiproc to a given gpu on a given node https://slurm.schedmd.com/sbatch.html
# # https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/user-guide.html#gpu-enumeration
# # One may check that it works: $ ROCR_VISIBLE_DEVICES=none rocminfo
# # rocm-smi ignores the env variable
#
# export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID
# exec $*

# Omp tuning
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_SCHEDULE=dynamic
# You may want to change "cores", to "threads". But hyperthreading, for an well 
# optimized apps is generally not something you want (ROB buffer should be full).
export OMP_PLACES=cores

# Omp target debug
# export CRAY_ACC_DEBUG=3

LaunchSRun() {{
    module list

    srun "$@" > {the_output_file} 2>&1
}}

# You must have built smilei with the 'perftools' module loaded!
LaunchSRunPatProfile() {{
    module load perftools-base/21.12.0
    module load perftools

    # # Enable extra verbose tracing, can be very useful, produces a lot of data
    # export PAT_RT_SUMMARY=0
    export PAT_RT_MPI_THREAD_REQUIRED=3

    # Assuming "$1" is an executable
    pat_build -g hip,io,mpi,cuda -w -f $1 -o instrumented_executable

    LaunchSRun instrumented_executable ${{@:2}}
}}

# Try to use this profiling on only one GPU
LaunchRocmProfile() {{
    # Kernel stats
    echo 'pmc : VALUUtilization VALUBusy SALUBusy MemUnitBusy MemUnitStalled L2CacheHit Wavefronts' > hw_counters.txt
    LaunchSRun bash -c "rocprof -i hw_counters.txt --stats -o hw_counters_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"

    # hip RT tracing, --trace-period 30s:30s:1m is bugged
    # LaunchSRun bash -c "rocprof --hip-trace --stats -o hip_trace_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"
}}

LaunchSRun {a_task_command} {a_task_command_arguments}
# LaunchSRunPatProfile {a_task_command} {a_task_command_arguments}
# LaunchRocmProfile {a_task_command} {a_task_command_arguments}

kRETVAL=$?

# Put the result in the slurm output file.
cat {the_output_file}

echo "The task ended at = $(date)"

echo -n $kRETVAL > exit_status_file
exit $kRETVAL
"""

    def __init__(self, smilei_path, options):
        print('##############################################')
        print('# Preparing validation for the Adastra machine')
        print('##############################################')

        from math import ceil

        self.smilei_path = smilei_path
        self.options = options

        # Use the config flag "gpu_amd" to compile for GPU or remove it for CPU only
        the_make_command = 'make machine="adastra" config="' + (self.options.compile_mode if self.options.compile_mode else '') + '" -j'
        self.COMPILE_COMMAND = the_make_command + ' > ' + self.smilei_path.COMPILE_OUT + ' 2> '+self.smilei_path.COMPILE_ERRORS
        self.CLEAN_COMMAND = 'make clean > /dev/null 2>&1'

        # This will describe and schedule the tasks
        self.JOB = 'sbatch ' + self.smilei_path.exec_script
        # This'll start the tasks
        self.RUN_COMMAND = self.smilei_path.workdirs + '/smilei'
        
        if self.options.nodes:
            self.NODES = self.options.nodes
        else:
            # 4 MI200~GPUs per node. Each GPU contains 2 GCD
            kGPUPerNode = 4 * 2
            self.NODES = int(ceil(self.options.mpi / kGPUPerNode))

    def run(self, arguments, dir):
        """
        Run a simulation
        """

        # Write the slurm script
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write(self.the_slurm_script.format(a_task_command=self.RUN_COMMAND, 
                                                 a_task_command_arguments=arguments, 
                                                 the_node_count=self.NODES, 
                                                 the_mpi_process_count=self.options.mpi, 
                                                 the_maximum_task_duration=self.options.max_time, 
                                                 the_omp_thread_count=self.options.omp,
                                                 the_output_file=self.smilei_path.output_file,
                                                 the_gpu_count='1' if 'gpu_amd' in self.options.compile_mode else '0'))

        # Schedule the task(s)
        self.launch_job(self.RUN_COMMAND + ' ' + arguments, 
                        self.JOB, 
                        dir, 
                        self.options.max_time_seconds, 
                        self.smilei_path.output_file, 
                        repeat=2)
