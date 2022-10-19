from . import Machine

class MachineAdastra(Machine):
    """
    As of 22/03/17, this class, while named Adastra, instead targets the 
    Adastra porting machines. The real Adastra environment should not change much.
    """

    # If you are editing this file, be carful with the python format brackets '{' '}'.
    # You may need to escape some.
    the_slurm_script = """#!/bin/bash
#SBATCH --job-name=smilei_validation
#SBATCH --nodes={the_node_count}               # Number of nodes
#SBATCH --ntasks={the_mpi_process_count}       # Number of MPI ranks
#SBATCH --cpus-per-task={the_omp_thread_count} # Number of cores per MPI rank
# #SBATCH --gpus-per-node=8
# #SBATCH --gres=gpu:8
#SBATCH --gpus-per-task={the_gpu_count}        # Number of gpu per MPI rank, Dont use that, it create problems with MPICH_GPU_SUPPORT_ENABLED=1 on intra node GPU to GPU coms
# #SBATCH --ntasks-per-gpu={the_gpu_count}       # Number of MPI rank per task (may be useful to oversubscribe, if you cant fill the whole gpu)
#SBATCH --gpu-bind=closest                     # For a given task and its associated numa, bind the closest GPU(s) (maybe more than one) to the numa. As of 2022/06, breaks script GPU visibility of the task, ROCR_VISIBLE_DEVICES must be used to counter the effect
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --hint=memory_bound                    # Maximize memory bandwidth by spreading the task on the numa and physical cores. Note: https://slurm.schedmd.com/mc_support.html
#SBATCH --distribution=block:cyclic:cyclic     # Spread linearly (stride 1) across nodes, round robin across numa of a node and cores of the numas of a node (stride of the number of core in a numa) : Node0, Node1, Node2.. then Core0 of Numa0, Core0 of Numa1 etc..
#SBATCH --output=output
#SBATCH --error=output                         # stderr and stdout in the same file
#SBATCH --time={the_maximum_task_duration}

# Dump all executed commands (very, VERY verbose)
# set +x

# set -e

# TODO(Etienne M): dunno the partition/feature/constraints for adastra yet
# # --feature=MI200
# #SBATCH --gpus=mi100:1 or mi200:1

echo "Date              = $(date -R) | $(date +%s)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

# Build  the environment (delegating this to a script would be better)
module purge
module load craype-network-ofi craype-x86-rome libfabric/1.13.1
module load PrgEnv-cray/8.3.3 cce/13.0.1
module load cray-mpich/8.1.13
module load rocm/4.5.0
module load craype-accel-amd-gfx908 # MI100
# module load craype-accel-amd-gfx90a # MI250X
module load cray-hdf5-parallel/1.12.0.6 cray-python/3.9.7.1

# Info on the node
rocm-smi
rocminfo

# Omp tuning
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_SCHEDULE=dynamic
export OMP_PLACES=cores

# These variables may trigger a huge flood of information when using multiple 
# threads per MPI + GPU support
export OMP_DISPLAY_AFFINITY=TRUE # Unused by the CCE omp runtime
export CRAY_OMP_CHECK_AFFINITY=TRUE

# MPICH info
export MPICH_VERSION_DISPLAY=1
export MPICH_ENV_DISPLAY=1
export MPICH_CPUMASK_DISPLAY=1
export MPICH_MPIIO_HINTS_DISPLAY=1

# MPICH general
export MPICH_ABORT_ON_ERROR=1 # Errors are not checked by Smilei, they must not happen

# MPICH GPU support (support for functionality allowing mpi to understand GPU buffers)
export MPICH_GPU_SUPPORT_ENABLED=1

# Omp target debug
# export CRAY_ACC_DEBUG=3
# export CRAY_ACC_TIME=1 # Not viable, the reported times are wrong ?

# Amd runtime debug
# export AMD_LOG_LEVEL=4

LaunchSRun() {{
    module list

    # Task/GPU binding. This mapping can be automated using Slurm's 
    # --gpus-per-task but this lead a GPU memory accessibility/visibility  
    # during MPI coms (process_vm_readv error). The workaround is:
    # 1): To allocate GPUs using --gres/--gpus-per-node and not define any 
    # mapping. Or to allocate GPUs using --gpus-per-task and break the mapping
    # using --gpu-bind=closest.
    # 2): Use ROCR_VISIBLE_DEVICES to do the mapping, associating a GPU to a 
    # given, user specified task/rank. This may make things harder when we wanna
    # map a task to a given GPU based on hardware/NUMA characteristics.
    # 
    # The following script could be useful to bind an mpiproc to a given gpu on 
    # a given node https://slurm.schedmd.com/sbatch.html
    # https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/user-guide.html#gpu-enumeration
    # One may check that it works: $ ROCR_VISIBLE_DEVICES=none rocminfo
    # rocm-smi ignores the env variable
    #
    # Note: When not using ROCm, this script should not change the behavior of 
    # any software component.
    #
    cat <<EOF > soft_gpu_visibility_restrict.sh
#!/bin/bash
# Note: When not using ROCm, this script should not change the behavior of any 
# software component.

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec "\$@"
EOF

    chmod +x soft_gpu_visibility_restrict.sh

    srun --cpu-bind=verbose soft_gpu_visibility_restrict.sh "$@" > {the_output_file} 2>&1
    # srun strace "$@" > {the_output_file} 2>&1
    # kCmd="if [ \${{SLURM_PROCID}} -eq 0 ]; then strace $@; else $@; fi"
    # srun bash -c "$kCmd" > {the_output_file} 2>&1
}}

# You must have built smilei with the 'perftools' module loaded!
LaunchSRunPatProfile() {{
    module load perftools-base/21.12.0
    module load perftools

    # # Enable extra verbose tracing, can be very useful, produces a lot of data
    # export PAT_RT_SUMMARY=0
    export PAT_RT_MPI_THREAD_REQUIRED=3

    # Assuming "$1" is an executable

    # GPU profiling
    pat_build -g hip,io,mpi -w -f $1 -o instrumented_executable

    # CPU profiling
    # export PAT_RT_SAMPLING_INTERVAL=1000; # 1ms interval
    # export PAT_RT_SAMPLING_MODE=3
    # export PAT_RT_EXPERIMENT=samp_cs_time
    # # export PAT_RT_PERFCTR=default # PAPI_get_component_info segfault in rsmi_func_iter_value_get
    # # -u make the program crash (abi break) or silently corrupts it state
    # pat_build -g mpi,syscall,io,omp,hdf5 -w -f $1 -o instrumented_executable

    LaunchSRun ./instrumented_executable ${{@:2}}
}}

# Try to use this profiling on only one GPU
LaunchRocmProfile() {{
    # Basic kernel dump ("tid","grd","wgr","lds","scr","vgpr","sgpr","fbar","sig","obj","DispatchNs","BeginNs","EndNs","CompleteNs","DurationNs") + consolidated kernel stats
    # Low overhead
    # LaunchSRun bash -c "rocprof --stats -o stats_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"

    #  VALUUtilization : The percentage of active vector ALU threads in a wave. A lower number can mean either more thread divergence in a wave or that the work-group size is not a multiple of 64. Value range: 0% (bad), 100% (ideal - no thread divergence).
    #         VALUBusy : The percentage of GPUTime vector ALU instructions are processed. Value range: 0% (bad) to 100% (optimal).
    #         SALUBusy : The percentage of GPUTime scalar ALU instructions are processed. Value range: 0% (bad) to 100% (optimal).
    #        FetchSize : The total kilobytes fetched from the video memory. This is measured with all extra fetches and any cache or memory effects taken into account.
    #        WriteSize : The total kilobytes written to the video memory. This is measured with all extra fetches and any cache or memory effects taken into account.
    #       L2CacheHit : The percentage of fetch, write, atomic, and other instructions that hit the data in L2 cache. Value range: 0% (no hit) to 100% (optimal).
    #      MemUnitBusy : The percentage of GPUTime the memory unit is active. The result includes the stall time (MemUnitStalled). This is measured with all extra fetches and writes and any cache or memory effects taken into account. Value range: 0% to 100% (fetch-bound).
    #   MemUnitStalled : The percentage of GPUTime the memory unit is stalled. Try reducing the number or size of fetches and writes if possible. Value range: 0% (optimal) to 100% (bad).
    # WriteUnitStalled : The percentage of GPUTime the Write unit is stalled. Value range: 0% to 100% (bad).
    #  ALUStalledByLDS : The percentage of GPUTime ALU units are stalled by the LDS input queue being full or the output queue being not ready. If there are LDS bank conflicts, reduce them. Otherwise, try reducing the number of LDS accesses if possible. Value range: 0% (optimal) to 100% (bad).
    #  LDSBankConflict : The percentage of GPUTime LDS is stalled by bank conflicts. Value range: 0% (optimal) to 100% (bad).

    # Basic kernel dump + consolidated kernel stats + hw counters (kernel duration/stats may be completely broken)
    # High overhead (~15%)
    # echo 'pmc : VALUUtilization VALUBusy SALUBusy GPUBusy MemUnitBusy MemUnitStalled Wavefronts FetchSize WriteSize' > hw_counters.txt
    #
    # echo 'pmc : VALUUtilization VALUBusy SALUBusy L2CacheHit MemUnitBusy MemUnitStalled WriteUnitStalled ALUStalledByLDS LDSBankConflict' > hw_counters.txt
    #
    # How to do an AMD roofline:
    # https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#roofline-profiling
    # echo 'pmc : SQ_INSTS_VALU_ADD_F16 SQ_INSTS_VALU_MUL_F16 SQ_INSTS_VALU_FMA_F16 SQ_INSTS_VALU_TRANS_F16' > hw_counters.txt
    # echo 'pmc : SQ_INSTS_VALU_ADD_F32 SQ_INSTS_VALU_MUL_F32 SQ_INSTS_VALU_FMA_F32 SQ_INSTS_VALU_TRANS_F32' > hw_counters.txt
    # echo 'pmc : SQ_INSTS_VALU_ADD_F64 SQ_INSTS_VALU_MUL_F64 SQ_INSTS_VALU_FMA_F64 SQ_INSTS_VALU_TRANS_F64' > hw_counters.txt
    # echo 'pmc : SQ_INSTS_VALU_MFMA_MOPS_F16 SQ_INSTS_VALU_MFMA_MOPS_BF16 SQ_INSTS_VALU_MFMA_MOPS_F32 SQ_INSTS_VALU_MFMA_MOPS_F64' > hw_counters.txt
    # echo 'pmc : TCC_EA_RDREQ_32B_sum TCC_EA_RDREQ_sum TCC_EA_WRREQ_sum TCC_EA_WRREQ_64B_sum' > hw_counters.txt
    #
    # LaunchSRun bash -c "rocprof -i hw_counters.txt --stats -o hw_counters_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"

    # HIP RT tracing + consolidated kernel stats
    # LaunchSRun bash -c "rocprof --hip-trace --stats -o hip_trace_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"
    # ROCm RT (low level) tracing + consolidated kernel stats
    # LaunchSRun bash -c "rocprof --hsa-trace --stats -o hsa_trace_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"
    # HIP/HSA tracing + consolidated kernel stats | After the program exits, rocprof takes a long time to process the generated data !
    LaunchSRun bash -c "rocprof --sys-trace --stats -o sys_trace_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"
}}

LaunchSRun {a_task_command} {a_task_command_arguments}
# LaunchSRunPatProfile {a_task_command} {a_task_command_arguments}
# LaunchRocmProfile {a_task_command} {a_task_command_arguments}

kRETVAL=$?

# Put the result in the slurm output file.
cat {the_output_file}

echo "The task ended at = $(date -R) | $(date +%s)"

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
        the_make_command = 'make machine="adastra" config="' + (self.options.compile_mode if self.options.compile_mode else '') + (' verbose' if self.options.verbose else '') + '" -k -j'
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
            kMaxGPUPerNode = 4 * 2
            # self.gpu_per_node = self.options.mpi if self.options.mpi <= kMaxGPUPerNode else kMaxGPUPerNode
            self.NODES = int(ceil(self.options.mpi / kMaxGPUPerNode))

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
