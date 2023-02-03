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
#SBATCH --account={the_account}
#SBATCH --constraint={the_partition}&turbo
#SBATCH --nodes={the_node_count} --exclusive
#SBATCH --ntasks={the_mpi_process_count}
#SBATCH --cpus-per-task={the_reserved_thread_count} --gpus-per-node={the_gpu_per_node_count}
#SBATCH --threads-per-core=1 # --hint=nomultithread
#SBATCH --output=output --error=output
#SBATCH --time={the_maximum_task_duration}

################################################################################
# Time of day
################################################################################

echo "Date              = $(date -R) | $(date +%s)"
echo "Hostname          = $(hostname -s)"
echo "Working directory = $(pwd)"
echo ""
echo "Number of nodes allocated       = $SLURM_JOB_NUM_NODES"
echo "Number of tasks allocated       = $SLURM_NTASKS"
echo "Number of cores/tasks allocated = $SLURM_CPUS_PER_TASK"

################################################################################
# Environment setup
################################################################################

set -eu

# Build the environment (delegating this to a script would be better)
module purge

module load craype-accel-amd-gfx90a craype-x86-trento
module load PrgEnv-cray/8.3.3
module load cray-mpich/8.1.21 cray-hdf5-parallel/1.12.2.1 cray-python/3.9.13.1
module load amd-mixed/5.2.3

# OpenMP
# export OMP_DISPLAY_ENV=VERBOSE
export OMP_DISPLAY_AFFINITY=TRUE

export OMP_NUM_THREADS={the_omp_thread_count}
export OMP_SCHEDULE=DYNAMIC
export OMP_PLACES=THREADS
export OMP_PROC_BIND=FALSE

# export CRAY_ACC_DEBUG=3
# export CRAY_ACC_TIME=1 # Not viable. Are the reported times wrong ?

# AMD runtime
# export AMD_LOG_LEVEL=4

# MPICH info
export MPICH_VERSION_DISPLAY=1
export MPICH_ENV_DISPLAY=1
export MPICH_CPUMASK_DISPLAY=1
export MPICH_MPIIO_HINTS_DISPLAY=1
# export MPICH_OFI_VERBOSE=1
# export MPICH_OFI_NIC_VERBOSE=2

export MPICH_GPU_SUPPORT_ENABLED=1

export MPICH_ABORT_ON_ERROR=1 # Errors are not checked by Smilei, they must not happen

# Workaround to the deadlock we get in MPI communication
export FI_MR_CACHE_MONITOR=memhooks
# Ou:
# export FI_MR_CACHE_MAX_COUNT=0 

set +x

################################################################################
# Node information
################################################################################

ulimit -c unlimited
# ulimit -a

# NOTE: rocm-smi often fails.. and set -e ends the script. 'true' prevents that.
rocm-smi || true
# rocminfo
lscpu
numactl -H

################################################################################
# Helpers
################################################################################

LaunchSRun() {{
    module list

    SRUN_FLAGS_BASE="--mpi=cray_shasta --kill-on-bad-exit=1"

    SRUN_FLAGS_BINDING="--cpu-bind=none --mem-bind=none"
    # SRUN_FLAGS_BINDING="--cpus-per-task={the_reserved_thread_count} --gpu-bind=closest"

    SRUN_FLAGS_EXTRA="--label"

    cat <<EOF > soft_gpu_visibility_restrict.sh && chmod +x soft_gpu_visibility_restrict.sh
#!/bin/bash
export HIP_VISIBLE_DEVICES=\$SLURM_LOCALID
exec "\$@"
EOF

    SRUN_WRAPPER_SCRIPT="soft_gpu_visibility_restrict.sh"
    # SRUN_WRAPPER_SCRIPT="embind"
    # SRUN_WRAPPER_SCRIPT="strace"

    set -x
    srun ${{SRUN_FLAGS_BASE}} ${{SRUN_FLAGS_BINDING}} ${{SRUN_FLAGS_EXTRA}} -- ${{SRUN_WRAPPER_SCRIPT}} "$@" > {the_output_file} 2>&1
    set +x
}}

# You must have built smilei with the 'perftools' module loaded!
#
LaunchSRunPATProfile() {{
    module load perftools-base/22.09.0
    # module load perftools

    # # Enable extra verbose tracing, can be very useful, produces a lot of data
    # export PAT_RT_SUMMARY=0
    export PAT_RT_MPI_THREAD_REQUIRED=3

    # Assuming "$1" is an executable

    # GPU profiling
    # pat_build -g hip,io,mpi -w -f $1 -o instrumented_executable

    # CPU profiling
    # export PAT_RT_SAMPLING_INTERVAL=1000; # 1ms interval
    # export PAT_RT_SAMPLING_MODE=3
    # export PAT_RT_EXPERIMENT=samp_cs_time
    # # export PAT_RT_PERFCTR=default # PAPI_get_component_info segfault in rsmi_func_iter_value_get
    # # -u make the program crash (abi break) or silently corrupts it state
    # pat_build -g mpi,syscall,io,omp,hdf5 -w -f $1 -o instrumented_executable

    # LaunchSRun ./instrumented_executable ${{@:2}}

    LaunchSRun "$@"
}}

LaunchROCmProfile() {{
    # Basic kernel dump ("tid","grd","wgr","lds","scr","vgpr","sgpr","fbar","sig","obj","DispatchNs","BeginNs","EndNs","CompleteNs","DurationNs") + consolidated kernel stats
    # Low overhead
    LaunchSRun bash -c "rocprof --stats -o stats_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"

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
    # echo 'pmc : VALUUtilization VALUBusy L2CacheHit LDSBankConflict ALUStalledByLDS SALUBusy MemUnitBusy MemUnitStalled WriteUnitStalled' > hw_counters.txt
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
    # LaunchSRun bash -c "rocprof --sys-trace --stats -o sys_trace_\${{SLURM_JOBID}}-\${{SLURM_PROCID}}.csv $1 ${{@:2}}"
}}

LaunchSRun {a_task_command} {a_task_command_arguments}
# LaunchSRunPATProfile {a_task_command} {a_task_command_arguments}
# LaunchROCmProfile {a_task_command} {a_task_command_arguments}

kRETVAL=$?

# Put the result in the SLURM output file.
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
        
        if 'gpu_amd' in self.options.compile_mode:
            self.options.reserved_thread_per_task = 8
            self.options.openmp_thread_per_task = 1
            # We need 8 to get a proper binding. We reserve the node in exclusive mode anyway.
            self.options.gpu_per_node_count = 8
        else:
            self.options.reserved_thread_per_task = self.options.omp
            self.options.openmp_thread_per_task = self.options.omp
            self.options.gpu_per_node_count = 0

        if self.options.nodes:
            pass
        else:
            # 4 MI200~GPUs per node. Each GPU contains 2 GCD
            kMaxGPUPerNode = 4 * 2
            self.options.nodes = int(ceil(self.options.mpi / kMaxGPUPerNode))

    def run(self, arguments, dir):
        """
        Run a simulation
        """

        # Write the slurm script
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write(self.the_slurm_script.format(a_task_command=self.RUN_COMMAND,
                                                 a_task_command_arguments=arguments,
                                                 the_account=self.options.account,
                                                 the_partition=self.options.partition,
                                                 the_node_count=self.options.nodes,
                                                 the_mpi_process_count=self.options.mpi,
                                                 the_reserved_thread_count=self.options.reserved_thread_per_task,
                                                 the_gpu_per_node_count=self.options.gpu_per_node_count,
                                                 the_omp_thread_count=self.options.openmp_thread_per_task,
                                                 the_output_file=self.smilei_path.output_file,
                                                 the_maximum_task_duration=self.options.max_time))

        # Schedule the task(s)
        self.launch_job(self.RUN_COMMAND + ' ' + arguments, 
                        self.JOB, 
                        dir, 
                        self.options.max_time_seconds, 
                        self.smilei_path.output_file, 
                        repeat=2)
