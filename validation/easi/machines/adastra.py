from . import Machine


class MachineAdastra(Machine):
    """
    As of 17 of march 22, this class, while named Adastra, instead targets the 
    Adastra porting machines. The real Adastra environment should not change much.
    """

    the_slurm_script = """#!/bin/bash
#SBATCH --job-name=smilei_validation
#SBATCH --nodes={the_node_count}               # Number of nodes
#SBATCH --ntasks={the_mpi_process_count}       # Number of MPI ranks
#SBATCH --cpus-per-task={the_omp_thread_count} # Number of cores per MPI rank 
#SBATCH --gpus-per-task={the_gpu_count}        # Number of cores per MPI rank 
#SBATCH --output=output
#SBATCH --error=error
#SBATCH --time={the_maximum_task_duration}

# TODO(Etienne M): dunno the partition/feature/constraints for adastra yet
# # --feature=MI200
# #SBATCH --gpus=mi100:1 or mi200:1

echo "Date              = $(date)";
echo "Hostname          = $(hostname -s)";
echo "Working Directory = $(pwd)";
echo "";
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES";
echo "Number of Tasks Allocated      = $SLURM_NTASKS";
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK";

# We should need only that to run the rest is loaded by default
module load rocm/4.5.0;

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
export OMP_SCHEDULE=dynamic;
# You may want to change "cores", to "threads". But hyperthreading, for an well 
# optimized apps is generally not something you want (ROB buffer should be full).
export OMP_PLACES=cores;

rocm-smi;
rocminfo;

export CRAY_ACC_DEBUG=0;


{a_task_command} > {the_output_file} 2>&1;

kRETVAL=$?;

# Put the result in the slurm output file.
cat {the_output_file};

echo "The task ended at = $(date)";

echo -n $kRETVAL > exit_status_file;
exit $kRETVAL;
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
        self.RUN_COMMAND = 'srun ' + self.smilei_path.workdirs + '/smilei {the_arguments}'
        
        if self.options.nodes:
            self.NODES = self.options.nodes
        else:
            self.NODES = int(ceil(self.options.mpi / 4))

    def run(self, arguments, dir):
        """
        Run a simulation
        """

        the_command = self.RUN_COMMAND.format(the_arguments=arguments)

        # Write the slurm script
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write(self.the_slurm_script.format(a_task_command=the_command, 
                                                 the_node_count=self.NODES, 
                                                 the_mpi_process_count=self.options.mpi, 
                                                 the_maximum_task_duration=self.options.max_time, 
                                                 the_omp_thread_count=self.options.omp,
                                                 the_output_file=self.smilei_path.output_file,
                                                 the_gpu_count='1' if 'gpu_amd' in self.options.compile_mode else '0'))

        # Schedule the task(s)
        self.launch_job(the_command, 
                        self.JOB, 
                        dir, 
                        self.options.max_time_seconds, 
                        self.smilei_path.output_file, 
                        repeat=2)
