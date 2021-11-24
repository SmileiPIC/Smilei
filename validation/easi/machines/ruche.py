from . import Machine

class MachineRuche(Machine):
    script = """\
#!/bin/bash
#SBATCH --job-name=smilei
#SBATCH --nodes={NODES}
#SBATCH --ntasks="+{mpi}
#SBATCH --cpus-per-task={omp}
#SBATCH --output=output
#SBATCH --error=error
#SBATCH --time={max_time}
#SBATCH --partition=cpu_short
module purge
module load zlib/1.2.9/gcc-9.2.0
module load anaconda3/2020.02/gcc-9.2.0
module load intel/19.0.3/gcc-4.8.5
module load intel-mpi/2019.3.199/intel-19.0.3.199
module load hdf5/1.10.6/intel-19.0.3.199-intel-mpi
export OMP_NUM_THREADS={omp}
export OMP_SCHEDULE=DYNAMIC
export OMP_PLACES=cores
export KMP_AFFINITY=verbose
export I_MPI_CXX=icpc
export HDF5_ROOT_DIR=/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/intel-19.0.3.199/hdf5-1.10.6-na3ilncuwbx2pdim2xaqwf23sgqza6qo
ulimit -s unlimited
cd ${SLURM_SUBMIT_DIR}
#Specify the number of sockets per node in -mca orte_num_sockets
#Specify the number of cores per sockets in -mca orte_num_cores
cd {dir}
module list 2> module.log
{command}
echo $? > exit_status_file"""
    
    def __init__(self, smilei_path, options):
        from math import ceil
        from os import exit
        
        self.smilei_path = smilei_path
        self.options = options
        
        self.MAKE = "make" + (" config=%s"%self.options.compile_mode if self.options.compile_mode else "")
        
        self.JOB = "sbatch "+self.smilei_path.exec_script
        
        self.ppn = 20
        if ppn < self.options.omp :
            print("Smilei cannot be run with "+str(self.options.omp)+" threads on Ruche and partition "+self.options.partition)
            exit(4)
        if self.options.nodes:
            self.NODES = self.options.nodes
            MPI_PER_SOCKET = int(ceil(self.options.mpi/(self.NODES*2)))
        else:
            self.NODES = int(ceil(self.options.mpi/2.))
            MPI_PER_SOCKET = 1
        self.COMPILE_COMMAND = self.MAKE+' -j 20 machine=ruche > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        # self.COMPILE_TOOLS_COMMAND = 'make tables > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        self.CLEAN_COMMAND = 'make clean > /dev/null 2>&1'
        self.RUN_COMMAND ="srun "+self.smilei_path.workdirs+"smilei %s >"+self.smilei_path.output_file+" 2>&1"
    
    
    def compile(self, dir):
        """
        Compile Smilei
        """
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write( self.script.format(command=self.COMPILE_COMMAND, nodes=self.NODES, mpi=self.options.mpi, max_time=self.options.max_time, omp=self.options.omp, dir=dir) )
        
        self.launch_job(self.COMPILE_COMMAND, self.JOB, dir, self.options.max_time_seconds, self.smilei_path.output_file, repeat=2)
    
    
    def run(self, arguments, dir):
        """
        Run a simulation
        """
        command = self.RUN_COMMAND % arguments
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write( self.script.format(command=command, nodes=self.NODES, mpi=self.options.mpi, max_time=self.options.max_time, omp=self.options.omp, dir=dir) )
        
        self.launch_job(command, self.JOB, dir, self.options.max_time_seconds, self.smilei_path.output_file, repeat=2)
    
