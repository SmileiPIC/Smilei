from . import Machine

class MachineLLR(Machine):
    script = """\
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -q default
#PBS -j oe
#PBS -l walltime={max_time}
module purge
unset MODULEPATH
module use /opt/exp_soft/vo.llr.in2p3.fr/modulefiles_el7
module load hdf5/1.10.5-icc-omp4.1.1
module load python/3.7.0
module load h5py/hdf5_1.10.5-icc-omp4.1.1-py3.7.0
module load mpi4py/omp4.1.1-ib-icc_py3.7.0
module load compilers/gcc/9.x.x
module load fftw/3.3.10-omp-4.1.1-icc-19
export OMP_NUM_THREADS={omp}
export OMP_SCHEDULE=DYNAMIC
export KMP_AFFINITY=verbose
export PATH=$PATH:/opt/exp_soft/vo.llr.in2p3.fr/GALOP/beck
export LIBPXR=/home/llr/galop/derouil/applications.ompi216.Py3/picsar/lib
export LD_LIBRARY_PATH=$LIBPXR:$LD_LIBRARY_PATH
export FFTW3_LIB=/opt/exp_soft/vo.llr.in2p3.fr/fftw/3.3.10/opm-4.1.1-intel-19-el7/lib
export FFTW3_INC=/opt/exp_soft/vo.llr.in2p3.fr/fftw/3.3.10/opm-4.1.1-intel-19-el7/include
ulimit -s unlimited
#Specify the number of sockets per node in -mca orte_num_sockets
#Specify the number of cores per sockets in -mca orte_num_cores
cd {dir}
module list 2> module.log
{command}
echo $? > exit_status_file"""
    
    def __init__(self, smilei_path, options):
        from math import ceil
        
        self.smilei_path = smilei_path
        self.options = options
        
        self.MAKE = "make" + (" config=%s"%self.options.compile_mode if self.options.compile_mode else "")
        
        if self.options.nodes:
            self.NODES = self.options.nodes
            MPI_PER_SOCKET = int(ceil(self.options.mpi/(self.NODES*2)))
        else:
            self.NODES = int(ceil(self.options.mpi/2.))
            MPI_PER_SOCKET = 1
        self.ppn = {"jollyjumper":24, "tornado":36}[self.options.partition]
        
        if self.options.partition == "jollyjumper":
            self.JOB = "PBS_DEFAULT=llrlsi-jj.in2p3.fr qsub  "+self.smilei_path.exec_script
        elif self.options.partition == "tornado":
            self.JOB = "PBS_DEFAULT=poltrnd.in2p3.fr qsub  "+self.smilei_path.exec_script
        
        self.COMPILE_COMMAND = self.MAKE+' -j '+str(self.ppn)+' > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        self.COMPILE_TOOLS_COMMAND = 'make tables > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        self.CLEAN_COMMAND = 'make clean > /dev/null 2>&1'
        self.RUN_COMMAND = "mpirun --mca mpi_warn_on_fork 0 -mca orte_num_sockets 2 -mca orte_num_cores "+str(self.ppn) + " -map-by ppr:"+str(MPI_PER_SOCKET)+":socket:"+"pe="+str(self.options.omp) + " -n "+str(self.options.mpi)+" -x OMP_NUM_THREADS -x OMP_SCHEDULE "+self.smilei_path.workdirs+"smilei %s >"+self.smilei_path.output_file+" 2>&1"
    
    
    def compile(self, dir):
        """
        Compile Smilei
        """
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write( self.script.format(command=self.COMPILE_COMMAND, nodes=self.NODES, ppn=self.ppn, max_time=self.options.max_time, omp=self.options.omp, dir=dir) )
        
        self.launch_job(self.COMPILE_COMMAND, self.JOB, dir, self.options.max_time_seconds, self.smilei_path.COMPILE_ERRORS, repeat=2)
    
    
    def run(self, arguments, dir):
        """
        Run a simulation
        """
        command = self.RUN_COMMAND % arguments
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write( self.script.format(command=command, nodes=self.NODES, ppn=self.ppn, max_time=self.options.max_time, omp=self.options.omp, dir=dir) )
        
        self.launch_job(command, self.JOB, dir, self.options.max_time_seconds, self.smilei_path.output_file, repeat=2)
    