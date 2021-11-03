from . import Machine

class MachineIrene(Machine):
    
    environments = {
        "skylake":"""
module purge
module load intel/20.0.4
module load mpi/intelmpi/20.0.4
module load flavor/hdf5/parallel hdf5/1.8.20
module load python3
export PATH=${HDF5_ROOT}/bin:${PATH}
export LD_LIBRARY_PATH=${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}
export HDF5_ROOT_DIR=${HDF5_ROOT}""",
        "a64fx":"""
module purge
module load fujitsu mpi
module load python3/3.8.10
export HDF5_ROOT=/ccc/work/cont003/mds/lobetmat/Libraries/hdf5-1.12_fujitsu_a64fx/install/
export PATH=${HDF5_ROOT}/bin:${PATH}
export LD_LIBRARY_PATH=${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}
export HDF5_ROOT_DIR=${HDF5_ROOT}
export LD_LIBRARY_PATH=/ccc/products/ucx-1.10.1/system/default/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/ccc/products2/python3-3.8.10/Rhel_8__aarch64-a64fx/system/default/install_tree/python/3.8.10/s2azw3pgbfzhfcf44tvnh652pju2vtyj/lib:/ccc/products2/python3-3.8.10/Rhel_8__aarch64-a64fx/system/default/install_tree/gettext/0.21/i4aacmkl6pqxumiqfa36455yfhoojidl/lib"""
    }
    
    compilation = """\
#!/bin/bash
#MSUB --job-name=smilei
#MSUB -N 1
#MSUB -n 48
#MSUB --output=output
#MSUB --error=error
#MSUB -q {partition}
#MSUB --time={max_time}
#MSUB -A {account}
#MSUB -m work,scratch
{env}
set -x
cd {dir}
cd ${BRIDGE_MSUB_PWD} 
module list 2> module.log
{command}
echo $? > exit_status_file"""
    
    script = """\
#!/bin/bash
#MSUB --job-name=smilei
#MSUB -N {NODES}
#MSUB -n {mpi}
#MSUB -c {omp}
#MSUB --output=output
#MSUB --error=error
#MSUB -q {partition}
#MSUB --time={max_time}
#MSUB -A {account}
#MSUB -m work,scratch
{env}
export OMP_NUM_THREADS={omp}
export OMP_SCHEDULE=DYNAMIC 
export OMP_PROC_BIND=true 
export KMP_AFFINITY=verbose 
set -x
ulimit -s unlimited 
cd ${BRIDGE_MSUB_PWD} 
cd "+dir+" 
module list 2> module.log
{command}
echo $? > exit_status_file"""
    
    def __init__(self, smilei_path, options):
        from math import ceil
        
        self.smilei_path = smilei_path
        self.options = options
        
        self.MAKE = "make" + (" config=%s"%self.options.compile_mode if self.options.compile_mode else "")
        
        self.JOB = "ccc_msub  "+self.smilei_path.exec_script
        
        self.COMPILE_COMMAND = self.MAKE+' -j4 > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        # self.COMPILE_TOOLS_COMMAND = 'make tables > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        self.CLEAN_COMMAND = 'make clean > /dev/null 2>&1'
        self.RUN_COMMAND = "export OMP_NUM_THREADS="+str(self.options.omp)+"; mpirun -np "+str(self.options.mpi)+" "+self.smilei_path.workdirs+"smilei %s >"+self.smilei_path.output_file
        
        self.env = environments[self.options.partition]
    
    def compile(self, dir):
        """
        Compile Smilei
        """
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write( compilation.format(command=self.COMPILE_COMMAND, env=self.env, account=self.options.account, max_time=self.options.max_time, dir=dir) )
        
        self.launch_job(self.COMPILE_COMMAND, self.JOB, dir, self.options.max_time_seconds, self.smilei_path.output_file, repeat=2)
    
    
    def run(self, arguments, dir):
        """
        Run a simulation
        """
        from math import ceil
        command = self.RUN_COMMAND % arguments
        NODES = int(ceil(self.options.mpi/2.))
        ppn = 24
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write( script.format(command=command, env=self.env, account=self.options.account, nodes=NODES, ppn=ppn, max_time=self.options.max_time, mpi=self.options.mpi, omp=self.options.omp, dir=dir) )
        
        self.launch_job(command, self.JOB, dir, self.options.max_time_seconds, self.smilei_path.output_file, repeat=2)
