from . import Machine

class MachinePoincare(Machine):
    script = """\
# environnement
module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu 2>&1 > /dev/null
unset LD_PRELOAD
export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0
# 
# execution
export OMP_NUM_THREADS={omp}
{command}
echo $? > exit_status_file 
exit $?"""
    
    def __init__(self, smilei_path, options):
        from math import ceil
        
        self.smilei_path = smilei_path
        self.options = options
        
        self.MAKE = "make" + (" config=%s"%self.options.compile_mode if self.options.compile_mode else "")
        
        self.JOB = "/bin/bash "+self.smilei_path.exec_script+" > "+self.smilei_path.exec_script_out+" 2>&1"
        
        self.COMPILE_COMMAND = self.MAKE+' -j 6 > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        self.COMPILE_TOOLS_COMMAND = 'make tables > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        self.CLEAN_COMMAND = 'module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0 gnu gnu ; unset LD_PRELOAD ; export PYTHONHOME=/gpfslocal/pub/python/anaconda/Anaconda-2.1.0 > /dev/null 2>&1;make clean > /dev/null 2>&1'
        self.RUN_COMMAND = "mpirun -np "+str(self.options['mpi'])+" "+self.smilei_path.wordirs+"smilei %s >"+self.smilei_path.output_file
    
    def compile(self, dir):
        """
        Compile Smilei
        """
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write( self.script.format(command=self.COMPILE_COMMAND, omp=self.options.omp) )
        
        self.launch_job(self.COMPILE_COMMAND, self.JOB, dir, self.options.max_time_seconds, self.smilei_path.output_file, repeat=2)
    
    
    def run(self, arguments, dir):
        """
        Run a simulation
        """
        command = self.RUN_COMMAND % arguments
        with open(self.smilei_path.exec_script, 'w') as f:
            f.write( self.script.format(command=command, omp=self.options.omp) )
            
        self.launch_job(command, self.JOB, dir, self.options.max_time_seconds, self.smilei_path.output_file, repeat=2)
    