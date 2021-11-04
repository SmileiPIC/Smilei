class Machine(object):
    
    def __init__(self, smilei_path, options):
        from subprocess import check_output
        import re
        
        self.smilei_path = smilei_path
        self.options = options
        
        self.MAKE = "make" + (" config=%s"%self.options.compile_mode if self.options.compile_mode else "")
        
        mpi_version = str(check_output("mpirun --version", shell=True))
        if re.search("Open MPI", mpi_version, re.I):
            v = re.search("\d\d?\.\d\d?\.\d\d?", mpi_version).group() # Full version number
            v = int(v.split(".")[0]) # Major version number
            if v > 1:
                MPIRUN = "mpirun --oversubscribe -np "
            else:
                MPIRUN = "mpirun -mca btl tcp,sm,self -np "
        else:
            MPIRUN = "mpirun -np "
        
        self.COMPILE_COMMAND = self.MAKE+' -j4 > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        # self.COMPILE_TOOLS_COMMAND = 'make tables > '+self.smilei_path.COMPILE_OUT+' 2>'+self.smilei_path.COMPILE_ERRORS
        self.CLEAN_COMMAND = 'make clean > /dev/null 2>&1'
        self.RUN_COMMAND = "export OMP_NUM_THREADS="+str(self.options.omp)+"; "+MPIRUN+str(self.options.mpi)+" "+self.smilei_path.workdirs+"smilei %s >"+self.smilei_path.output_file
    
    
    def clean(self):
        from subprocess import call
        
        call(self.CLEAN_COMMAND , shell=True)
    
    
    def compile(self, dir):
        """
        Compile Smilei
        """
        self.shell( self.COMPILE_COMMAND )
    
    
    def run(self, arguments, dir):
        """
        Run a simulation
        """
        self.shell( self.RUN_COMMAND % arguments )
    
    
    def shell(self, command):
        """
        Run a command directly in the shell
        """
        from subprocess import check_call, CalledProcessError
        from sys import exit
        
        try:
            if self.options.verbose:
                print("")
                print("Trying command `"+command+"`")
            check_call(command, shell=True)
        except CalledProcessError:
            if  self.options.verbose:
                print("")
                print("Execution failed for command `"+command+"`")
            exit(2)
    
    
    def launch_job(self, base_command, sub_command, dir, max_time_seconds, error_file, repeat=1):
        """
        Submit/launch a job based on the given options
        """
        from subprocess import check_call, CalledProcessError
        from os import sep
        from sys import exit
        from time import sleep
        
        # Make a file containing temporary exit status
        EXIT_STATUS = "100"
        with open(dir+sep+"exit_status_file", "w") as f:
            f.write(str(EXIT_STATUS))
        
        # Run the job several times if requested
        for n in range(repeat):
            try:
                check_call(sub_command, shell=True)
                break
            except CalledProcessError:
                if self.options.verbose:
                    print("")
                    print("Command failed #%d: `%s`"%(n, sub_command))
                    if n < repeat:
                        print("Wait and retry")
                sleep(10)
        # Exit if unsuccesful
        else:
            if self.options.verbose:
                print("Exit")
            exit(2)
        
        # Otherwise job is running
        if self.options.verbose:
            print("")
            print("Submitted job with command `"+base_command+"`")
            print("\tmax duration: %d s"%max_time_seconds)
        # Wait for the exit status to be set
        # current_time = 0
        while EXIT_STATUS == "100":# and current_time < options['max_time_seconds']):
            sleep(5)
            # current_time += 5
            with open(dir+sep+"exit_status_file", "r+") as f:
                EXIT_STATUS = f.readline()
            # if current_time > options['max_time_seconds']:
            #     print("Max time exceeded for command `"+command+"`")
            #     exit(2)
        # Check that the run succeeded
        if int(EXIT_STATUS) != 0:
            if self.options.verbose:
                print("")
                print("Execution failed for command `"+base_command+"`")
                COMMAND = "cat "+error_file
                try:
                    check_call(COMMAND, shell=True)
                except CalledProcessError:
                    print("")
                    print("Failed to print file `%s`"%error_file)
            exit(2)
