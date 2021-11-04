
class SmileiPath(object):
    def __init__(self):
        from os import environ, path, sep
        from inspect import stack
        
        if "self.smilei_path.root" in environ :
            self.root = environ["self.smilei_path.root"]+sep
        else:
            self.root = path.dirname(path.abspath(stack()[0][1]))+sep+".."+sep+".."+sep
        self.root = path.abspath(self.root)+sep
        self.benchmarks = self.root+"benchmarks"+sep
        self.scrips = self.root+"scripts"+sep
        self.validation = self.root+"validation"+sep
        self.references = self.validation+"references"+sep
        self.analyses = self.validation+"analyses"+sep
        
        self.workdirs = self.root+"validation"+sep+"workdirs"+sep
        
        self.SMILEI_TOOLS_W = self.workdirs+"smilei_tables"
        self.SMILEI_TOOLS_R = self.root+"smilei_tables"
        
        self.COMPILE_ERRORS = self.workdirs+'compilation_errors'
        self.COMPILE_OUT = self.workdirs+'compilation_out'
        
        self.exec_script = 'exec_script.sh'
        self.exec_script_output = 'exec_script.out'
        self.output_file = 'smilei_exe.out'


class ValidationOptions(object):
    def __init__(self, **kwargs):
        # Get general parameters from kwargs
        self.verbose      = kwargs.pop( "verbose"     , False         )
        self.compile_only = kwargs.pop( "compile_only", False         )
        self.bench        = kwargs.pop( "bench"       , ""            )
        self.omp          = kwargs.pop( "omp"         , 12            )
        self.mpi          = kwargs.pop( "mpi"         , 4             )
        self.generate     = kwargs.pop( "generate"    , False         )
        self.showdiff     = kwargs.pop( "showdiff"    , False         )
        self.nb_restarts  = kwargs.pop( "nb_restarts" , 0             )
        self.max_time     = kwargs.pop( "max_time"    , "00:10:00"    )
        self.compile_mode = kwargs.pop( "compile_mode", ""            )
        self.log          = kwargs.pop( "log"         , ""            )
        self.partition    = kwargs.pop( "partition"   , "jollyjumper" )
        self.account      = kwargs.pop( "account"     , ""            )
        
        if kwargs:
            raise Exception("Unknown options for validation: "+", ".join(kwargs))
        
        from numpy import array, sum
        self.max_time_seconds = sum(array(self.max_time.split(":"),dtype=int)*array([3600,60,1]))



def loadReference(references_path, bench_name):
    import pickle
    from sys import exit
    try:
        try:
            with open(references_path + bench_name + ".txt", 'rb') as f:
                return pickle.load(f, fix_imports=True, encoding='latin1')
        except:
            with open(references_path + bench_name + ".txt", 'r') as f:
                return pickle.load(f)
    except:
        print("Unable to find the reference data for "+bench_name)
        exit(1)

def matchesWithReference(data, expected_data, data_name, precision, error_type="absolute_error"):
    from numpy import array, double, abs, unravel_index, argmax, all, flatnonzero, isnan
    # ok if exactly equal (including strings or lists of strings)
    try:
        if expected_data == data:
            return True
    except:
        pass
    # If numbers:
    try:
        double_data = array(double(data), ndmin=1)
        if precision is not None:
            error = abs( double_data-array(double(expected_data), ndmin=1) )
            if error_type == "absolute_error":
                pass
            elif error_type == "relative_error":
                try:
                    error /= double_data
                    if isnan( error ).any():
                        raise
                except Exception as e:
                    print( "Error in comparing with reference: division by zero (relative error)" )
                    return False
            else:
                print( "Unknown error_type = `"+error_type+"`" )
                return False
            if type(precision) in [int, float]:
                max_error_location = unravel_index(argmax(error), error.shape)
                max_error = error[max_error_location]
                if max_error < precision:
                    return True
                print("Reference quantity '"+data_name+"' does not match the data (required precision "+str(precision)+")")
                print("Max error = "+str(max_error)+" at index "+str(max_error_location))
            else:
                try:
                    precision = array(precision)
                    if (error <= precision).all():
                        return True
                    print("Reference quantity '"+data_name+"' does not match the data")
                    print("Error = ")
                    print(error)
                    print("Precision = ")
                    print(precision)
                    print("Failure at indices "+", ".join([str(a) for a in flatnonzero(error > precision)]))
                except Exception as e:
                    print( "Error with requested precision (of type %s). Cannot be compared to the data (of type %s)"%(type(precision), type(error)) )
                    return False
        else:
            if all(double_data == double(expected_data)):
                return True
            print("Reference quantity '"+data_name+"' does not match the data")
    except Exception as e:
        print("Reference quantity '"+data_name+"': unable to compare to data")
        print( e )
    return False


class Validation(object):
    
    _dataNotMatching = False
    
    def __init__(self, **kwargs):
        # Obtain options
        self.options = ValidationOptions(**kwargs)
        
        # Find smilei folders
        self.smilei_path = SmileiPath()
        
        # Get the current version of Smilei
        from os import environ
        from subprocess import check_output
        self.git_version = check_output( "cd "+self.smilei_path.root+" && echo `git log -n 1 --format=%h`-", shell=True ).decode()[:-1]
        if 'CI_COMMIT_BRANCH' in environ:
            self.git_version += environ['CI_COMMIT_BRANCH']
        else:
            self.git_version += check_output("cd "+self.smilei_path.root+" && echo `git rev-parse --abbrev-ref HEAD`", shell=True ).decode()[:-1]
        
        # Define commands depending on host
        from socket import gethostname
        from .machines import Machine, MachineLLR, MachinePoincare, MachineRuche, MachineIrene
        self.HOSTNAME = gethostname()
        if "llrlsi-gw" in self.HOSTNAME:
            self.machine = MachineLLR( self.smilei_path, self.options )
        elif "poincare" in self.HOSTNAME:
            self.machine = MachinePoincare( self.smilei_path, self.options )
        elif "ruche" in self.HOSTNAME:
            self.machine = MachineRuche( self.smilei_path, self.options )
        elif "irene" in self.HOSTNAME:
            self.machine = MachineIrene( self.smilei_path, self.options )
        else:
            self.machine = Machine( self.smilei_path, self.options )
    
    def compile(self):
        from sys import exit
        from os import chdir, sep, stat, remove, rename
        from os.path import exists
        from .tools import mkdir, date, date_string
        from shutil import copy2
        from subprocess import CalledProcessError
        
        if self.options.verbose:
            print("---------------------------")
            print("Compiling Smilei")
            print("---------------------------")
        
        SMILEI_W = self.smilei_path.workdirs + "smilei"
        SMILEI_R = self.smilei_path.root + "smilei"
        
        # Get state of smilei bin in root folder
        chdir(self.smilei_path.root)
        STAT_SMILEI_R_OLD = stat(SMILEI_R) if exists(SMILEI_R) else ' '
        
        # CLEAN
        # If no smilei bin in the workdir, or it is older than the one in smilei directory,
        # clean to force compilation
        mkdir(self.smilei_path.workdirs)
        if exists(SMILEI_R) and (not exists(SMILEI_W) or date(SMILEI_W)<date(SMILEI_R)):
            self.machine.clean()
        
        def workdir_archiv() :
            # Creates an archives of the workdir directory
            exe_path = self.smilei_path.workdirs+"smilei"
            if exists(exe_path):
                ARCH_WORKDIR = self.smilei_path.workdirs+'_'+date_string(exe_path)
                rename(self.smilei_path.workdirs, ARCH_WORKDIR)
                mkdir(self.smilei_path.workdirs)
        
        try:
            # Remove the compiling errors files
            if exists(self.smilei_path.COMPILE_ERRORS):
                remove(self.smilei_path.COMPILE_ERRORS)
            # Compile
            self.machine.compile( self.smilei_path.root )
            if STAT_SMILEI_R_OLD!=stat(SMILEI_R) or date(SMILEI_W)<date(SMILEI_R): # or date(SMILEI_TOOLS_W)<date(SMILEI_TOOLS_R) :
                # if new bin, archive the workdir (if it contains a smilei bin)
                # and create a new one with new smilei and compilation_out inside
                if exists(SMILEI_W): # and path.exists(SMILEI_TOOLS_W):
                    workdir_archiv()
                copy2(SMILEI_R, SMILEI_W)
                #copy2(SMILEI_TOOLS_R,SMILEI_TOOLS_W)
                if self.options.verbose:
                    print("Smilei compilation succeed.")
            else:
                if self.options.verbose:
                    print("Smilei compilation not needed.")
        
        except CalledProcessError as e:
            # if compiling errors, archive the workdir (if it contains a smilei bin),
            # create a new one with compilation_errors inside and exit with error code
            workdir_archiv()
            if self.options.verbose:
                print("Smilei compilation failed. " + str(e.returncode))
            exit(3)
        
        if self.options.verbose:
            print()
    
    
    def run_all(self):
        from sys import exit
        from os import sep, chdir, getcwd
        from os.path import basename, splitext, exists, isabs
        from shutil import rmtree
        from .tools import mkdir, execfile
        from .log import Log
        import re, sys
        
        # Load the happi module
        sys.path.insert(0, self.smilei_path.root)
        import happi
        
        INITIAL_DIRECTORY = getcwd()
        
        _dataNotMatching = False
        for BENCH in self.list_benchmarks():
            
            SMILEI_BENCH = self.smilei_path.benchmarks + BENCH
            
            # Create the workdir path
            WORKDIR = self.smilei_path.workdirs + 'wd_'+basename(splitext(BENCH)[0]) + sep
            mkdir(WORKDIR)
            WORKDIR += str(self.options.mpi) + sep
            mkdir(WORKDIR)
            WORKDIR += str(self.options.omp) + sep
            mkdir(WORKDIR)
            
            # If there are restarts, prepare a Checkpoints block in the namelist
            RESTART_INFO = ""
            if self.options.nb_restarts > 0:
                # Load the namelist
                namelist = happi.openNamelist(SMILEI_BENCH)
                niter = namelist.Main.simulation_time / namelist.Main.timestep
                # If the simulation does not have enough timesteps, change the number of restarts
                if self.options.nb_restarts > niter - 4:
                    self.options.nb_restarts = max(0, niter - 4)
                    if self.options.verbose:
                        print("Not enough timesteps for restarts. Changed to "+str(self.options.nb_restarts)+" restarts")
                if self.options.nb_restarts > 0:
                    # Find out the optimal dump_step
                    dump_step = int( (niter+3.) / (self.options.nb_restarts+1) )
                    # Prepare block
                    if len(namelist.Checkpoints) > 0:
                        RESTART_INFO = (" \""
                            + "Checkpoints.keep_n_dumps="+str(self.options.nb_restarts)+";"
                            + "Checkpoints.dump_minutes=0.;"
                            + "Checkpoints.dump_step="+str(dump_step)+";"
                            + "Checkpoints.exit_after_dump=True;"
                            + "Checkpoints.restart_dir=%s;"
                            + "\""
                        )
                    else:
                        RESTART_INFO = (" \"Checkpoints("
                            + "    keep_n_dumps="+str(self.options.nb_restarts)+","
                            + "    dump_minutes=0.,"
                            + "    dump_step="+str(dump_step)+","
                            + "    exit_after_dump=True,"
                            + "    restart_dir=%s,"
                            + ")\""
                        )
                del namelist
            
            # Prepare logging
            if self.options.log:
                log_dir = ("" if isabs(self.options.log) else INITIAL_DIRECTORY + sep) + self.options.log + sep
                log = Log(log_dir, log_dir + BENCH + ".log")
            
            # Loop restarts
            for irestart in range(self.options.nb_restarts+1):
                
                RESTART_WORKDIR = WORKDIR + "restart%03d"%irestart + sep

                execution = True
                if not exists(RESTART_WORKDIR):
                    mkdir(RESTART_WORKDIR)
                elif self.options.generate:
                    execution = False

                chdir(RESTART_WORKDIR)

                # Copy of the databases
                # For the cases that need a database
                # if BENCH in [
                #         "tst1d_09_rad_electron_laser_collision.py",
                #         "tst1d_10_pair_electron_laser_collision.py",
                #         "tst2d_08_synchrotron_chi1.py",
                #         "tst2d_09_synchrotron_chi0.1.py",
                #         "tst2d_v_09_synchrotron_chi0.1.py",
                #         "tst2d_v_10_multiphoton_Breit_Wheeler.py",
                #         "tst2d_10_multiphoton_Breit_Wheeler.py",
                #         "tst2d_15_qed_cascade_particle_merging.py",
                #         "tst3d_15_magnetic_shower_particle_merging.py"
                #     ]:
                #     try :
                #         # Copy the database
                #         check_call(['cp '+SMILEI_DATABASE+'/*.h5 '+RESTART_WORKDIR], shell=True)
                #     except CalledProcessError:
                #         if self.options.verbose :
                #             print(  "Execution failed to copy databases in ",RESTART_WORKDIR)
                #         sys.exit(2)
                
                # If there are restarts, adds the Checkpoints block
                arguments = SMILEI_BENCH
                if self.options.nb_restarts > 0:
                    if irestart == 0:
                        RESTART_DIR = "None"
                    else:
                        RESTART_DIR = "'"+WORKDIR+("restart%03d"%(irestart-1))+sep+"'"
                    arguments += RESTART_INFO % RESTART_DIR
                
                # Run smilei
                if execution:
                    if self.options.verbose:
                        print("---------------------------")
                        print('Running '+BENCH+' on '+self.HOSTNAME+' with '+str(self.options.omp)+'x'+str(self.options.mpi)+' OMPxMPI' + ((", restart #"+str(irestart)) if irestart>0 else ""))
                        print("---------------------------")
                    self.machine.run( arguments, RESTART_WORKDIR )
                
                # Check the output for errors
                errors = []
                search_error = re.compile('error', re.IGNORECASE)
                with open(self.smilei_path.output_file,"r") as fout:
                    errors = [line for line in fout if search_error.search(line)]
                if errors:
                    if self.options.verbose:
                        print("")
                        print("Errors appeared while running the simulation:")
                        print("---------------------------------------------")
                        for error in errors:
                            print(error)
                    exit(2)
                
                # Scan some info for logging
                if self.options.log:
                    log.scan(self.smilei_path.output_file)
            
            # Append info in log file
            if self.options.log:
                log.append(git_version)
            
            # Find the validation script for this bench
            validation_script = self.smilei_path.analyses + "validate_" + BENCH
            if self.options.verbose:
                print("")
            if not exists(validation_script):
                print("Unable to find the validation script "+validation_script)
                exit(1)
            
            chdir(WORKDIR)
            
            # If required, generate the references
            if self.options.generate:
                if self.options.verbose:
                    print( '----------------------------------------------------')
                    print( 'Generating reference for '+BENCH)
                    print( '----------------------------------------------------')
                Validate = self.CreateReference(self.smilei_path.references, BENCH)
                execfile(validation_script, globals(), {"Validate":Validate})
                Validate.write()
            
            # Or plot differences with respect to existing references
            elif self.options.showdiff:
                if self.options.verbose:
                    print( '----------------------------------------------------')
                    print( 'Viewing differences for '+BENCH)
                    print( '----------------------------------------------------')
                Validate = self.ShowDiffWithReference(self.smilei_path.references, BENCH)
                execfile(validation_script, globals(), {"Validate":Validate})
                if _dataNotMatching:
                    print("Benchmark "+BENCH+" did NOT pass")
            
            # Otherwise, compare to the existing references
            else:
                if self.options.verbose:
                    print( '----------------------------------------------------')
                    print( 'Validating '+BENCH)
                    print( '----------------------------------------------------')
                Validate = self.CompareToReference(self.smilei_path.references, BENCH)
                execfile(validation_script, globals(), {"Validate":Validate})
                if _dataNotMatching:
                    chdir(INITIAL_DIRECTORY)
                    exit(1)
            
            # Clean workdirs, goes here only if succeeded
            chdir(self.smilei_path.workdirs)
            rmtree(WORKDIR, True)
            if self.options.verbose:
                print( "")
        
        if _dataNotMatching:
            print( "Errors detected")
        else:
            print( "Everything passed")
        chdir(INITIAL_DIRECTORY)
    
    
    def list_benchmarks(self):
        from os.path import basename
        from glob import glob
        
        # Build the list of the requested input files
        list_validation = [basename(b) for b in glob(self.smilei_path.analyses+"validate_tst*py")]
        if self.options.bench == "":
            benchmarks = [basename(b) for b in glob(self.smilei_path.benchmarks+"tst*py")]
        else:
            benchmarks = glob( self.smilei_path.benchmarks + self.options.bench )
            benchmarks = [b.replace(self.smilei_path.benchmarks,'') for b in benchmarks]
        benchmarks = [b for b in benchmarks if "validate_"+b in list_validation]
        if not benchmarks:
            raise Exception("Input file(s) "+self.options.bench+" not found, or without validation file")
        
        if self.options.verbose:
            print("")
            print("The list of input files to be validated is:\n\t"+"\n\t".join(benchmarks))
            print("")
        
        return benchmarks
    
    
    # DEFINE A CLASS TO CREATE A REFERENCE
    class CreateReference(object):
        def __init__(self, references_path, bench_name):
            self.reference_file = references_path+bench_name+".txt"
            self.data = {}
        
        def __call__(self, data_name, data, precision=None, error_type="absolute_error"):
            self.data[data_name] = data
        
        def write(self):
            import pickle
            from os.path import getsize
            from os import remove
            from sys import exit
            with open(reference_file, "wb") as f:
                pickle.dump(data, f, protocol=2)
            size = getsize(reference_file)
            if size > 1000000:
                print("Reference file is too large ("+str(size)+"B) - suppressing ...")
                remove(reference_file)
            elif self.options.verbose:
                print("Created reference file "+reference_file)

    # DEFINE A CLASS TO COMPARE A SIMULATION TO A REFERENCE
    class CompareToReference(object):
        def __init__(self, references_path, bench_name):
            self.ref_data = loadReference(references_path, bench_name)
        
        def __call__(self, data_name, data, precision=None, error_type="absolute_error"):
            from sys import exit
            # verify the name is in the reference
            if data_name not in self.ref_data.keys():
                print("Reference quantity '"+data_name+"' not found")
                _dataNotMatching = True
                return
            expected_data = self.ref_data[data_name]
            if not matchesWithReference(data, expected_data, data_name, precision, error_type):
                print("Reference data:")
                print(expected_data)
                print("New data:")
                print(data)
                print()
                _dataNotMatching = True

    # DEFINE A CLASS TO VIEW DIFFERENCES BETWEEN A SIMULATION AND A REFERENCE
    class ShowDiffWithReference(object):
        def __init__(self, references_path, bench_name):
            self.ref_data = loadReference(references_path, bench_name)
        
        def __call__(self, data_name, data, precision=None, error_type="absolute_error"):
            import matplotlib.pyplot as plt
            from numpy import array
            plt.ion()
            print("Showing differences about '"+data_name+"'")
            print("--------------------------")
            # verify the name is in the reference
            if data_name not in self.ref_data.keys():
                print("\tReference quantity not found")
                expected_data = None
            else:
                expected_data = self.ref_data[data_name]
            print_data = False
            # First, check whether the data matches
            if not matchesWithReference(data, expected_data, data_name, precision, error_type):
                _dataNotMatching = True
            # try to convert to array
            try:
                data_float = array(data, dtype=float)
                expected_data_float = array(expected_data, dtype=float)
            # Otherwise, simply print the result
            except:
                print("\tQuantity cannot be plotted")
                print_data = True
                data_float = None
            # Manage array plotting
            if data_float is not None:
                if expected_data is not None and data_float.shape != expected_data_float.shape:
                    print("\tReference and new data do not have the same shape: "+str(expected_data_float.shape)+" vs. "+str(data_float.shape))
                if expected_data is not None and data_float.ndim != expected_data_float.ndim:
                    print("\tReference and new data do not have the same dimension: "+str(expected_data_float.ndim)+" vs. "+str(data_float.ndim))
                    print_data = True
                elif data_float.size == 0:
                    print("\t0D quantity cannot be plotted")
                    print_data = True
                elif data_float.ndim == 1:
                    nplots = 2
                    if expected_data is None or data_float.shape != expected_data_float.shape:
                        nplots = 1
                    fig = plt.figure()
                    fig.suptitle(data_name)
                    print("\tPlotting in figure "+str(fig.number))
                    ax1 = fig.add_subplot(nplots,1,1)
                    ax1.plot( data_float, label="new data" )
                    ax1.plot( expected_data_float, label="reference data" )
                    ax1.legend()
                    if nplots == 2:
                        ax2 = fig.add_subplot(nplots,1,2)
                        ax2.plot( data_float-expected_data_float )
                        ax2.set_title("difference")
                elif data_float.ndim == 2:
                    nplots = 3
                    if expected_data is None:
                        nplots = 1
                    elif data_float.shape != expected_data_float.shape:
                        nplots = 2
                    fig = plt.figure()
                    fig.suptitle(data_name)
                    print("\tPlotting in figure "+str(fig.number))
                    ax1 = fig.add_subplot(1,nplots,1)
                    im = ax1.imshow( data_float )
                    ax1.set_title("new data")
                    plt.colorbar(im)
                    if nplots > 1:
                        ax2 = fig.add_subplot(1,nplots,2)
                        im = ax2.imshow( expected_data_float )
                        ax2.set_title("reference data")
                        plt.colorbar( im )
                    if nplots > 2:
                        ax3 = fig.add_subplot(1,nplots,nplots)
                        im = ax3.imshow( data_float-expected_data_float )
                        ax3.set_title("difference")
                        plt.colorbar( im )
                    plt.draw()
                    plt.show()
                else:
                    print("\t"+str(data_float.ndim)+"D quantity cannot be plotted")
                    print_data = True
            # Print data if necessary
            if print_data:
                if expected_data is not None:
                    print("\tReference data:")
                    print(expected_data)
                print("\tNew data:")
                print(data)
    
    

