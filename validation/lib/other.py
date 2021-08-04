from launch_job import *

def run_other(command, dir, mode, options, parameters):
    try:
        if options['verbose']:
            print()
            print("Trying command `"+command+"`")
        check_call(command, shell=True)
    except CalledProcessError:
        if options['verbose']:
            print()
            print("Execution failed for command `"+command+"`")
        sys.exit(2)
