# --------------------
# Launch job function
# --------------------

from subprocess import check_call

def launch_job(base_command, job_command, dir, max_time, output, repeat=1):
    """
    Submit/launch a job based on the given options
    """
    # Make a file containing temporary exit status
    EXIT_STATUS = "100"
    with open(dir+s+"exit_status_file", "w") as f:
        f.write(str(EXIT_STATUS))
    # Run the job several times if requested
    for n in range(repeat):
        try:
            check_call(job_command, shell=True)
            break
        except CalledProcessError:
            if options['verbose']:
                print()
                print("Command failed #%d: `%s`"%(n, job_command))
                if n < repeat:
                    print("Wait and retry")
            sleep(10)
    # Exit if unsuccesful
    else:
        if options['verbose']:
            print("Exit")
        sys.exit(2)
    # Otherwise job is running
    if options['verbose']:
        print()
        print("Submitted job with command `"+base_command+"`")
        print("\tmax duration: %d s"%max_time)
    # Wait for the exit status to be set
    # current_time = 0
    while EXIT_STATUS == "100":# and current_time < options['max_time_seconds']):
        sleep(5)
        # current_time += 5
        with open(dir+s+"exit_status_file", "r+") as f:
            EXIT_STATUS = f.readline()
        # if current_time > options['max_time_seconds']:
        #     print("Max time exceeded for command `"+command+"`")
        #     sys.exit(2)
    # Check that the run succeeded
    if int(EXIT_STATUS) != 0:
        if options['verbose']:
            print()
            print("Execution failed for command `"+base_command+"`")
            COMMAND = "/bin/bash cat "+output
            try:
                check_call(COMMAND, shell=True)
            except CalledProcessError:
                print()
                print("Failed to print file `%s`"%output)
        sys.exit(2)
