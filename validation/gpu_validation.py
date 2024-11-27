from easi import Validation

import subprocess, sys

# TODO: get list of all scripts starting with gpu_*.py instead
list_test=["tst2d_01_plasma_mirror", "tst1d_03_thermal_expansion"]

# create workdir for each test that includes: slurm submission script, smilei executable and the inputfile

#import os
#path = os.getcwd()

for test_name in list_test:
    test_folder = "/gpfs/workdir/prouveurc/runs/ci_runs/" + test_name
    subprocess.run("mkdir " + test_folder, shell = True, executable="/bin/bash")
    subprocess.run("cp ../smilei " + test_folder + "/", shell = True, executable="/bin/bash")
    subprocess.run("cp ../benchmarks/gpu/" + "gpu_" + test_name + ".py " + test_folder + "/input.py", shell = True, executable="/bin/bash")
    subprocess.run("cp easi/machines/ruche_gpu.sh " + test_folder + "/", shell = True, executable="/bin/bash")
    subprocess.run("cd " + test_folder +" ; sbatch -W ruche_gpu.sh", shell = True, executable="/bin/bash")
    #subprocess.run("cd " + path, shell = True, executable="/bin/bash") useless

# there is a real question on what is the best strategy... wait after all jobs are submitted, maybe wait only for the last if we assume it will be the last to finish
# of wait after each job is finished. One is faster but assume the sleep time is a good approximation of the total time, the other makes sure all tests are done

# after launching all test cases, wait appropriate amount of time:
#import time
#time.sleep(600) 

# Validate the results with the reference results
V = Validation()
test_passed = True
for test_name in list_test:
    #copying test_name.py as temp.py to cp gpu_test_name.py test_name.py in order to use the same reference (to be changed once fusion in validation.py)
    #subprocess.run("cp analyses/validate_" + test_name + ".py analyses/temp.py" , shell = True, executable="/bin/bash")
    subprocess.run("cp analyses/gpu_validate_" + test_name + ".py analyses/validate_" + test_name + ".py" , shell = True, executable="/bin/bash")
    #time.sleep(5)

# WARINING: fusing this 2 for loops will result in an error on ruche, presumably because of filesystem delay (?)

for test_name in list_test:
    print(test_name)
    result_correct = V.compare( test_name + ".py", "/gpfs/workdir/prouveurc/runs/ci_runs/" + test_name)
    #subprocess.run("mv analyses/temp.py analyses/validate_" + test_name + ".py" , shell = True, executable="/bin/bash")
    # if need be, we can do a git restore here / or use the cp to  temporary file and back
    if not result_correct: #result!="PASS":
        print("ERROR FOR: " + test_name)
        test_passed = False
        break 

if (test_passed):
    print("CI GPU PASSED")

# clean up
subprocess.run("rm -rf /gpfs/workdir/prouveurc/runs/ci_runs/*", shell = True, executable="/bin/bash")

