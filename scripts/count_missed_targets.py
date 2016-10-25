#Compute how many load balancing targets were missed.

import scipy
from pylab import *
directories = [""]
MPIloadfiles = "MPI_load.txt"


MPIloadfile = open(MPIloadfiles)
keepreading = True
total = 0
total_time = 0
nmpi = 0
time = 0
#Compute number of MPI ranks
for iskip in range(3):  #Skip 2 lines
    MPIbuffer = MPIloadfile.next()

while (MPIbuffer.split()[0] == "patch_count"):
    MPIbuffer = MPIloadfile.next()
    nmpi += 1

print nmpi, "MPI ranks were used in the simulation." 
#Skip 1 line
MPIbuffer = MPIloadfile.next()

#Compute number of missed target
while (keepreading):
    
    patch_count = scipy.zeros((nmpi),dtype='int64')
    patch_achieved = scipy.zeros((nmpi),dtype='int64')
    patch_target = scipy.zeros((nmpi),dtype='int64')

    for impi in range(nmpi):
        MPIbuffer = MPIloadfile.next()
        MPIbuffer = MPIbuffer.split()
        patch_achieved[impi] = int (MPIbuffer[2])
        patch_target[impi] = int (MPIbuffer[6])
    diff = scipy.where(patch_achieved != patch_target)[0]
    if len(diff) > 0:
        print "Target was missed at load balancing call number ", time
        print "MPI ranks which missed: ", diff
        total += len(diff)
        total_time += 1
    time += 1
    for iskip in range(2):  #Skip 2 lines
        try:
            MPIbuffer = MPIloadfile.next()
        except StopIteration:
            keepreading = False
MPIloadfile.close()
print "Target was missed ", total_time, " times for a total of ", total, " targets missed."
    
