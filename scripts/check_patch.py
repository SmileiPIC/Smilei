#This script checks the number of patches in patch_load.txt.
#It prints the difference between the total number of patch and the theoretical number of patch.

import scipy
import matplotlib.pyplot as plt

#   USer defined parameters
nmpi = 400
npatchx = 256
npatchy = 128
npatchz = 128
###########################
filename = "patch_load.txt"

file = open(filename)

data = file.readlines()
noutput = len(data)/(nmpi+1)

for iout in range(noutput):
    data_first =data[(nmpi+1)*iout+1:(nmpi+1)*iout+nmpi+1] 
    for j in range(nmpi):
        data_first[j]= scipy.int32(data_first[j].split()[-1])
    
    arrfirst = scipy.array(data_first)
    print arrfirst.sum()-(npatchx*npatchy*npatchz)
