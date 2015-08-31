#Plot the evolution of the time for 10 iterations. It measures the evolution of the load imbalance.

import scipy
from pylab import *

def treatline(b, timearray):
    b = b.split()
    if(len(b) > 7):
        if (b[6] == "sec") and (b[7] == "="):
            timearray.append(float(b[8]))

Ndpi = 200 #Quality and size of the plots. Set to 1000 for publication quality.

timepertenit=[]; it_list=[]

directories = ["/sps/beck/smilei/tst_patchexchange_noexchange/","/sps/beck/smilei/tst_patchexchange_withexchanges/"]
logfiles = ["sub_tst2d_accelec.sh.o4894","sub_tst2d_accelec.sh.o4895"]

for ifile in range(len(directories)):

    timearray = []
    datafile = open(directories[ifile]+logfiles[ifile])
    
    location = 0
    location_old = -1
    nline = 0
    
    while (location != location_old):
        location_old = location
        buffer = datafile.readline()
        nline = nline + 1
        location = datafile.tell()
        treatline(buffer, timearray)
    
    datafile.close()
    
    print nline, "lines were read, ", len(timearray), " times were collected."
    
    timearray = scipy.array(timearray)
    print "array timearray mean = ", timearray.mean() 

 
    timepertenit.append( timearray[1:]-timearray[:-1])
    it_list.append((scipy.arange(size(timepertenit[ifile]))+1)*10)
   
 
fig=figure(figsize=(20,10))
line0 = plot(it_list[0],timepertenit[0],lw=5)
line1 = plot(it_list[1],timepertenit[1],lw=5)
    
legend( (line0, line1), ('No exchange', 'With exchanges'), loc='upper left')
xlabel("Number of iterations",fontsize=35)
ylabel("Time for 10 iterations [s]",fontsize=35)
xticks(fontsize=22)
yticks(fontsize=22)
#x=scipy.ones(100)*1010
#y=scipy.linspace(0,20,100)
#plot(x,y,linewidth=3) #Start moving window
title("Domain size = 6912 X 128 cells")
show()


