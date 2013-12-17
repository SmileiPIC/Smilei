#!/opt/local/bin/python

import tables
import numpy as np
import matplotlib.pyplot as plt

num_probe = 3
fields = ['Ex','Ey','Ez','Bx','By','Bz']


for p in range(0,num_probe):
    name="probe_%d.h5" %(p)
    f=tables.openFile(name)
    d=[]
    for node in f.root:
        d.append(node.read())
    d1=np.array(d)
    for l in range(0,len(fields)):
        plt.figure()
        plt.plot(d1[:,0],d1[:,l+1])
        plt.grid(True)
        plt.xlabel('Time')
        plt.ylabel('%s'%(fields[l]))
        plt.title('Probe0D_%d'%(p))
        plt.savefig('Probe0D_%d'%(p)+'_%s'%(fields[l]))
    	
