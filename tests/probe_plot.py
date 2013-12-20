# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/tommaso/.spyder2-27/.temp.py
"""

import tables 
from pylab import *

clf()
f=tables.openFile("/Users/tommaso/local/src/smilei/tests/Probes0D.h5")

d=[]
pos=[]
for n in f.root:
    d.append(n.read())
    pos.append(n.name+" "+str(n.attrs.Position[0]))
    
d1=np.array(d)

print np.amin(d1[:,:,1:]), np.amax(d1[:,:,1:])

cols=[3,5]

for col in cols:
    for probe in range(0,d1.shape[0]):
        p,= plot(d1[0,:,0],probe+d1[probe,:,col], label=pos[probe]+" "+str(col))

#    for i in range(1,d1.shape[2]):
#        d2=d1[1,:,i]
#        diff=np.amax(d2)-np.amin(d2)
#        print probe,diff

legend(loc = 'upper left')

ion()
show()

f.close()