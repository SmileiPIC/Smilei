#!/opt/local/bin/python

import tables
import numpy as np
import matplotlib.pyplot as plt

f=tables.openFile("Probe0d.h5")


for n in f.root.Data:
    np.savetxt(n.name+".txt", n.read())

        
f.close()
