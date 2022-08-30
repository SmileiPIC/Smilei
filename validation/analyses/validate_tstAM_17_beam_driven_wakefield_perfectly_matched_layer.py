import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)
um = S.namelist.um
dt = S.namelist.dt

# COMPARE THE Ex FIELD after 1 um 
Ex = S.Probe(0, "Ex", timesteps=int(1000*um/dt)).getData()[0]
Validate("Ex field after 1 mm", Ex, 0.01)

