import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# 0-D PROBE IN 3D
Ey = S.Probe.Probe0.Ey().getData()
Validate("0-D probe Ey vs time", Ey, 0.01)
# 1-D PROBE IN 3D
Ey = S.Probe.Probe1.Ey(timesteps=50).getData()[0]
Validate("1-D probe Ey at last iteration", Ey, 0.01)
# 2-D PROBE IN 3D
Ey = S.Probe.Probe2.Ey(timesteps=50).getData()[0]
Validate("2-D probe Ey at last iteration", Ey, 0.01)
Ex = S.Probe.Probe2.Ey(timesteps=50).getData()[0]
Validate("2-D probe Ey at last iteration", Ex, 0.01)
# 3-D PROBE IN 3D
#Ey = S.Probe.Probe3.Ey(timesteps=100).getData()[0]
#Validate("3-D probe Ey at last iteration", Ey, 0.01)


## TEST THE GRID PARAMETERS
#with h5py.File("./restart000/Fields0.h5", "r") as f:
#	dt = f["data/0000000000"].attrs["dt"]
#	dx = f["data/0000000000/Ex"].attrs["gridSpacing"]
#	patchSize = f["data/0000000000"].attrs["patchSize"]
#Validate("Value of the timestep" , dt, 1e-6)
#Validate("Value of the grid step", dx, 1e-6)
#Validate("Patch size", patchSize)
