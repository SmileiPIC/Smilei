import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(average={"z":S.namelist.Main.grid_length[2]*0.3}, timesteps=90).getData()[0]
Validate("Ey field at iteration 90", Ey, 0.01)

# SUBGRID OF FIELD DIAG
Ey  = S.Field.Field0.Ey(timesteps=90).getData()[0]
Ey1 = S.Field.Field1.Ey(timesteps=90).getData()[0]
subgrid = S.namelist.DiagFields[1].subgrid
Validate("Field subgrid works", (Ey1==Ey[subgrid]).all())

# 0-D PROBE IN 3D
Ey = S.Probe.Probe0.Ey().getData()
Validate("0-D probe Ey vs time", Ey, 0.01)
# 1-D PROBE IN 3D
Ey = S.Probe.Probe1.Ey(timesteps=90).getData()[0]
Validate("1-D probe Ey at last iteration", Ey, 0.01)
# 2-D PROBE IN 3D
Ey = S.Probe.Probe2.Ey(timesteps=90).getData()[0]
Validate("2-D probe Ey at last iteration", Ey, 0.01)
# 3-D PROBE IN 3D
Ey = S.Probe.Probe3.Ey(timesteps=90).getData()[0]
Validate("3-D probe Ey at last iteration", Ey, 0.01)


# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5") as f:
	dt = f["data/0000000000"].attrs["dt"]
	dx = f["data/0000000000/Ex"].attrs["gridSpacing"]
	patchSize = f["data/0000000000"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)
