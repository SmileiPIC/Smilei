import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE all the electromagnetic FIELD
Ex = S.Field.Field0.Ex(subset={"z":S.namelist.Main.grid_length[2]/2,"x":[100,160,2],"y":[50,150,2]}, timesteps=0.).getData()[0]
Validate("Ex field at iteration 0", Ex, 0.01)

Ey = S.Field.Field0.Ey(subset={"z":S.namelist.Main.grid_length[2]/2,"x":[100,160,2],"y":[50,150,2]}, timesteps=0.).getData()[0]
Validate("Ey field at iteration 0", Ey, 0.01)

Ez = S.Field.Field0.Ez(subset={"z":S.namelist.Main.grid_length[2]/2,"x":[100,160,2],"y":[50,150,2]}, timesteps=0.).getData()[0]
Validate("Ez field at iteration 0", Ez, 0.01)

Bx = S.Field.Field0.Bx(subset={"z":S.namelist.Main.grid_length[2]/2,"x":[100,160,2],"y":[50,150,2]}, timesteps=0.).getData()[0]
Validate("Bx field at iteration 0", Bx, 0.01)

By = S.Field.Field0.By(subset={"z":S.namelist.Main.grid_length[2]/2,"x":[100,160,2],"y":[50,150,2]}, timesteps=0.).getData()[0]
Validate("By field at iteration 0", By, 0.01)

Bz = S.Field.Field0.Bz(subset={"z":S.namelist.Main.grid_length[2]/2,"x":[100,160,2],"y":[50,150,2]}, timesteps=0.).getData()[0]
Validate("Bz field at iteration 0", Bz, 0.01)

# 1-D PROBE IN 3D
Ex = S.Probe.Probe0.Ex(timesteps=0.).getData()[0]
Validate("1-D probe for initialized Ex at first iteration", Ex, 0.01)

# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5") as f:
	dt = f["data/0000000000"].attrs["dt"]
	dx = f["data/0000000000/Ex"].attrs["gridSpacing"]
	patchSize = f["data/0000000000"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)
