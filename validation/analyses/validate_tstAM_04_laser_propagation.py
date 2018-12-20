import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Ey FIELD in polarization direction
Ey = S.Field(0, "Et", theta=math.pi/2., timesteps=1000., modes=1, subset={"x":[0,S.namelist.Main.grid_length[0]/4.],"y":[0,S.namelist.Main.grid_length[1]/2.]}).getData()[0]
Validate("Ey field at iteration 700", Ey, 0.01)

# COMPARE THE Jy FIELD in polarization direction
Jy = S.Field(0, "Jt", theta=math.pi/2., timesteps=1000., modes=1, subset={"x":[0,S.namelist.Main.grid_length[0]/4.],"y":[0,S.namelist.Main.grid_length[1]/2.]}).getData()[0]
Validate("Jy field at iteration 700", Jy, 0.0005)

# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5") as f:
	dt = f["data/0000000000"].attrs["dt"]
	dx = f["data/0000000000/Er_mode_0"].attrs["gridSpacing"]
	patchSize = f["data/0000000000"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)
