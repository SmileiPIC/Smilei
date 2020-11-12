import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Ey FIELD in polarization plane
Ey = S.Field(0, "Er", theta=0, timesteps=1700.).getData()[0][::4,::4]
Validate("Ey field at iteration 1700", Ey, 0.01)

# COMPARE THE Ey FIELD in polarization plane + SUBGRID
Ey = S.Field(0, "Er", theta=0, timesteps=1700.).getData()[0][::4,::4]
Validate("Ey field subgrid at iteration 1700", Ey, 0.01)

# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5", "r") as f:
	dt = f["data/0000000000"].attrs["dt"]
	dx = f["data/0000000000/Er_mode_0"].attrs["gridSpacing"]
	patchSize = f["data/0000000000"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)
