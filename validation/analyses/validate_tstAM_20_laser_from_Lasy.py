import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Er FIELD in polarization plane
Er = S.Field(0, "Er", theta=0, timesteps=1700.).getData()[0][::4,::4]
Validate("Er field subgrid at iteration 1700", Er, 0.01)

# COMPARE THE Et FIELD in polarization plane + SUBGRID
Et = S.Field(0, "Et", theta=0, timesteps=1700.).getData()[0][::4,::4]
Validate("Et field subgrid at iteration 1700", Et, 0.01)

# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5", "r") as f:
	dt = f["data/0"].attrs["dt"]
	dx = f["data/0/Er_mode_0"].attrs["gridSpacing"]
	patchSize = f["data/0"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)
