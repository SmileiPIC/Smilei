import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(timesteps=1500, subset={"x":[0,10000,8], "y":[0,10000,8]}).getData()[0]
Validate("Ey field at iteration 1500", Ey, 0.01)

# SUBSET OF FIELD DIAG
Ey  = S.Field.Field0.Ey(timesteps=1500).getData()[0]
Ey1 = S.Field.Field1.Ey(timesteps=1500).getData()[0]
subset = S.namelist.DiagFields[1].subset
Validate("Field subset works", (Ey1==Ey[subset]).all())

# 2-D probe in 2D
Ey = S.Probe.Probe0.Ey(timesteps=1500).getData()[0]
Validate("Ey probe at iteration 1500", Ey, 0.01)

# 0-D probe in 2D
Ey = S.Probe.Probe1.Ey().getData()
Validate("0-D probe Ey vs time", Ey, 0.01)

# TEST THAT Ubal_norm STAYS OK
max_ubal_norm = np.max( np.abs(S.Scalar.Ubal_norm().getData()) )
Validate("Max Ubal_norm is below 10%", max_ubal_norm<0.1 )

# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5") as f:
	dt = f["data/0000000000"].attrs["dt"]
	dx = f["data/0000000000/Ex"].attrs["gridSpacing"]
	patchSize = f["data/0000000000"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)
