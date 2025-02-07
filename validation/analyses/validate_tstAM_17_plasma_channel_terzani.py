import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Ey FIELD in polarization plane
Ey      = S.Probe(1, "Ey",timestep_indices=-1).getData()[0][::4,::4]
Validate("Ey field at last iteration", Ey, 0.01)

# COMPARE THE Ex FIELD in polarization plane
Ex      = S.Probe(1, "Ex",timestep_indices=-1).getData()[0][::4,::4]
Validate("Ex field at last iteration", Ex, 0.01)

# COMPARE THE BzBTIS3 FIELD in polarization plane
BzBTIS3 = S.Probe(1, "BzBTIS3",timestep_indices=-1).getData()[0][::4,::4]
Validate("BzBTIS3 field at last iteration", BzBTIS3, 0.01)

