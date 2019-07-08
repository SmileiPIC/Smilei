import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Ey FIELD in polarization direction
Ey = S.Probe(0, "Ey", timesteps=2000.).getData()[0]
Validate("Ey field at iteration 2000", Ey, 0.01)

# COMPARE THE Jy FIELD in polarization direction
Jy = S.Probe(0, "Jy", timesteps=2000.).getData()[0]
Validate("Jy field at iteration 2000", Jy, 0.0005)

