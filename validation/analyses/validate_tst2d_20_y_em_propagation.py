import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(timesteps=1500, subset={"x":[0,10000,8], "y":[0,10000,8]}).getData()[0]
Validate("Ey field at iteration 1500", Ey, 0.01)
