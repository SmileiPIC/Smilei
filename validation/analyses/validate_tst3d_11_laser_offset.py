import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE EM ENERGY AT LASER FOCUS
U = S.Field(0,"Bz**2+Ey**2", subset={"x":10,"y":[21,40,3],"z":[21,40,3]}, timesteps=80).getData()[0]
Validate("EM energy profile at laser focus", U, 0.1)
