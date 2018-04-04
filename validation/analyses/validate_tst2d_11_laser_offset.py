import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE EM ENERGY AT LASER FOCUS
U = S.Field(0,"Bz**2+Ey**2", subset={"x":40,"y":[141,173]}, timesteps=675).getData()[0]
Validate("EM energy profile at laser focus", U, 0.1)
