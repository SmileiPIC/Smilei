import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE EM ENERGY AT LASER FOCUS
U = S.Field(0,"Ex**2+Ey**2", subset={"x":[9,25,2],"y":[176,190,2]}, timesteps=500).getData()[0]
Validate("EM energy profile at laser focus", U, 0.1)
