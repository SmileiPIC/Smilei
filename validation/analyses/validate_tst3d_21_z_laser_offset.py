import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE EM ENERGY AT LASER FOCUS
U = S.Field(0,"Ex**2+Ey**2+Ez**2", subset={"z":[5,20,2],"x":[32,45,2],"y":31.}, timesteps=150).getData()[0]
Validate("EM energy profile at laser focus", U, 0.1)
