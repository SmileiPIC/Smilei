import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE Ey FIELD
Ex = S.Field.Field0.Ex(average={"z":S.namelist.Main.grid_length[2]*0.3}, subset={"x":[20,50,2], "y":[0,10000,2]}, timesteps=240).getData()[0]
Validate("Ex field at iteration 240", Ex, 0.01)
