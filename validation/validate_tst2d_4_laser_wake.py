import os, re, numpy as np, math 
from Smilei import *

S = Smilei(".", verbose=False)

# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(timesteps=1600).getData()[0][::10,:]
Validate("Ey field at iteration 1600", Ey, 0.1)

# CHECK THE LOAD BALANCING
with open("patch_load.txt") as f:
	txt = f.read()
patch_count0 = re.findall(r"patch_count\[0\] = (\d+)",txt)
patch_count1 = re.findall(r"patch_count\[1\] = (\d+)",txt)
initial_balance = [int(patch_count0[0] ), int(patch_count1[0] )]
final_balance   = [int(patch_count0[-1]), int(patch_count1[-1])]
Validate("Initial load balance", initial_balance, 1)
Validate("Final load balance", final_balance, 1)

