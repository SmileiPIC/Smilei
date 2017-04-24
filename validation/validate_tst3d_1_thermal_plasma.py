import os, re, numpy as np, math 
from Smilei import *

S = Smilei(".", verbose=False)

# 3D SCREEN DIAGS
precision = [0.02, 0.05, 0.01, 0.05, 0.03, 0.1, 0.02, 0.1]
for i,d in enumerate(S.namelist.DiagScreen):
	last_data = S.Screen(i, timesteps=160).getData()[-1]
	Validate("Screen "+d.shape+" diag with "+d.direction+" direction", last_data, precision[i])
