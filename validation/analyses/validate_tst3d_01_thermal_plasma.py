import os, re, numpy as np, math 
import happi

S = happi.Open(["./restart*"], verbose=False)



# 3D SCREEN DIAGS
precision = [0.02, 0.06, 0.01, 0.06, 0.03, 0.1, 0.02, 0.1]
for i,d in enumerate(S.namelist.DiagScreen):
	last_data = S.Screen(i, timesteps=160).getData()[-1]
	Validate("Screen "+d.shape+" diag with "+d.direction+" direction", last_data, precision[i])
