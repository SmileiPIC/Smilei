import os, re, numpy as np, math 
from Smilei import *

S = Smilei(".", verbose=False)

ncel = [int(i) for i in S._ncels]

for field in ["Ex", "Ey", "Ez", "Bx_m", "By_m", "Bz_m", "Jx", "Jy", "Jz", "Rho"]:
	if field in ["Ex", "Rho"]:
		precision = 0.2
	else:
		precision = 0.02
	Validate("Minimum of scalar "+field, S.Scalar(field+"Min").getData()[-1], precision)
	Validate("Maximum of scalar "+field, S.Scalar(field+"Max").getData()[-1], precision)
	
	MinLoc = np.unravel_index(int(S.Scalar(field+"MinCell").getData()[-1]), ncel)
	MaxLoc = np.unravel_index(int(S.Scalar(field+"MaxCell").getData()[-1]), ncel)
	Validate("Location of minimum of scalar "+field, MinLoc, 6)
	Validate("Location of maximum of scalar "+field, MaxLoc, 6)

