import os, re, numpy as np, math 
from Smilei import *

S = Smilei(".", verbose=False)

ncel = [int(i) for i in S._ncels]

# SCALARS RELATED TO MIN/MAX OF FIELDS
for field in ["Ex", "Ey", "Ez", "Bx_m", "By_m", "Bz_m", "Jx", "Jy", "Jz", "Rho"]:
	if field in ["Ex", "Rho"]:
		precision = 0.25
	else:
		precision = 0.025
	Validate("Minimum of scalar "+field, S.Scalar(field+"Min").getData()[-1], precision)
	Validate("Maximum of scalar "+field, S.Scalar(field+"Max").getData()[-1], precision)
	
	MinLoc = np.unravel_index(int(S.Scalar(field+"MinCell").getData()[-1]), ncel)
	MaxLoc = np.unravel_index(int(S.Scalar(field+"MaxCell").getData()[-1]), ncel)
	Validate("Location of minimum of scalar "+field, MinLoc, 10)
	Validate("Location of maximum of scalar "+field, MaxLoc, 10)

# FIELD DIAGNOSTICS
fields     = ["Ex","Ey" ,"Ez" ,"Bx" ,"By" ,"Bz" ,"Bx_m","By_m","Bz_m","Jx" ,"Jy" ,"Jz" ,"Rho","Jx_test0","Jy_test0","Jz_test0","Rho_test0"]
precisions = [ 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 , 0.01 , 0.01 , 0.01, 0.01, 0.01, 0.1 , 3e-5     , 2e-5     , 0.001    , 0.1       ]
for field, precision in zip(fields, precisions) :
	Validate("Field "+field, S.Field.Field0(field, average={"z":"all"}, timesteps=40, subset={"x":[0,10,4], "y":[0,10,4]}).getData()[-1], precision)

# PROBE DIAGNOSTICS
fields     = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx", "Jy", "Jz", "Rho"]
precisions = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 , 0.3  ]
for field, precision in zip(fields, precisions):
	Validate("Probe field "+field, S.Probe(0, field).getData()[-1], precision)
