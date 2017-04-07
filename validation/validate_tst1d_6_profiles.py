import os, re, numpy as np
from Smilei import *

S = Smilei(".", verbose=False)

for name, profile in S.namelist.profiles.items():
	A = S.Field.Field0("Rho_"+name)
	data = A.get()
	values = data["data"][0]
	# x = data["x"]
	# v = [profile(xx) for xx in x]
	Validate("Profile "+name, values, 0.01 )

