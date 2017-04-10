import os, re, numpy as np
from Smilei import *

S = Smilei(".", verbose=False)

c = ["b","r","g","c","m","y"]

# Spatial profiles
for i, (name, profile) in enumerate(S.namelist.profiles.items()):
	A = S.Field.Field0("Rho_"+name)
	data = A.get()
	values = data["data"][0]
#	x = data["x"]
#	v = [profile(xx) for xx in x]
#	plt.plot(x, values, color=c[i])
#	plt.plot(x, v, "o", markeredgecolor=c[i], markerfacecolor="none")
	Validate("Profile "+name, values, 0.01 )
	

# Temporal profiles
Jz = S.Field.Field0("Jz").get()
x = Jz["x"]
#t = Jz["times"]
Jz = np.array(Jz["data"])
w = S.namelist.antenna_width
for i, (name, tprofile) in enumerate(S.namelist.tprofiles):
	x0 = i*w + 0.1*w
	x1 = x0 + 0.8*w
	thisJz = Jz[:, (x>x0)*(x<x1)].mean(axis=1)
#	v = [tprofile(tt) for tt in t*S.namelist.Main.timestep]
#	plt.plot(t, thisJz, color=c[i])
#	plt.plot(t, v, "o", markeredgecolor=c[i], markerfacecolor="none")
	Validate("Temporal profile "+name, thisJz, 0.01)
