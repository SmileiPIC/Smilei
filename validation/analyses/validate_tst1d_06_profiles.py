import os, re, numpy as np
from scipy.signal import butter, filtfilt
import happi

S = happi.Open(["./restart*"], verbose=False)



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
	Validate("Temporal profile "+name, thisJz[1:], 0.01)

# Verify that "time_fields_frozen" and external fields works
#for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
	F = S.Field.Field0(field, timesteps=190).getData()[0]
	Validate(field+" field at late time", F, 0.01)

# Verify that "time_frozen" works
Late_density = S.Field.Field0.Rho_trapezoidal(timesteps=190).getData()[0]
Validate("Density of a species at a late time", Late_density, 0.01)

# Maxwell-Juttner initialization
i = -1
drift = ["no","px","py","pz"]
for eon in S.namelist.mj_species:
	i+=1
	p = S.ParticleBinning(i).get()
	p_distr = p["data"][0]
	p = p[drift[i]] if i>0 else p["px"]
#	# theory
#	g0 = S.namelist.g0 if i>0 else 1.
#	Te = S.namelist.Te
#	fth = (g0 * np.sqrt(1.+p**2)+Te) * np.exp( -g0/Te* (np.sqrt(1.+p**2) - np.sqrt(1.-g0**-2)*p) )
#	itg = (p[1]-p[0])*np.sum(fth)
#	fth = fth/itg
#	plt.plot(px, p_distr, 'o', markeredgecolor=c[i], markerfacecolor="none")
#	plt.plot(px, fth, '-', color=c[i])
	b, a = butter(8, 0.15, btype='low', analog=False)
	p_filt = filtfilt(b, a, p_distr)
	Validate("Maxwell-Juttner Momentum distribution ("+drift[i]+" drift)", p_filt, p_filt.max()*1e-2)

