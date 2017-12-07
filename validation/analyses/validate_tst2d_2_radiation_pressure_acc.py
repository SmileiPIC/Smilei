import os, re, numpy as np, math
from scipy.ndimage import gaussian_filter as gfilt
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE INITIAL POSITIONS OF ELECTRONS
particles = S.TrackParticles.eon(axes=["x","y"], select="any(t==0, (x<7)*(y<0.5))").getData()
Validate("Regular electron x at iteration 0", particles["x"][0])
Validate("Regular electron y at iteration 0", particles["y"][0])

# COMPARE THE ELECTRON FRONT POSITION
def getFrontPosition(timestep):
	ne = S.Field.Field0.Rho_eon(timesteps=timestep).getData()[0] # electron density
	ne = gfilt(ne, 30) # smoothing
	ne = np.abs(ne[:,490:510].mean(axis=1)) # profile
	return np.argmax(ne)
front_position = [getFrontPosition(t) for t in [500, 1000, 1500]]
Validate("Electron front vs time", front_position, 20)

# 2D SCREEN DIAGS
for i,d in enumerate(S.namelist.DiagScreen):
	last_data = S.Screen(i, timesteps=1400).getData()[-1]
	if d.direction in ["backward", "canceling"]:
		precision = 8
	else:
		precision = 20
	if d.shape == "sphere": precision *= 0.02
	Validate("Screen "+d.shape+" diag with "+d.direction+" direction", last_data, precision)
