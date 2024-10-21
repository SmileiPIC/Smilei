import os, re, numpy as np, math
from scipy.ndimage import gaussian_filter as gfilt
import happi

#S = happi.Open(["./restart*"], verbose=False)

S = happi.Open(verbose=False)

# COMPARE THE TOTAL NUMBER OF CREATED IONS
nion = S.Scalar.Ntot_ion(timesteps=0).getData()[0]
Validate("Number of ions at iteration 0", nion)

# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(timesteps=800).getData()[0]
Ey = gfilt(Ey, 6) # smoothing
Ey = Ey[50:200:4, 300::4] # zoom on the reflected laser
Validate("Ey field at iteration 800", Ey, 0.02)
