import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE PLASMA DENSITY @ t=0
N0 = S.ParticleBinning(0,sum={"z":"all"}, timesteps=0).getData()[0]
Validate("Plasma density at iteration 0", N0, 0.001)

# COMPARE THE PLASMA DENSITY @ t=1500
N1 = S.ParticleBinning(0,sum={"z":"all"}, timesteps=1500).getData()[0]
Validate("Plasma density at iteration 1500", N1, 0.001)
