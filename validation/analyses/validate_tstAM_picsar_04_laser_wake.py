import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Ey FIELD in polarization direction
rho = S.Field(0, "Rho",theta=0., timesteps=0).getData()[0]
rho = np.abs(rho).sum
Validate("Rho field at iteration 0", rho == 0.)

