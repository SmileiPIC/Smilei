import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Er FIELD in polarization direction
Er = S.Field(0, "Er", timesteps=200., theta=0.).getData()[0]
Validate("Er field at iteration 200", Er, 0.01)

