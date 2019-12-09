import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Ey FIELD in polarization direction
Ey = S.Probe(0, "Ey", timesteps=200.).getData()[0]
Validate("Ey field at iteration 200", Ey, 0.01)

