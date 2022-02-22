import os, re, numpy as np, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

Utot = S.Scalar.Utot().getData()
Ukin = S.Scalar.Ukin().getData()
Validate("Scalar Utot", Utot / Utot[0], 1e-6)
Validate("Scalar Ukin", Ukin / Ukin[0], 1e-6)
