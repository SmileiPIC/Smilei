# ____________________________________________________________________________
#
# This script validates the adaptive vectorization in mode 1
#
# _____________________________________________________________________________

import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# Scalars
ukin = S.Scalar("Ukin").getData()
utot = S.Scalar("Utot").getData()

Validate("Total kinetic energy evolution: ", ukin / ukin[0], 1e-3 )
Validate("Total energy evolution: ", utot / utot[0], 1e-3 )
