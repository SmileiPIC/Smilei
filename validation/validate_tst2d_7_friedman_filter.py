import os, re, numpy as np, math 
from Smilei import *

S = Smilei(".", verbose=False)

# MINIMAL VALIDATION
Validate("Number of timesteps", len( S.Scalar.Ubal_norm().getData()), 1)
Validate("Scalar Ubal_norm", S.Scalar.Ubal_norm().getData()[500], 0.01)

