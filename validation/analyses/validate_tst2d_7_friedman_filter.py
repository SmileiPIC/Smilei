import os, re, numpy as np, math 
import happi

S = happi.Open(".", verbose=False)

# MINIMAL VALIDATION
Validate("Number of timesteps", len( S.Scalar.Ubal_norm().getData()), 1)
Validate("Scalar Ubal_norm", S.Scalar.Ubal_norm().getData()[500], 0.01)

