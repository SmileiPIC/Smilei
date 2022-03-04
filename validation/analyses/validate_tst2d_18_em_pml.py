import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



# COMPARE THE Ez FIELD
Ez = S.Field.Field0.Ez(timesteps=300).getData()[0]
Validate("Ez field at iteration 300", Ez, 0.001)

# TEST THAT Ubal_norm STAYS OK
uelm = S.Scalar.Uelm().getData(timestep=300)
Validate("Uelm is below 0.2% after absorption", uelm[0]<0.002 )

