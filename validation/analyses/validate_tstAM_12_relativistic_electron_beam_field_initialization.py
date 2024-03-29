import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)


# E field on grid
Ex = S.Probe(0, "Ex",timesteps=0.).getData()
Validate("Field Ex" , Ex, 0.01)

Ey = S.Probe(0, "Ey",timesteps=0.).getData()
Validate("Field Ey" , Ey, 0.01)

Ez = S.Probe(0, "Ez",timesteps=0.).getData()
Validate("Field Ez" , Ez, 0.01)


# B field on grid
Bx = S.Probe(0, "Bx",timesteps=0.).getData()
Validate("Field Bx" , Bx, 0.01)

By = S.Probe(0, "By",timesteps=0.).getData()
Validate("Field By" , By, 0.01)

Bz = S.Probe(0, "Bz",timesteps=0.).getData()
Validate("Field Bz" , Bz, 0.01)


# The relativistic field depends mainly on the energy 
# of the relativistic bunch, ~px, so this is validated
# through the validation on the EM fields

# sum of py
Py = S.ParticleBinning(0, timesteps=0).getData()[0]
Validate("Weighted sum of py at iteration 0", Py, 0.001)

# sum of pz
Pz = S.ParticleBinning(1, timesteps=0).getData()[0]
Validate("Weighted sum of pz at iteration 0", Pz, 0.001)
