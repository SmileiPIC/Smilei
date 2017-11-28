import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)



# COMPARE THE Ez FIELD
Ez = S.Field.Field0.Ez(timesteps=55).getData()[0]
Ez = Ez[70:160:2,50:150:2]
Validate("Ez field at iteration 55", Ez, 0.02)

# VERIFY THE PRESENCE OF ELECTRONS BEFORE AND AFTER THE WALL
ne = S.Field.Field0.Rho_electron(timesteps=300).getData()[0]
before_wall = ne[:152,:].sum() / ne.sum()
after_wall  = ne[152:,:].sum() / ne.sum()
Validate("Electrons before and after wall", [before_wall, after_wall], 1e-7)