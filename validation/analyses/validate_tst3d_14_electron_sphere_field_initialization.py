import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)


# E field on grid
Ex = S.Probe(1, "Ex",timesteps=0.).getData()[0]
Validate("Probe Ex" , Ex, 0.01)

Ey = S.Probe(1, "Ey",timesteps=0.).getData()[0]
Validate("Probe Ey" , Ey, 0.01)

Ez = S.Probe(1, "Ez",timesteps=0.).getData()[0]
Validate("Probe Ez" , Ez, 0.01)


