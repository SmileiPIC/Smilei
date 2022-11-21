import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Ex FIELD
Ex = S.Probe(0, "Ex", timesteps=800.).getData()[0]
Validate("Ex field at iteration 800", Ex, 0.01)

# COMPARE THE Rho FIELD
Rho = S.Probe(0, "Rho", timesteps=800.).getData()[0]
Validate("Rho field at iteration 800", Rho, 0.01)

# COMPARE THE Rho_electronfromion FIELD
Rho_el = S.Probe(0, "Rho_electronfromion", timesteps=800.).getData()[0]
Validate("Rho_electronfromion field at iteration 800", Rho_el, 0.01)







