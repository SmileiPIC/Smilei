import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Env_E_abs FIELD
Env_E_abs = S.Probe(0, "Env_E_abs", timesteps=1000.).getData()[0]
Validate("Env_E_abs field at iteration 1000", Env_E_abs, 0.01)

# COMPARE THE Ex FIELD
Ex = S.Probe(0, "Ex", timesteps=1000.).getData()[0]
Validate("Ex field at iteration 1000", Ex, 0.01)

# COMPARE THE Rho FIELD
Rho = S.Probe(0, "Rho", timesteps=1000.).getData()[0]
Validate("Rho field at iteration 1000", Rho, 0.01)

# COMPARE THE Rho_electronfromion FIELD
Rho_el = S.Probe(0, "Rho_electronfromion", timesteps=1000.).getData()[0]
Validate("Rho_electronfromion field at iteration 1000", Rho_el, 0.01)







