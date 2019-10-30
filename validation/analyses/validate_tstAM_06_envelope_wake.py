import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE Env_A_abs FIELD
Env_A_abs = S.Field(0, "Env_A_abs", theta=0, timesteps=1700.).getData()[0]
Validate("Env_A_abs field at iteration 1700", Env_A_abs, 0.01)

# COMPARE THE Ex FIELD
El = S.Field(1, "El", theta=0, timesteps=1700.).getData()[0]
Validate("El field at iteration 1700", El, 0.01)

# COMPARE THE Env_Chi FIELD
Env_Chi = S.Field(0, "Env_Chi", theta=0, timesteps=1700.).getData()[0]
Validate("Env_Chi field at iteration 1700", Env_Chi, 0.01)


# 1-D PROBE IN AM
Env_A_abs = S.Probe.Probe0.Env_A_abs(timesteps=1700).getData()[0]
Validate("1-D probe Env_A_abs at iteration 1700", Env_A_abs, 0.01)

Ex = S.Probe.Probe0.Ex(timesteps=1700).getData()[0]
Validate("1-D probe Ex at iteration 1700", Ex, 0.01)

Env_Chi = S.Probe.Probe0.Env_Chi(timesteps=1700).getData()[0]
Validate("1-D probe Env_Chi at iteration 1700", Env_Chi, 0.01)


