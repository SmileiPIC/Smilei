import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Envelope FIELD, absolute value
Env_A_abs = S.Field(0, "Env_A_abs_mode_0", theta=0, timesteps=300.).getData()[0]
Validate("Env_A_abs field at iteration 300", Env_A_abs, 0.01)

# 1-D PROBE IN AM
Env_A_abs = S.Probe.Probe0.Env_A_abs(timesteps=300).getData()[0]
Validate("1-D probe Env_A_abs at iteration 300", Env_A_abs, 0.01)


