import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# COMPARE THE Ey FIELD in polarization direction
Ey = S.Probe(0, "Ey", timesteps=2000.).getData()[0]
Validate("Ey field at iteration 2000", Ey, 0.01)

# COMPARE THE Jy FIELD in polarization direction
Jy = S.Probe(0, "Jy", timesteps=2000.).getData()[0]
Validate("Jy field at iteration 2000", Jy, 0.0005)

## Performances non regression
#timer_particle = S.Performances(raw="timer_particles").getData()[-1].mean()
#Validate("Mean time spent in particles", timer_particle, 2.)

#timer_mw = S.Performances(raw="timer_movWindow").getData()[-1].mean()
#Validate("Mean time spent in moving windows", timer_mw, 0.5)
#
#timer_syncdens = S.Performances(raw="timer_syncDens").getData()[-1].mean()
#Validate("Mean time spent in sync densities", timer_syncdens, 2.)

