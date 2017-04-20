import os, re, numpy as np
from Smilei import *

S = Smilei(".", verbose=False)

some_particles_x = S.TrackParticles.ion(axes=["x"],select="any(t==0, x<0.02)").getData()["x"][0]
Validate("Regularly spaced particles", some_particles_x )

momentum_distribution = S.ParticleDiagnostic.Diag0(timesteps=0, slice={"x":"all"}).getData()[0]
Validate("Electron momentum distribution", momentum_distribution, 1e-7 )

itimes = S.ParticleDiagnostic.Diag0().getAvailableTimesteps()
Validate("Timesteps in particle diagnostic", itimes )

phase_space = S.ParticleDiagnostic.Diag0(timesteps=2000).getData()[0]
Validate("Electron phase-space", phase_space, 1.e-3)