import os, re, numpy as np
from Smilei import *

S = Smilei(".", verbose=False)

some_particles_x = S.TrackParticles.ion(axes=["x"]).getData()["x"][0]
Validate("Regularly spaced particles", some_particles_x, 1e-7)

itimes = S.ParticleBinning.Diag0().getAvailableTimesteps()
Validate("Timesteps in particle diagnostic", itimes )

px = S.ParticleBinning.Diag0(slice={"x":"all"}, timesteps=0   ).getData()[0]
Validate("Initial electron momentum distribution", px, 0.1 )

px = S.ParticleBinning.Diag0(slice={"x":"all"}, timesteps=2000).getData()[0]
Validate("Final electron momentum distribution", px, 0.1)