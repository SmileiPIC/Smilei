import os, re, numpy as np
from Smilei import *

S = Smilei(".", verbose=False)

ion_profiles = S.Field.Field0.Rho_eon().getData()
ion_front = []
for p in ion_profiles: ion_front += [np.flatnonzero(p<-10)[0]]
ion_front = np.array(ion_front)
Validate("Ion front position vs time", ion_front, 100.)

