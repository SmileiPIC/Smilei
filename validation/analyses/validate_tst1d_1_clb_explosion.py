import os, re, numpy as np
import happi

S = happi.Open(["./restart*"], verbose=False)



Ntot_charges = S.Scalar.Ntot_charges(timesteps=0).getData()[0]
Validate("Initial number of particles", Ntot_charges )

momentum_distribution = S.ParticleBinning.Diag0(timesteps=0).getData()[0]
Validate("Initial momentum distribution", momentum_distribution, 0.001 )

max_ubal = np.max( np.abs(S.Scalar.Ubal().getData()) )
Validate("Max Ubal is below 2%", max_ubal<0.02 )

Validate("List of fields in Field0", S.Field.Field0().getFields() )
