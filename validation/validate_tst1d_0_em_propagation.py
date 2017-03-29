import os, re, numpy as np
from Smilei import *

S = Smilei(".", verbose=False)

Validate("List of fields in Field0", S.Field.Field0().getFields() )

timesteps = list(S.Field.Field0().getAvailableTimesteps())
Validate("List of timesteps in Field0", timesteps )

Ez = S.Field.Field0("Ez", timesteps=timesteps[-1]).getData()[0]
Validate("Last Ez profile in Field0", Ez, 1e-7 )

Validate("List of fields in Probe1", S.Probe(1).getFields() )

timesteps = S.Probe(1,"Ez").getAvailableTimesteps()
Validate("List of timesteps in Probe1", timesteps )

Ez = S.Probe(1,"Ez", timesteps=timesteps[-1]).getData()[0]
Validate("Last Ez profile in Probe1", Ez, 1e-7 )

max_ubal = np.max( np.abs(S.Scalar.Ubal().getData()) )
Validate("Max Ubal is below 2%", max_ubal<0.02 )

