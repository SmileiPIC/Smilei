import os, re, numpy as np, math
from Smilei import *

S = Smilei(".", verbose=False)

# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(timesteps=1500, stride=8).getData()[0]
Validate("Ey field at iteration 1500", Ey, 0.01)

# COMPARE THE Ey PROBE
Ey = S.Probe.Probe0.Ey(timesteps=1500).getData()[0]
Validate("Ey probe at iteration 1500", Ey, 0.01)

# TEST THAT Ubal_norm STAYS OK
max_ubal_norm = np.max( np.abs(S.Scalar.Ubal_norm().getData()) )
Validate("Max Ubal_norm is below 10%", max_ubal_norm<0.1 )

