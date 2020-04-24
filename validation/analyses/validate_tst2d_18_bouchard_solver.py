import os,re
import numpy as np
import happi

s0 = happi.Open(["./restart*"], verbose=False)

# COMPARE THE bz VIA PROBE
bz = np.array(s0.Field.Field0.Bz(timesteps=1500).getData()[0]).T
Validate("Bz field at iteration 1500", bz, 0.01)

# TEST THAT ubal_norm STAYS OK
ubal_norm = s0.Scalar.Ubal_norm(timesteps=(0,1500)).getData()
Validate("Ubal_norm comparison", ubal_norm, 0.01)
