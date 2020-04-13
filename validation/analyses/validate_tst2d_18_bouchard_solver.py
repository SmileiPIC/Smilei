import numpy as np
import happi

s = happi.Open(["./restart*"], verbose=False)

# COMPARE THE bz VIA PROBE
# bz = np.array(s.Probe.Probe0.Bz(timesteps=2000).getData()[0])
# Validate("Bz field at iteration 2000", bz, 0.01)

# TEST THAT ubal_norm STAYS OK
ubal_norm = s.Scalar.Ubal_norm(timesteps=(0,3000)).getData()
Validate("Ubal_norm comparison", ubal_norm, 0.01)
