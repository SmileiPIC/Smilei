import happi
import numpy as np
import os, re

S = happi.Open(["./restart*"], verbose=False)

### Compare the carbon charge density & corresponding electron charge density at last timestep
Rho_carbon = S.Field.Field0.Rho_carbon(timesteps=1280).getData()[0]
Validate("Carbon charge density at last timestep", Rho_carbon, 0.075)
Rho_electronC = S.Field.Field0.Rho_electronC(timesteps=1280).getData()[0]
Validate("Carbon-electron charge density at last timestep", Rho_electronC, 0.075)

### Compare the time-evolution of the Carbon charge state from ParticleBinning
n  = np.array( S.ParticleBinning(0).getData() )
n /= n[0,0]
n0 = n[:,0]
n1 = n[:,1]
n2 = n[:,2]
Validate("Time evolution of Carbon charge state Z=0", n0, 0.02)
Validate("Time evolution of Carbon charge state Z=1", n1, 0.02)
Validate("Time evolution of Carbon charge state Z=2", n2, 0.02)

### Compare the hydrogen charge density & corresponding electron charge density at last timestep
Rho_hydrogen = S.Field.Field0.Rho_hydrogen(timesteps=1280).getData()[0]
Validate("Hydrogen charge density at last timestep", Rho_hydrogen, 0.075)
Rho_electronH = S.Field.Field0.Rho_electronH(timesteps=1280).getData()[0]
Validate("Hydrogen-electron charge density at last timestep", Rho_electronH, 0.075)
