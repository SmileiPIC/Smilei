import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter1d
import happi
import sys

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(timesteps=1000).getData()[0]
Validate("Ey field at iteration 1000", Ey, 0.001)

# COMPARE THE Ez FIELD
Ez = S.Field.Field0.Ez(timesteps=1000).getData()[0]
Validate("Ez field at iteration 1000", Ez, 1e-9)

# hydrogen
charge_distribution = S.ParticleBinning.Diag0().getData()
charge_distribution /= charge_distribution[0].sum()
n = charge_distribution[0].size
mean_charge = [np.sum(d*np.arange(n)) for d in charge_distribution]
Validate("Hydrogen mean charge vs time", mean_charge, 0.03)

# carbon
charge_distribution = S.ParticleBinning.Diag1().getData()
charge_distribution /= charge_distribution[0].sum()
n = charge_distribution[0].size
mean_charge = [np.sum(d*np.arange(n)) for d in charge_distribution]
Validate("Carbon mean charge vs time", mean_charge, 0.03)

# SCALARS RELATED TO SPECIES
Validate("Scalar Dens_electron", S.Scalar.Dens_electron().getData(), 0.003)
Validate("Scalar Ntot_electron", S.Scalar.Ntot_electron().getData(), 100.)
Validate("Scalar Zavg_carbon"  , S.Scalar.Zavg_carbon  ().getData(), 0.2)

# TRACKING DIAGNOSTIC
d = S.TrackParticles("electron", axes=["Id","x","Wy"], timesteps=1000).getData()
keep = d["Id"] > 0
order = np.argsort(d["x"][keep])
Validate("Track electron x", d["x"][keep][order][::200], 1e-4)
Validate("Track electron Wy", gaussian_filter(maximum_filter1d(d["Wy"][keep][order],20),200)[::200], 1e-5)

# NEW PARTICLES DIAGNOSTIC
d = S.NewParticles.electron().get()
t = d["t"]
q = d["q"]
Validate("DiagNewParticles: number of particles", t.size, 5. )
tavg = [np.mean(t[q==i]) for i in [0,1,2,3]]
Validate("DiagNewParticles: time vs ionization state", tavg, 0.01 )

