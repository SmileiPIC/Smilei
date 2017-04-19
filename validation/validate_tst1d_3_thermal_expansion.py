import os, re, numpy as np
from Smilei import *
from scipy.signal import butter, filtfilt

S = Smilei(".", verbose=False)

eon_spectrum = S.ParticleDiagnostic.Diag2().get()
ekin = eon_spectrum["ekin"]
eon_spectrum = np.mean(eon_spectrum["data"], axis=0)
Te = S.namelist.Species["eon"].temperature[0]
factor = S.namelist.Species["eon"].nb_density.xplateau / S.namelist.Main.sim_length[0]
theoretical_spectrum = factor*2./Te * (ekin/np.pi/Te)**0.5 * np.exp(-ekin/Te)
ok = (np.abs(eon_spectrum-theoretical_spectrum) < 20.).all()
Validate("Error on electron spectrum", ok )


rho = S.Field.Field0.Rho_ion(timesteps=11800).getData()[0]
b, a = butter(5, 0.2, btype='low', analog=False)
rho_filt = filtfilt(b, a, rho)
Validate("Final ion profile", rho_filt[::10], 0.15)