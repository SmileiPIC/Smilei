import os, re, numpy as np
import happi
from scipy.signal import butter, filtfilt
b, a = butter(5, 0.2, btype='low', analog=False)


S = happi.Open("./restart*", verbose=False)

eon_spectrum = S.ParticleBinning.Diag2().get()
ekin = eon_spectrum["ekin"]
eon_spectrum = np.mean(eon_spectrum["data"], axis=0)
eon_spectrum_filt = filtfilt(b, a, eon_spectrum)
# # theory
# Te = S.namelist.Species["eon"].temperature[0]
# factor = S.namelist.Species["eon"].number_density.xplateau / S.namelist.Main.grid_length[0]
# theoretical_spectrum = factor*2./Te * (ekin/np.pi/Te)**0.5 * np.exp(-ekin/Te)
# plt.plot(ekin, eon_spectrum_filt, '.-')
# plt.plot(ekin, theoretical_spectrum, '-')
Validate("Electron spectrum", eon_spectrum_filt, 20. )


rho = S.Field.Field0.Rho_ion(timesteps=11800).getData()[0]
rho_filt = filtfilt(b, a, rho)
Validate("Final ion profile", rho_filt[::10], 0.15)
