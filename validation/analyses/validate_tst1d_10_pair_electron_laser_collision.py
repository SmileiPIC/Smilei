# ______________________________________________________________________________
#
# Validation script for the benchmark tst1d_10_pair_counter_prop
#
# Descrition:
# Head-on collision between an electron bunch and a counter-propagating laser wave.
# The electron beam has an initial energy of 4 GeV.
# The laser has an intensity close to 10^23 W/cm^2
#
# Purpose:
# During the interaction,the electrons will radiate high-energy photons
# that in turn will decay onto pairs via the multiphoton Breit-Wheeler process.
#
# Validation:
# - Radiation reaction Monte-Carlo process
# - Emission and propagation of macro-photons
# - multiphoton Breit-Wheeler pair creation process
#
# ______________________________________________________________________________

import os, re, numpy as np, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



# ______________________________________________________________________________
# Read scalar diagnostics

times = np.array(S.Scalar("Ukin_electron").get()["times"])
ukin_electron = np.array(S.Scalar("Ukin_electron").get()["data"])
ukin_positron = np.array(S.Scalar("Ukin_positron").get()["data"])
ukin_photon = np.array(S.Scalar("Ukin_photon").get()["data"])
urad = np.array(S.Scalar("Urad").get()["data"])

ntot_electron = np.array(S.Scalar("Ntot_electron").get()["data"])
ntot_positron = np.array(S.Scalar("Ntot_positron").get()["data"])
ntot_photon = np.array(S.Scalar("Ntot_photon").get()["data"])

print(' Final electron energy / initial electron energy: '+str(ukin_electron[-1] / ukin_electron[0]))
print(' Final positron energy / initial electron energy: '+str(ukin_positron[-1] / ukin_electron[0]))
print(' Final photon energy / initial electron energy: '+str(ukin_photon[-1] / ukin_electron[0]))
print(' Final radiated energy / initial electron energy: '+str(urad[-1] / ukin_electron[0]))

print(' Final number of electrons:'+str(ntot_electron[-1]))
print(' Final number of positrons: '+str(ntot_positron[-1]))
print(' Final number of photons: '+str(ntot_photon[-1]))

# Validation of the kinetic energy
Validate("Electron kinetic energy evolution: ", ukin_electron/ukin_electron[0], 3e-2 )
Validate("Positron kinetic energy evolution: ", ukin_positron/ukin_electron[0], 3e-2 )
Validate("Photon kinetic energy evolution: ", ukin_photon/ukin_electron[0], 3e-2 )
Validate("Radiated energy evolution: ", urad/ukin_electron[0], 3e-2 )

Validate("Evolution of the number of electrons: ", ntot_electron, 120 )
Validate("Evolution of the number of positrons: ", ntot_positron, 120 )
Validate("Evolution of the number of photons: ", ntot_photon, 3000 )
