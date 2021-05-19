# ______________________________________________________________________________
#
# Validation script for the benchmark tst1d_10_multiphoton_Breit_Wheeler
#
# Descrition:
# Interaction between a bunch of photons and a constant magnetic field.
#
# Purpose:
# During the interaction,the photons will progressively decay into pairs
# of electron-positron via the multiphoton Breit-Wheeler process.
# The radiation reaction of particles is not activated.
#
# Validation:
# - Propagation of macro-photons
# - multiphoton Breit-Wheeler pair creation process
#
# ______________________________________________________________________________

import os, re, numpy as np, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



dt = S.namelist.Main.timestep
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]

# ______________________________________________________________________________
# Read scalar diagnostics

print("")
print(" 1) Analyze of scalar diags")
print("")

times = np.array(S.Scalar("Ukin_electron").get()["times"])
ukin_electron = np.array(S.Scalar("Ukin_electron").get()["data"])
ukin_positron = np.array(S.Scalar("Ukin_positron").get()["data"])
ukin_photon = np.array(S.Scalar("Ukin_photon").get()["data"])

utot = ukin_electron + ukin_positron + ukin_photon

ntot_electron = np.array(S.Scalar("Ntot_electron").get()["data"])
ntot_positron = np.array(S.Scalar("Ntot_positron").get()["data"])
ntot_photon = np.array(S.Scalar("Ntot_photon").get()["data"])

print( ' Final electron energy / total energy: ',ukin_electron[-1] / utot[0] )
print( ' Final positron energy / total energy: ',ukin_positron[-1] / utot[0] )
print( ' Final photon energy / total energy: ',ukin_photon[-1] / utot[0] )
print( ' Final pair energy / total energy: ',UmBWpairs[-1] / utot[0] )
print( ' Maximal relative error total energy: ', max(abs(utot[:] - utot[0]))/utot[0] )

print( ' Final number of electrons:',ntot_electron[-1] )
print( ' Final number of positrons: ',ntot_positron[-1] )
print( ' Final number of photons: ',ntot_photon[-1] )

Validate("Electron kinetic energy evolution: ", ukin_electron/utot[0], 1e-2 )
Validate("Positron kinetic energy evolution: ", ukin_positron/utot[0], 2e-2 )
Validate("Photon kinetic energy evolution: ", ukin_photon/utot[0], 2e-2 )
Validate("Maximal relative error total energy: ", max(abs(utot[:] - utot[0]))/utot[0], 5e-3 )

# ______________________________________________________________________________
# Read energy spectrum

print("")
print(" 2) Analyze of the gamma distribution (particle binning)")
print("")

species_list = ["electron","positron","photon"]

integrated_gamma_spectrum = {}
max_gamma_spectrum = {}

for ispecies,species in enumerate(species_list):

    integrated_gamma_spectrum[species] = np.zeros(8)
    max_gamma_spectrum[species] = np.zeros(8)

    for diag_it in range(8):

        PartDiag = S.ParticleBinning(diagNumber=ispecies,timesteps = diag_it*50)
        gamma = np.array(PartDiag.get()["gamma"])
        density = np.array(PartDiag.get()["data"][0])

        log10_gamma = np.log10(gamma)
        delta = log10_gamma[1] - log10_gamma[0]
        bins =  np.power(10.,log10_gamma + 0.5*delta) - np.power(10.,log10_gamma - 0.5*delta)

        integrated_gamma_spectrum[species][diag_it] = np.sum(bins*density)
        imax = np.argmax(density)
        max_gamma_spectrum[species][diag_it] = gamma[imax]

print(" ---------------------------------------------------------")
print(" Integrated Gamma distribution                           |")
line = "                  |"
for ispecies,species in enumerate(species_list):
    line += " {0:<9}  |".format(species)
print(line)
print(" ---------------------------------------------------------")
# Loop over the timesteps
for diag_it in range(8):
    line = " Iteration {0:5d}  |".format(diag_it*50)
    for ispecies,species in enumerate(species_list):
        line += " {0:.4e} |".format(integrated_gamma_spectrum[species][diag_it])
    print(line)
# Validation
for diag_it in range(8):
    Validate("Integrated Gamma distribution for species electron at iteration {}".format(diag_it),integrated_gamma_spectrum["electron"][diag_it],integrated_gamma_spectrum["electron"][diag_it]*0.1)
    Validate("Integrated Gamma distribution for species positron at iteration {}".format(diag_it),integrated_gamma_spectrum["positron"][diag_it],integrated_gamma_spectrum["positron"][diag_it]*0.1)
    # Validate("Integrated Gamma distribution for species photon at iteration {}".format(diag_it),integrated_gamma_spectrum["photon"][diag_it],integrated_gamma_spectrum["photon"][diag_it]*0.2)
