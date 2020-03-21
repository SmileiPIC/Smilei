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

# ______________________________________________________________________________
# Useful functions

def index_of_last_nonzero(lst):
    for i, value in enumerate(reversed(lst)):
        if value != 0:
            return len(lst)-i-1
    return -1

# ______________________________________________________________________________

S = happi.Open(["./restart*"], verbose=False)

# ______________________________________________________________________________
# Read scalar diagnostics

print("")
print(" 1) Analyze of scalar diags")
print("")

times = np.array(S.Scalar("Ukin_electron").get()["times"])
ukin_electron = np.array(S.Scalar("Ukin_electron").get()["data"])
ukin_positron = np.array(S.Scalar("Ukin_positron").get()["data"])
ukin_photon = np.array(S.Scalar("Ukin_photon").get()["data"])
urad = np.array(S.Scalar("Urad").get()["data"])

utot = ukin_electron + ukin_positron + ukin_photon + urad

ntot_electron = np.array(S.Scalar("Ntot_electron").get()["data"])
ntot_positron = np.array(S.Scalar("Ntot_positron").get()["data"])
ntot_photon = np.array(S.Scalar("Ntot_photon").get()["data"])

ntot = ntot_electron + ntot_positron + ntot_photon

print(' Final electron energy / initial electron energy: '+str(ukin_electron[-1] / ukin_electron[0]))
print(' Final positron energy / initial electron energy: '+str(ukin_positron[-1] / ukin_electron[0]))
print(' Final photon energy / initial electron energy: '+str(ukin_photon[-1] / ukin_electron[0]))
print(' Final radiated energy / initial electron energy: '+str(urad[-1] / ukin_electron[0]))

print(' Final number of electrons: {}'.format(ntot_electron[-1]))
print(' Final number of positrons: '+str(ntot_positron[-1]))
print(' Final number of photons: '+str(ntot_photon[-1]))

print("")
print(" ---------------------------------------------------------------------------|")
print(" Diag scalars (Kinetic energy)                                              |")
print(" iteration | ukin e-    | ukin e+    | ukin ph    | urad       | ukin tot   |")
print(" ---------------------------------------------------------------------------|")

for it,time in enumerate(times):

    if (it > 39):

        print(" {0:5d}     | {1:.4e} | {2:.4e} | {3:.4e} | {4:.4e} | {5:.4e} |".format(it,ukin_electron[it],ukin_positron[it],ukin_photon[it],urad[it], utot[it]))

        # Validation of the kinetic energy
        Validate("Electron kinetic energy at {}".format(it), ukin_electron[it]/utot[it], ukin_electron[it]/utot[it]*0.1)
        Validate("Positron kinetic energy at {}".format(it), ukin_positron[it]/utot[it], ukin_positron[it]/utot[it]*0.1)
        Validate("Photon kinetic energy at {}".format(it), ukin_photon[it]/utot[it], ukin_photon[it]/utot[it]*0.1)
        Validate("Radiated energy at {}".format(it), urad[it]/utot[it], urad[it]/utot[it]*0.1)

print(" ---------------------------------------------------------------------------|")
print(" Diag scalars (Number of particles)                                         |")
print(" iteration | num e-     | num e+     | num  ph    | num tot    |            |")
print(" ---------------------------------------------------------------------------|")

for it,time in enumerate(times):

    if (it > 39):
        
        print(" {0:5d}     | {1:.4e} | {2:.4e} | {3:.4e} | {4:.4e} |".format(it,ntot_electron[it],ntot_positron[it],ntot_photon[it],ntot[it]))

        Validate("Number of electrons at {}".format(it), ntot_electron[it], ntot_electron[it]*0.1 )
        Validate("Number of positrons at {}".format(it), ntot_positron[it], ntot_positron[it]*0.1 )
        Validate("Number of photons at {}".format(it), ntot_photon[it], ntot_photon[it]*0.1 )

# ______________________________________________________________________________
# Read particle binning

species_list = ["electron","positron","photon"]

minimal_iteration = 4000
maximal_iteration = 5200
period = 100
number_of_files = int((maximal_iteration - minimal_iteration)/period) + 1

print("")
print(" 2) Analyze of chi spectrum")
print("")

average_chi = {}
integrated_chi = {}
max_chi = {}

# Loop over the species/radiation models
for ispecies,species in enumerate(species_list):
    
    average_chi[species] = np.zeros(number_of_files)
    integrated_chi[species] = np.zeros(number_of_files)
    max_chi[species] = np.zeros(number_of_files)
        
    # Loop over the timesteps
    for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration,period)):

        # Chi distribution
        chi_distribution = S.ParticleBinning(diagNumber=ispecies+1*len(species_list),timesteps=timestep).get()
        weight_values = np.array(chi_distribution["data"][0])
        chi_values = np.array(chi_distribution["chi"])
        log10_chi_values = np.log10(chi_values)
        delta = log10_chi_values[1] - log10_chi_values[0]
        bins =  np.power(10.,log10_chi_values + 0.5*delta) - np.power(10.,log10_chi_values - 0.5*delta)

        # Average chi value from chi spectrum
        total_weight = np.sum(weight_values)
        integrated_chi[species][itimestep] = np.sum (bins * weight_values)
        if (total_weight > 0):
            average_chi[species][itimestep] = integrated_chi[species][itimestep] / total_weight
            k = index_of_last_nonzero(weight_values)
            max_chi[species][itimestep] = chi_values[k]

print(" ---------------------------------------------------------------------------|")
print(" Maximal quantum parameter                                                  |")
print(" iteration | electrons  | positrons  | photons    |                         |")
print(" ---------------------------------------------------------------------------|")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    line = "  {0:5d}    |".format(timestep)
    for ispecies,species in enumerate(species_list):
        line += " {0:.4e} |".format(max_chi[species][itimestep])
    print(line)
    for ispecies,species in enumerate(species_list):
        Validate("Maximal quantum parameter for the {} model at iteration {}".format(species,timestep),max_chi[species][itimestep],max_chi[species][itimestep]*0.5)
print(" ---------------------------------------------------------------------------|")
print(" Average quantum parameter                                                  |")
print(" iteration | electrons  | positrons  | photons    |                         |")
print(" ---------------------------------------------------------------------------|")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    line = "  {0:5d}    |".format(timestep)
    for ispecies,species in enumerate(species_list):
        line += " {0:.4e} |".format(average_chi[species][itimestep])
    print(line)
    for ispecies,species in enumerate(species_list):
        Validate("Average quantum parameter for the {} model at iteration {}".format(species,timestep),average_chi[species][itimestep],average_chi[species][itimestep]*0.5)

print("")
print(" 3) Analyze of gamma spectrum")
print("")

average_gamma = {}
integrated_gamma = {}
max_gamma = {}

# Loop over the species/radiation models
for ispecies,species in enumerate(species_list):
    
    average_gamma[species] = np.zeros(number_of_files)
    integrated_gamma[species] = np.zeros(number_of_files)
    max_gamma[species] = np.zeros(number_of_files)
        
    # Loop over the timesteps
    for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration,period)):

        # Chi distribution
        gamma_distribution = S.ParticleBinning(diagNumber=ispecies+2*len(species_list),timesteps=timestep).get()
        weight_values = np.array(gamma_distribution["data"][0])
        gamma_values = np.array(gamma_distribution["gamma"])
        log10_gamma_values = np.log10(gamma_values)
        delta = log10_gamma_values[1] - log10_gamma_values[0]
        bins =  np.power(10.,log10_gamma_values + 0.5*delta) - np.power(10.,log10_gamma_values - 0.5*delta)

        total_weight = np.sum(weight_values)
        integrated_gamma[species][itimestep] = np.sum (bins * weight_values)
        if (total_weight > 0):
            average_gamma[species][itimestep] = integrated_gamma[species][itimestep] / total_weight
            k = index_of_last_nonzero(weight_values)
            max_gamma[species][itimestep] = gamma_values[k]
            
print(" ---------------------------------------------------------------------------|")
print(" Maximal gamma                                                              |")
print(" iteration | electrons  | positrons  | photons    |                         |")
print(" ---------------------------------------------------------------------------|")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    line = "  {0:5d}    |".format(timestep)
    for ispecies,species in enumerate(species_list):
        line += " {0:.4e} |".format(max_gamma[species][itimestep])
    print(line)
    for ispecies,species in enumerate(species_list):
        Validate("Maximal gamma for the {} model at iteration {}".format(species,timestep),max_gamma[species][itimestep],max_gamma[species][itimestep]*0.5)
print(" ---------------------------------------------------------------------------|")
print(" Average gamma                                                              |")
print(" iteration | electrons  | positrons  | photons    |                         |")
print(" ---------------------------------------------------------------------------|")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    line = "  {0:5d}    |".format(timestep)
    for ispecies,species in enumerate(species_list):
        line += " {0:.4e} |".format(average_gamma[species][itimestep])
    print(line)
    for ispecies,species in enumerate(species_list):
        Validate("Average gamma for the {} model at iteration {}".format(species,timestep),average_gamma[species][itimestep],average_gamma[species][itimestep]*0.1)
