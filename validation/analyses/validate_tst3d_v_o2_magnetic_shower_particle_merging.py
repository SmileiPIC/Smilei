import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)

# Get some Smilei parameters
timestep = S.namelist.Main.timestep
simulation_time = S.namelist.Main.simulation_time

# List of species
species_list = ["electron","positron","photon"]
#
species_first_index = [0,0,5]
# Large relative error because of stochasticity
relative_error = 0.1

Scalar = {}

for ispecies,species in enumerate(species_list):
    for name in ["Ntot_", "Ukin_", "Dens_"]:
        sc = np.array(S.Scalar(name+species).getData())
        Validate("Scalar "+name+species, sc[species_first_index[ispecies]:], relative_error, "relative_error")
        
        if name == "Ntot_":
            print("Number of {}s at beginning: {}".format(species,sc[0]))
            print("Number of {}s at the end: {}".format(species,sc[-1]))

# Energy _________________________________________________________________

relative_error = 0.1

gamma_spectrum_photon = S.ParticleBinning(diagNumber=2,timesteps=50)
data = np.array(gamma_spectrum_photon.getData()[0])
gamma = np.array(gamma_spectrum_photon.get()["gamma"])
integration = np.sum(data*gamma)*(gamma[1] - gamma[0])

Validate("Integrated photon gamma spectrum at the end: ",integration,integration*relative_error)
print("Integrated photon gamma spectrum at the end: {}.".format(integration))
#
# for index,value in enumerate(data):
#     Validate("Gamma spectrum for photons at {}: ".format(index) , value, value*relative_error)
