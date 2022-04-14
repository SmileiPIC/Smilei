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
    name = "Ntot_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name][species_first_index[ispecies]:]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

    print("Number of {}s at beginning: {}".format(species,Scalar[name][0]))
    print("Number of {}s at the end: {}".format(species,Scalar[name][-1]))

    name = "Ukin_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name][species_first_index[ispecies]:]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

    name = "Dens_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name][species_first_index[ispecies]:]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

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
