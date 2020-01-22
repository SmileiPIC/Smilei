import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)

# Get some Smilei parameters
timestep = S.namelist.Main.timestep
simulation_time = S.namelist.Main.simulation_time

# List of species
species_list = ["electron","positron","photon"]
# First iteration used to validate the diag scalars
species_first_iteration = [150,150,150]
# Large relative error because of stochasticity
relative_error = 0.5

Scalar = {}

for ispecies,species in enumerate(species_list):
    name = "Ntot_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name][species_first_iteration[ispecies]:]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)
#
#     name = "Ukin_{}".format(species)
#     Scalar[name] = np.array(S.Scalar(name).getData())
#     for index,value in enumerate(Scalar[name][species_first_iteration[ispecies]:]):
#         Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)
#
#     name = "Dens_{}".format(species)
#     Scalar[name] = np.array(S.Scalar(name).getData())
#     for index,value in enumerate(Scalar[name]):
#         Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

# Energy _________________________________________________________________

# relative_error = 0.25
#
# gamma_spectrum_photon = S.ParticleBinning(diagNumber=2,timesteps=4021)
# data = np.array(gamma_spectrum_photon.getData()[0])
#
# for index,value in enumerate(data):
#     Validate("Gamma spectrum for photons at {}: ".format(index) , value, value*relative_error)
