import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)

# Get some Smilei parameters
timestep = S.namelist.Main.timestep
simulation_time = S.namelist.Main.simulation_time

diag_every = int(simulation_time / timestep)

species_list = ["eon1", "pon1", "eon2", "pon2"]

Scalar = {}

for species in species_list:
    name = "Ntot_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    Validate("Scalar {}".format(name) , Scalar[name], 0.01)
    
    name = "Ukin_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    Validate("Scalar {}".format(name) , Scalar[name], 0.01)

    name = "Dens_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    Validate("Scalar {}".format(name) , Scalar[name], 0.01)

# Energy _________________________________________________________________

for i in range(4):
    particle_binning_initial = S.ParticleBinning(diagNumber=i,timesteps=0)
    data_initial = np.array(particle_binning_initial.getData()[0])
    
    particle_binning = S.ParticleBinning(diagNumber=i,timesteps=diag_every)
    data_final = np.array(particle_binning.getData()[0])
    
    Validate("Gamma spectrum for {}".format(species_list[i]) , data_final, 0.05)
    
    sum_initial = np.sum(data_initial)
    sum_final = np.sum(data_final)
    
    error = (np.abs(np.subtract(data_initial / sum_initial, data_final / sum_final)))
    
    print(' Gamma spectrum max error for {}: {}'.format(species_list[i], np.max(error)))
    
    Validate("Gamma spectrum max error for {}".format(species_list[i]) , np.max(error), 0.01)
