import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)

# Get some Smilei parameters
timestep = S.namelist.Main.timestep
simulation_time = S.namelist.Main.simulation_time

diag_every = int(simulation_time / timestep)

species_list = ["eon1", "pon1", "eon2", "pon2"]

Scalar = {}

relative_error = 0.05

for species in species_list:
    name = "Ntot_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)
    
    name = "Ukin_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

    name = "Dens_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

# Energy _________________________________________________________________

for i in range(4):
    particle_binning_initial = S.ParticleBinning(diagNumber=i,timesteps=0)
    data_initial = np.array(particle_binning_initial.getData()[0])
    
    particle_binning = S.ParticleBinning(diagNumber=i,timesteps=diag_every)
    data_final = np.array(particle_binning.getData()[0])
    
    for index,value in enumerate(data_final):
        Validate("Gamma spectrum for {} at {}".format(species_list[i],index) , value, value*relative_error)
    
    sum_initial = np.sum(data_initial)
    sum_final = np.sum(data_final)
    
    error = (np.abs(np.subtract(data_initial / sum_initial, data_final / sum_final)))
    
    print(' Gamma spectrum max error for {}: {}'.format(species_list[i], np.max(error)))
    
    Validate("Gamma spectrum max error for {}".format(species_list[i]) , np.max(error), np.max(error)*relative_error)
