import os, re, numpy as np, math
import happi

# ______________________________________________________________________________
# Useful functions

def index_of_last_nonzero(lst):
    for i, value in enumerate(reversed(lst)):
        if value != 0:
            return len(lst)-i-1
    return -1

def adaptive_error(value, number_of_points, thresholds):
    """
    This function return an error that depends on the statistic.
    """
    
    # We eliminate the case where there is no data
    if (number_of_points <= 0):
        return thresholds["factor"][0]
    
    flag = True
    i_threshold = 0
    while(flag):
        if (number_of_points < thresholds["points"][i_threshold]):
            flag = False
        else:
            i_threshold+=1
        if (i_threshold >= np.size(thresholds["points"])):
            flag = False
    if ((i_threshold == 0) or (i_threshold >= np.size(thresholds["points"]))):
        return thresholds["factor"][i_threshold]*value
    else:
        i_threshold -= 1
        d = (number_of_points - thresholds["points"][i_threshold]) / (thresholds["points"][i_threshold+1] - thresholds["points"][i_threshold])
        return value*(thresholds["factor"][i_threshold]*(1-d) + d*thresholds["factor"][i_threshold+1])

# ______________________________________________________________________________

S = happi.Open(["./restart*"], verbose=False)

print("")
print(" 1) Analyze of Scalars")
print("")

# Get some Smilei parameters
timestep = S.namelist.Main.timestep
simulation_time = S.namelist.Main.simulation_time

# List of species
species_list = ["electron","positron","photon"]
# First iteration used to validate the diag scalars
species_first_iteration = [200,200,200]
# Large relative error because of stochasticity
relative_error = 0.5

Scalar = {}

Scalar["times"] = np.array(S.Scalar("Dens_electron").get()["times"])

for ispecies,species in enumerate(species_list):
    name = "Ntot_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    # for index,value in enumerate(Scalar[name][species_first_iteration[ispecies]:]):
    #     Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

    name = "Ukin_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    # for index,value in enumerate(Scalar[name][species_first_iteration[ispecies]:]):
    #     Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

    name = "Dens_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    # for index,value in enumerate(Scalar[name]):
    #     Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

Scalar["Ntot"] = Scalar["Ntot_electron"] + Scalar["Ntot_positron"] + Scalar["Ntot_photon"]
Scalar["Ukin_tot"] = Scalar["Ukin_electron"] + Scalar["Ukin_positron"] + Scalar["Ukin_photon"]

scalar_period = 100

print(" ---------------------------------------------------------------------------|")
print(" Diag scalars (Number of particles)                                         |")
print(" iteration | num e-     | num e+     | num  ph    | num tot    |            |")
print(" ---------------------------------------------------------------------------|")

for it,time in enumerate(Scalar["times"]):
    if (it*scalar_period >= 2000 and it*scalar_period <= 6500 and np.mod(it,2) == 0):
        print(" {0:5d}     | {1:.4e} | {2:.4e} | {3:.4e} | {4:.4e} |".format(it*scalar_period,Scalar["Ntot_electron"][it],Scalar["Ntot_positron"][it],Scalar["Ntot_photon"][it],Scalar["Ntot"][it]))

thresholds = {}
thresholds["points"] = np.array([0. , 10, 100, 1000, 10000])
thresholds["factor"] = np.array([1e9,1.5, 0.6, 0.35,   0.3, 0.25])

for it,time in enumerate(Scalar["times"]):
    if (it*scalar_period >= 2000 and it*scalar_period <= 6500 and np.mod(it,2) == 0):
        Validate("Number of electrons at {}".format(it*scalar_period), Scalar["Ntot_electron"][it], adaptive_error(Scalar["Ntot_electron"][it],Scalar["Ntot_electron"][it],thresholds))
        Validate("Number of positrons at {}".format(it*scalar_period), Scalar["Ntot_positron"][it], adaptive_error(Scalar["Ntot_positron"][it],Scalar["Ntot_positron"][it],thresholds))
        Validate("Number of photons at {}".format(it*scalar_period), Scalar["Ntot_photon"][it], adaptive_error(Scalar["Ntot_photon"][it],Scalar["Ntot_photon"][it],thresholds))

print(" ---------------------------------------------------------------------------|")
print(" Diag scalars (Kinetic energy)                                              |")
print(" iteration | electrons  | positron   | photon     | tot        |            |")
print(" ---------------------------------------------------------------------------|")

for it,time in enumerate(Scalar["times"]):
    if (it*scalar_period >= 2000 and it*scalar_period <= 6500 and np.mod(it,2) == 0):
        print(" {0:5d}     | {1:.4e} | {2:.4e} | {3:.4e} | {4:.4e} |".format(it*scalar_period,Scalar["Ukin_electron"][it],Scalar["Ukin_positron"][it],Scalar["Ukin_photon"][it],Scalar["Ukin_tot"][it]))

thresholds = {}
thresholds["points"] = np.array([0. ,10 ,100 , 1000, 10000])
thresholds["factor"] = np.array([1e9, 1.,0.6 ,  0.4,   0.3, 0.2])

for it,time in enumerate(Scalar["times"]):
    if (it*scalar_period >= 2000 and it*scalar_period <= 6500 and np.mod(it,2) == 0):
        Validate("Electron energy at {}".format(it*scalar_period), Scalar["Ukin_electron"][it], adaptive_error(Scalar["Ukin_electron"][it],Scalar["Ntot_electron"][it],thresholds))
        Validate("Positron energy at {}".format(it*scalar_period), Scalar["Ukin_positron"][it], adaptive_error(Scalar["Ukin_positron"][it],Scalar["Ntot_positron"][it],thresholds))

thresholds = {}
thresholds["points"] = np.array([0. ,10 ,100 , 1000, 10000])
thresholds["factor"] = np.array([1e9, 1.,0.6 , 0.45,  0.35, 0.25])

for it,time in enumerate(Scalar["times"]):
    if (it*scalar_period >= 2000 and it*scalar_period <= 6500 and np.mod(it,2) == 0):
        Validate("Photon energy at {}".format(it*scalar_period), Scalar["Ukin_photon"][it], adaptive_error(Scalar["Ukin_photon"][it],Scalar["Ntot_photon"][it],thresholds))

# Energy _________________________________________________________________

species_list = ["electron","positron","photon"]

minimal_iteration = 2000
maximal_iteration = 6500
period = 500
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
        chi_distribution = S.ParticleBinning(diagNumber=ispecies+len(species_list),timesteps=timestep).get()
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
#    for ispecies,species in enumerate(species_list):
#        Validate("Maximal quantum parameter for the {} model at iteration {}".format(species,timestep),max_chi[species][itimestep],max_chi[species][itimestep]*0.5)
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
    
thresholds = {}
thresholds["points"] = np.array([0.,10,1000])
thresholds["factor"] = np.array([1e9, 1.,0.5,0.3])
    
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    Validate("Average quantum parameter for electrons at iteration {}".format(timestep),average_chi["electron"][itimestep],adaptive_error(average_chi["electron"][itimestep],Scalar["Ntot_electron"][itimestep],thresholds))
    Validate("Average quantum parameter for positrons at iteration {}".format(timestep),average_chi["positron"][itimestep],adaptive_error(average_chi["positron"][itimestep],Scalar["Ntot_positron"][itimestep],thresholds))
    Validate("Average quantum parameter for photons at iteration {}".format(timestep),average_chi["photon"][itimestep],adaptive_error(average_chi["photon"][itimestep],Scalar["Ntot_photon"][itimestep],thresholds))

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
        gamma_distribution = S.ParticleBinning(diagNumber=ispecies,timesteps=timestep).get()
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
    
# thresholds = {}
# thresholds["points"] = np.array([0.,10,100])
# thresholds["factor"] = np.array([1e9, 1.,0.9,0.9])
#
# for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
#     for ispecies,species in enumerate(species_list):
#         Validate("Maximal gamma for the {} model at iteration {}".format(species,timestep),max_gamma[species][itimestep],adaptive_error(max_gamma[species][itimestep],Scalar["Ntot_{}".format(species)][itimestep],thresholds))

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

thresholds = {}
thresholds["points"] = np.array([0.,10,100,1000])
thresholds["factor"] = np.array([1e9, 1.,0.5,0.2,0.1])

for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    Validate("Average gamma for the electrons at iteration {}".format(timestep),average_gamma["electron"][itimestep],adaptive_error(average_gamma["electron"][itimestep],Scalar["Ntot_electron"][itimestep],thresholds))
    Validate("Average gamma for the positrons at iteration {}".format(timestep),average_gamma["positron"][itimestep],adaptive_error(average_gamma["positron"][itimestep],Scalar["Ntot_positron"][itimestep],thresholds))
    Validate("Average gamma for the photons at iteration {}".format(timestep),average_gamma["photon"][itimestep],adaptive_error(average_gamma["photon"][itimestep],Scalar["Ntot_photon"][itimestep],thresholds))

# relative_error = 0.25
#
# gamma_spectrum_photon = S.ParticleBinning(diagNumber=2,timesteps=4021)
# data = np.array(gamma_spectrum_photon.getData()[0])
# gamma = np.array(gamma_spectrum_photon.get()["gamma"])
# integration = np.sum(data*gamma)*(gamma[1] - gamma[0])
#
# Validate("Integrated photon gamma spectrum at 4021: ",integration,integration*relative_error)
# print("Integrated photon gamma spectrum at 4021: {}.".format(integration))
#
# for index,value in enumerate(data):
#     Validate("Gamma spectrum for photons at {}: ".format(index) , value, value*relative_error)
