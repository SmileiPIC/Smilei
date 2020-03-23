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

# ______________________________________________________________________________
# Read scalar diagnostics

print("")
print(" 1) Analyze of scalar diags")
print("")

ukin = {}

times = np.array(S.Scalar("Ukin_electron").get()["times"])
ukin["electron"] = np.array(S.Scalar("Ukin_electron").get()["data"])
ukin["positron"] = np.array(S.Scalar("Ukin_positron").get()["data"])
ukin["photon"] = np.array(S.Scalar("Ukin_photon").get()["data"])
urad = np.array(S.Scalar("Urad").get()["data"])

ukin["tot"] = ukin["electron"] + ukin["positron"] + ukin["photon"] + urad

ntot = {}

ntot["electron"] = np.array(S.Scalar("Ntot_electron").get()["data"])
ntot["positron"] = np.array(S.Scalar("Ntot_positron").get()["data"])
ntot["photon"] = np.array(S.Scalar("Ntot_photon").get()["data"])

ntot["tot"] = ntot["electron"] + ntot["positron"] + ntot["photon"]
ntot["particles"] = ntot["electron"] + ntot["positron"]

print(' Final electron energy / initial electron energy: '+str(ukin["electron"][-1] / ukin["electron"][0]))
print(' Final positron energy / initial electron energy: '+str(ukin["positron"][-1] / ukin["electron"][0]))
print(' Final photon energy / initial electron energy: '+str(ukin["photon"][-1] / ukin["electron"][0]))
print(' Final radiated energy / initial electron energy: '+str(urad[-1] / ukin["electron"][0]))

print(' Final number of electrons: {}'.format(ntot["electron"][-1]))
print(' Final number of positrons: '+str(ntot["positron"][-1]))
print(' Final number of photons: '+str(ntot["photon"][-1]))

print("")
print(" ---------------------------------------------------------------------------|")
print(" Diag scalars (Kinetic energy)                                              |")
print(" iteration | ukin e-    | ukin e+    | ukin ph    | urad       | ukin tot   |")
print(" ---------------------------------------------------------------------------|")

for it,time in enumerate(times):
    if (it > 39):
        print(" {0:5d}     | {1:.4e} | {2:.4e} | {3:.4e} | {4:.4e} | {5:.4e} |".format(it*100,ukin["electron"][it],ukin["positron"][it],ukin["photon"][it],urad[it], ukin["tot"][it]))

thresholds = {}
thresholds["points"] = np.array([0.,10.,100.,1000.])
thresholds["factor"] = np.array([1e9,1.,0.5,0.2,0.1])

for it,time in enumerate(times):
    if (it > 39):
        # Validation of the kinetic energy
        Validate("Electron kinetic energy at {}".format(it), ukin["electron"][it]/ukin["tot"][it], adaptive_error(ukin["electron"][it]/ukin["tot"][it],ntot["electron"][it],thresholds))
        Validate("Positron kinetic energy at {}".format(it), ukin["positron"][it]/ukin["tot"][it], adaptive_error(ukin["positron"][it]/ukin["tot"][it],ntot["positron"][it],thresholds))
        Validate("Photon kinetic energy at {}".format(it), ukin["photon"][it]/ukin["tot"][it], adaptive_error(ukin["photon"][it]/ukin["tot"][it],ntot["photon"][it],thresholds))
        Validate("Radiated energy at {}".format(it), urad[it]/ukin["tot"][it], adaptive_error(urad[it]/ukin["tot"][it],ntot["particles"][it],thresholds))

print(" ---------------------------------------------------------------------------|")
print(" Diag scalars (Number of particles)                                         |")
print(" iteration | num e-     | num e+     | num  ph    | num tot    |            |")
print(" ---------------------------------------------------------------------------|")

for it,time in enumerate(times):
    if (it > 39):
        
        print(" {0:5d}     | {1:.4e} | {2:.4e} | {3:.4e} | {4:.4e} |".format(it*100,ntot["electron"][it],ntot["positron"][it],ntot["photon"][it],ntot["tot"][it]))

thresholds = {}
thresholds["points"] = np.array([0.,10,100,1000])
thresholds["factor"] = np.array([1e9, 1.,0.5,0.2,0.1])

for it,time in enumerate(times):
    if (it > 39):

        Validate("Number of electrons at {}".format(it), ntot["electron"][it], adaptive_error(ntot["electron"][it],ntot["electron"][it],thresholds) )
        Validate("Number of positrons at {}".format(it), ntot["positron"][it], adaptive_error(ntot["positron"][it],ntot["positron"][it],thresholds) )
        Validate("Number of photons at {}".format(it), ntot["photon"][it], adaptive_error(ntot["photon"][it],ntot["photon"][it],thresholds) )

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
    
thresholds = {}
thresholds["points"] = np.array([0.,10,100])
thresholds["factor"] = np.array([1e9,1.,0.5,0.5])
    
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    for ispecies,species in enumerate(species_list):
        Validate("Maximal quantum parameter for the {} model at iteration {}".format(species,timestep),max_chi[species][itimestep],adaptive_error(max_chi[species][itimestep],ntot[species][itimestep],thresholds))

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
thresholds["points"] = np.array([0.,10,100,1000])
thresholds["factor"] = np.array([1e9,1.,0.5,0.2,0.1])

for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    for ispecies,species in enumerate(species_list):
        Validate("Average quantum parameter for the {} model at iteration {}".format(species,timestep),average_chi[species][itimestep],adaptive_error(average_chi[species][itimestep],ntot[species][itimestep],thresholds))

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
    
thresholds = {}
thresholds["points"] = np.array([0.,10,100])
thresholds["factor"] = np.array([1e9,1.,0.5,0.5])
    
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    for ispecies,species in enumerate(species_list):
        Validate("Maximal gamma for the {} model at iteration {}".format(species,timestep),max_gamma[species][itimestep],adaptive_error(max_gamma[species][itimestep],ntot[species][itimestep],thresholds))

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
thresholds["factor"] = np.array([1e9,1.,0.5,0.2,0.1])

for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    for ispecies,species in enumerate(species_list):
        Validate("Average gamma for the {} model at iteration {}".format(species,timestep),average_gamma[species][itimestep],adaptive_error(average_gamma[species][itimestep],ntot[species][itimestep],thresholds))
