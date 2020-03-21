# ______________________________________________________________________________
#
# Validation script for the radiation in the collision of
# a GeV electron bunch with a counter-propagatin circularly polarized wave
#
# In this tests case, an electron bunch is initialized per radiation
# loss models at the same positions with an energy of 1 GeV near the right
# boundary of the box. They propagate to the left of the box where a circularly
# polarized laser plane wave is injected. This wave has an hyper-guassian
# profile of wavelength \lambda.
#
# Validation:
# - classical Landau-Lifshitz radiation model
# - corrected Landau-Lifshitz radiation model
# - Monte-Carlo radiation model
# - Niel stochastic radiation model
# ______________________________________________________________________________

import os, re, numpy as np, h5py
import happi

# Useful functions

def index_of_last_nonzero(lst):
    for i, value in enumerate(reversed(lst)):
        if value != 0:
            return len(lst)-i-1
    return -1


S = happi.Open(["./restart*"], verbose=False)

# List of relativistic pushers
radiation_list = ["LL","CLL","Niel","MC"]
thresholds = [0.1,0.1,0.2,0.15]

ukin_dict = {}
urad_dict = {}

# We load successively the particle track associated
# to each radiation algorithm
for i,radiation in enumerate(radiation_list):

  # Read scalar diagnostics
  ScalarUkinDiag = S.Scalar("Ukin_electron_"+radiation).get()
  ukin = np.array(ScalarUkinDiag["data"])
  times = np.array(ScalarUkinDiag["times"])

  ScalarUradDiag = S.Scalar("Urad_electron_"+radiation).get()
  urad = np.abs(np.array(ScalarUradDiag["data"]))

  utot = ukin+urad

  print(' Electron_'+radiation+':')
  print(' Final kinetic energy: '+str(ukin[-1]))
  print(' Maximum radiated energy: '+str(urad.max()))

  for it in range(53,len(ukin)):

    # Validation of the kinetic energy
    Validate("Kinetic energy evolution for {} at {}".format(radiation,it), ukin[it]/utot[it], ukin[it]/utot[it]*thresholds[i] )

    # Validation of the radiated energy
    Validate("Radiated energy evolution for {} at {}".format(radiation,it) , urad[it]/utot[it], urad[it]/utot[it]*thresholds[i] )

  # Validation of the total energy
  Validate("Total energy error (max - min)/uref for {} at {}".format(radiation,it) ,
               (utot.max() - utot.min())/utot[0], 1e-2)

  ukin_dict[radiation] = ukin
  urad_dict[radiation] = urad

# ______________________________________________________________________________
# Comparison CLL and Niel methods

urad_rel_err_Niel = abs(urad_dict["Niel"] - urad_dict["CLL"]) / urad_dict["CLL"].max()
ukin_rel_err_Niel = abs(ukin_dict["Niel"] - ukin_dict["CLL"]) / ukin_dict["CLL"][0]

print('')
print(' Comparison Landau-Lifshitz/Niel radiation model')
print(' Maximum relative error kinetic energy: {}'.format(ukin_rel_err_Niel.max()))
print(' Maximum relative error radiative energy: {}'.format(urad_rel_err_Niel.max()))

# Validation difference between Landau-Lifshitz and Niel methods
Validate("Relative error on the kinetic energy / ukin at t=0 (Niel/CLL)" , ukin_rel_err_Niel.max(), 0.05 )
Validate("Relative error on the radiative energy / urad max (Niel/CLL)" , urad_rel_err_Niel.max(), 0.05 )

# ______________________________________________________________________________
# Comparison CLL and MC methods

urad_rel_err_MC = abs(urad_dict["MC"] - urad_dict["CLL"]) / urad_dict["CLL"].max()
ukin_rel_err_MC = abs(ukin_dict["MC"] - ukin_dict["CLL"]) / ukin_dict["CLL"][0]

print('')
print(' Comparison Landau-Lifshitz/Monte-Carlo radiation model')
print(' Maximum relative error kinetic energy: {}'.format(ukin_rel_err_MC.max()))
print(' Maximum relative error radiative energy: {}'.format(urad_rel_err_MC.max()))

# Validation difference between Landau-Lifshitz and Monte-Carlo methods
Validate("Relative error on the kinetic energy / ukin at t=0 (MC/CLL)" , ukin_rel_err_MC.max(), 0.05 )
Validate("Relative error on the radiative energy / urad max (MC/CLL)" , urad_rel_err_MC.max(), 0.05 )

# ______________________________________________________________________________
# Checking of the particle binning

print("")
print(" Checking of the particle binning diagnostics")

minimal_iteration = 5000
maximal_iteration = 6200
period = 100
number_of_files = int((maximal_iteration - minimal_iteration)/period) + 1

average_chi = {}
integrated_chi = {}
max_chi = {}

# Loop over the species/radiation models
for i,radiation in enumerate(radiation_list):
    
    average_chi[radiation] = np.zeros(number_of_files)
    integrated_chi[radiation] = np.zeros(number_of_files)
    max_chi[radiation] = np.zeros(number_of_files)
        
    # Loop over the timesteps
    for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration,period)):
        
        # Weight
        weight_diag = S.ParticleBinning(diagNumber=i,timesteps=timestep).get()
        weight = np.array(weight_diag["data"][0])
        
        # Weight x chi
        weight_chi_diag = S.ParticleBinning(diagNumber=i+len(radiation_list),timesteps=timestep).get()
        weight_chi = np.array(weight_chi_diag["data"][0])
        
        # Chi distribution
        chi_distribution = S.ParticleBinning(diagNumber=i+2*len(radiation_list),timesteps=timestep).get()
        weight_values = np.array(chi_distribution["data"][0])
        chi_values = np.array(chi_distribution["chi"])
        log10_chi_values = np.log10(chi_values)
        delta = log10_chi_values[1] - log10_chi_values[0]
        bins =  np.power(10.,log10_chi_values + 0.5*delta) - np.power(10.,log10_chi_values - 0.5*delta)
                
        # Average chi value from chi spectrum
        total_weight = np.sum(weight_values)
        if (total_weight > 0):
            integrated_chi[radiation][itimestep] = np.sum (bins * weight_values)
            average_chi[radiation][itimestep] = integrated_chi[radiation][itimestep] / total_weight
            k = index_of_last_nonzero(weight_values)
            max_chi[radiation][itimestep] = chi_values[k]


print(" ---------------------------------------------------------|")
print(" Maximal quantum parameter                                |")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7} |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------------|")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,radiation in enumerate(radiation_list):
        line += " {0:.5f} |".format(max_chi[radiation][itimestep])
    print(line)
    for k,radiation in enumerate(radiation_list):
        Validate("Maximal quantum parameter for the {} model at iteration {}".format(radiation,timestep),max_chi[radiation][itimestep],max_chi[radiation][itimestep]*0.9)

print(" ---------------------------------------------------------|")
print(" Average quantum parameter                                |")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7} |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------------|")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,radiation in enumerate(radiation_list):
        line += " {0:.5f} |".format(average_chi[radiation][itimestep])
    print(line)
    for k,radiation in enumerate(radiation_list):
        # if 0, absolute error of 0.01
        if (average_chi[radiation][itimestep] == 0):
            Validate("Average quantum parameter for the {} model at iteration {}".format(radiation,timestep),average_chi[radiation][itimestep],0.01)
        # else validation with relative 150% error
        else:
            Validate("Average quantum parameter for the {} model at iteration {}".format(radiation,timestep),average_chi[radiation][itimestep],average_chi[radiation][itimestep]*0.1)
