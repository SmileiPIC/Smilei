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
thresholds = [0.3,0.05,0.05,0.07]

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
  urad = np.array(ScalarUradDiag["data"])

  utot = ukin+urad

  print(' Electron_'+radiation+':')
  print(' Final kinetic energy: '+str(ukin[-1]))
  print(' Maximum radiated energy: '+str(urad.max()))

  # Validation of the kinetic energy
  Validate("Kinetic energy evolution: ", ukin/ukin[0], thresholds[i] )

  # Validation of the radiated energy
  Validate("Radiated energy evolution: " , urad/ukin[0], thresholds[i] )

  # Validation of the total energy
  Validate("Total energy error (max - min)/uref: " ,
           (utot.max() - utot.min())/utot[0], 1e-2)

  ukin_dict[radiation] = ukin
  urad_dict[radiation] = urad

# ______________________________________________________________________________
# Comparison CLL and Niel methods

urad_rel_err_Niel = abs(urad_dict["Niel"] - urad_dict["CLL"]) / urad_dict["CLL"].max()
ukin_rel_err_Niel = abs(ukin_dict["Niel"] - ukin_dict["CLL"]) / ukin_dict["CLL"][0]

print('')
print(' Comparison Landau-Lifshitz/Niel radiation model')
print(' Maximum relative error kinetic energy'+str(ukin_rel_err_Niel.max()))
print (' Maximum relative error radiative energy'+str(urad_rel_err_Niel.max()))

# Validation difference between Landau-Lifshitz and Niel methods
Validate("Relative error on the kinetic energy / ukin at t=0: " , ukin_rel_err_Niel.max(), 0.05 )
Validate("Relative error on the radiative energy / urad max " , urad_rel_err_Niel.max(), 0.05 )

# ______________________________________________________________________________
# Comparison CLL and MC methods

urad_rel_err_MC = abs(urad_dict["MC"] - urad_dict["CLL"]) / urad_dict["CLL"].max()
ukin_rel_err_MC = abs(ukin_dict["MC"] - ukin_dict["CLL"]) / ukin_dict["CLL"][0]

print('')
print(' Comparison Landau-Lifshitz/Monte-Carlo radiation model')
print(' Maximum relative error kinetic energy'+str(ukin_rel_err_MC.max()))
print(' Maximum relative error radiative energy'+str(urad_rel_err_MC.max()))

# Validation difference between Landau-Lifshitz and Monte-Carlo methods
Validate("Relative error on the kinetic energy / ukin at t=0: " , ukin_rel_err_MC.max(), 0.05 )
Validate("Relative error on the radiative energy / urad max " , urad_rel_err_MC.max(), 0.05 )

# ______________________________________________________________________________
# Checking of the particle binning

print("")
print(" Checking of the particle binning diagnostics")

minimal_iteration = 5000
maximal_iteration = 6300
period = 100
number_of_files = int((maximal_iteration - minimal_iteration)/period) + 1

chi_max = np.zeros([number_of_files,len(radiation_list)])
chi_ave = np.zeros([number_of_files,len(radiation_list)])

# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration,period)):

    # Loop over the species/radiation models
    for i,radiation in enumerate(radiation_list):
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
        # Average chi value
        total_weight = np.sum(weight_values)
        if (total_weight > 0):
            chi_ave[itimestep,i] = np.sum (chi_values * weight_values) / np.sum(weight_values)
            k = index_of_last_nonzero(weight_values)
            chi_max[itimestep,i] = chi_values[k]


print(" ---------------------------------------------------------|")
print(" Maximal quantum parameter")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7} |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------------|")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.5f} |".format(chi_max[itimestep,k])
    print(line)

print(" ---------------------------------------------------------|")
print(" Average quantum parameter")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7} |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------------|")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(minimal_iteration,maximal_iteration+period,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.5f} |".format(chi_ave[itimestep,k])
    print(line)
    for k,model in enumerate(radiation_list):
        # if 0, absolute error of 0.01
        if (chi_ave[itimestep,k] == 0):
            Validate("Average quantum parameter for the {} model at iteration {}".format(model,timestep),chi_ave[itimestep,k],0.01)
        # else validation with relative 150% error
        else:
            Validate("Average quantum parameter for the {} model at iteration {}".format(model,timestep),chi_ave[itimestep,k],chi_ave[itimestep,k]*1.)
