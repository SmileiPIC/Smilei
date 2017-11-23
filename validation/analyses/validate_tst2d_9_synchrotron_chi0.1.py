# ______________________________________________________________________________
#
# Validation script for the synchotron case
#
# In this tests case, an electron bunch is initialized per radiation
# loss models at the same positions. The magnetic field and the initial energy
# is computed so that the initial quantum parameter is equal to 1.
#
# Validation:
# - Discontinuous radiation loss
# - Continuous radiation loss
# ______________________________________________________________________________


import os, re, numpy as np, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



# List of relativistic pushers
radiation_list = ["Niel","Landau_Lifshitz"]

ukin_dict = {}
urad_dict = {}

# We load successively the particle track associated
# to each radiation algorithm
for radiation in radiation_list:

  # Read scalar diagnostics
  ScalarUkinDiag = S.Scalar("Ukin_electron_"+radiation).get()
  ukin = np.array(ScalarUkinDiag["data"])
  times = np.array(ScalarUkinDiag["times"])

  ScalarUradDiag = S.Scalar("Urad_electron_"+radiation).get()
  urad = np.array(ScalarUradDiag["data"])

  utot = ukin+urad

  # Validation of the kinetic energy
  Validate("Kinetic energy evolution: ", ukin/ukin[0], 5.e-2 )

  # Validation of the radiated energy
  Validate("Radiated energy evolution: " , urad/ukin[0], 5.e-2 )

  # Validation of the total energy
  Validate("Total energy error (max - min)/uref: " ,
           (utot.max() - utot.min())/utot[0], 1e-2)

  ukin_dict[radiation] = ukin
  urad_dict[radiation] = urad

# ______________________________________________________________________________
# Comparision continuous and Niel model

urad_rel_err = abs(urad_dict["Niel"] - urad_dict["Landau_Lifshitz"]) / urad_dict["Landau_Lifshitz"].max()
ukin_rel_err = abs(ukin_dict["Niel"] - ukin_dict["Landau_Lifshitz"]) / ukin_dict["Landau_Lifshitz"][0]

print ' Comparision Laudau-Lifshitz/Niel methods'
print ' Maximum relative error kinetic energy',ukin_rel_err.max()
print ' Maximum relative error radiative energy',urad_rel_err.max()

# Validation difference between continuous and discontinuous methods
Validate("Relative error on the kinetic energy / ukin at t=0: " , ukin_rel_err.max(), 0.01 )
Validate("Relative error on the radiative energy / urad max " , urad_rel_err.max(), 0.01 )

# ______________________________________________________________________________
# Checking of the particle binning

print("")
print(" Checking of the particle binning diagnostics")

chi_max = np.zeros([10,len(radiation_list)])
chi_ave = np.zeros([10,len(radiation_list)])

# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,5000,500)):

    # Loop over the species/radiation models
    for i,radiation in enumerate(radiation_list):
        # Weight
        weight_diag = S.ParticleBinning(diagNumber=i,timesteps=timestep).get()
        weight = np.array(weight_diag["data"][0])
        # Weight x chi
        weight_chi_diag = S.ParticleBinning(diagNumber=i+len(radiation_list),timesteps=timestep).get()
        weight_chi = np.array(weight_chi_diag["data"][0])
        # Chi distribution
        chi_dist = S.ParticleBinning(diagNumber=i+2*len(radiation_list),timesteps=timestep).get()
        # Local average chi
        chi = weight_chi[weight>0] / weight[weight>0]
        # Maximal chi value
        chi_max[itimestep,i] = chi.max()
        chi_ave[itimestep,i] = np.sum(chi)/sum(list(np.shape(chi)))

print(" ---------------------------------------------------")
print(" Maximal quantum parameter")
print("                  | {0:<7} | {1:<7} |".format(radiation_list[0],radiation_list[1]))
print(" ---------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,5000,500)):
	print(" Iteration {0:5d}  | {1:.5f} | {2:.5f} |".format(timestep,chi_max[itimestep,0],chi_max[itimestep,1]))

print(" ---------------------------------------------------")
print(" Average quantum parameter")
print("                  | {0:<7} | {1:<7} |".format(radiation_list[0],radiation_list[1]))
print(" ---------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,5000,500)):
	print(" Iteration {0:5d}  | {1:.5f} | {2:.5f} |".format(timestep,chi_ave[itimestep,0],chi_ave[itimestep,1]))

Validate("Maximal quantum parameter",chi_max,0.05)
Validate("Average quantum parameter",chi_ave,0.05)


