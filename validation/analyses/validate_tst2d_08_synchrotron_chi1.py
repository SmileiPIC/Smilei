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
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]

# List of relativistic pushers
radiation_list = ["CLL","Niel","MC"]

ukin_dict = {}
urad_dict = {}

# ______________________________________________________________________________
# We load successively particle scalar diagnostics associated
# to each radiation algorithm
for radiation in radiation_list:

    # Read scalar diagnostics
    ScalarUkinDiag = S.Scalar("Ukin_electron_"+radiation).get()
    ukin = np.array(ScalarUkinDiag["data"])
    times = np.array(ScalarUkinDiag["times"])

    ScalarUradDiag = S.Scalar("Urad_electron_"+radiation).get()
    urad = np.array(ScalarUradDiag["data"])

    utot = ukin+urad

    print(' Final kinetic energy for {}: {}'.format(radiation,ukin[-1]))

    # Validation of the kinetic energy
    Validate("Kinetic energy evolution for {}".format(radiation), ukin/ukin[0], 5.e-2 )

    # Validation of the radiated energy
    Validate("Radiated energy evolution for {}".format(radiation) , urad/ukin[0], 5.e-2 )

    # Validation of the total energy
    Validate("Total energy error (max - min)/uref for {}".format(radiation) ,
           (utot.max() - utot.min())/utot[0], 1e-2)

    ukin_dict[radiation] = ukin
    urad_dict[radiation] = urad

# ______________________________________________________________________________
# Comparison Corrected Landau-Lifshitz and the Monte-Carlo method

for k,model in enumerate(radiation_list[1:]):

    urad_rel_err = abs(urad_dict[model] - urad_dict["CLL"]) / urad_dict["CLL"].max()
    ukin_rel_err = abs(ukin_dict[model] - ukin_dict["CLL"]) / ukin_dict["CLL"][0]

    print('')
    print(' Comparison Corrected Landau-Lifshitz / {} methods'.format(model))
    print(' Maximum relative error kinetic energy: {}'.format(ukin_rel_err.max()))
    print(' Maximum relative error radiative energy: {}'.format(urad_rel_err.max()))

    # Validation difference between Corrected Landau-Lifshitz and Monte-Carlo methods
    Validate("Relative error on the kinetic energy / ukin at t=0 ({})".format(model) , ukin_rel_err.max(), 0.02 )
    Validate("Relative error on the radiative energy / urad max ({})".format(model) , urad_rel_err.max(), 0.02 )

# ______________________________________________________________________________
# Checking of the particle binning

print("")
print(" Checking of the particle binning diagnostics")

maximal_iteration = 5500
period = 500
number_of_files = maximal_iteration//period

chi_max = np.zeros([number_of_files,len(radiation_list)])
chi_ave = np.zeros([number_of_files,len(radiation_list)])

chi_ave_from_dists = np.zeros([number_of_files,len(radiation_list)])
chi_max_from_dists = np.zeros([number_of_files,len(radiation_list)])

ekin_ave_from_dists = np.zeros([number_of_files,len(radiation_list)])
ekin_max_from_dists = np.zeros([number_of_files,len(radiation_list)])
ekin_from_dists = np.zeros([number_of_files,len(radiation_list)])

# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):

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
        chi_distr_axis = np.array(chi_dist["chi"])
        chi_distr_data = np.array(chi_dist["data"][0])
        log10_chi_distr_axis = np.log10(chi_distr_axis)
        delta = log10_chi_distr_axis[1] - log10_chi_distr_axis[0]
        bin_size =  np.power(10.,log10_chi_distr_axis + 0.5*delta) - np.power(10.,log10_chi_distr_axis - 0.5*delta)
        # Energy distribution
        ekin_dist = S.ParticleBinning(diagNumber=i+3*len(radiation_list),timesteps=timestep).get()
        ekin_distr_axis = np.array(ekin_dist["ekin"])
        ekin_distr_data = np.array(ekin_dist["data"][0])
        log10_ekin_distr_axis = np.log10(ekin_distr_axis)
        delta = log10_ekin_distr_axis[1] - log10_ekin_distr_axis[0]
        bin_size =  np.power(10.,log10_ekin_distr_axis + 0.5*delta) - np.power(10.,log10_ekin_distr_axis - 0.5*delta)
        # Local average chi
        chi = weight_chi[weight>0] / weight[weight>0]
        # Maximal and average chi value from spatial distribution
        chi_max[itimestep,i] = chi.max()
        chi_ave[itimestep,i] = np.sum(chi)/sum(list(np.shape(chi)))
        # Maximal and average chi value from chi distribution
        if (itimestep > 0):
            chi_ave_from_dists[itimestep,i] = np.sum(chi_distr_axis * chi_distr_data * bin_size) / np.sum(chi_distr_data * bin_size)
        imax = np.argmax(chi_distr_data)
        chi_max_from_dists[itimestep,i] = chi_distr_axis[imax]
        # Average kinetic energy
        if (itimestep > 0):
            ekin_ave_from_dists[itimestep,i] = dx*dy*np.sum(ekin_distr_axis * ekin_distr_data * bin_size) / np.sum(ekin_distr_data * bin_size)
        ekin_from_dists[itimestep,i] = dx*dy*np.sum(ekin_distr_axis * ekin_distr_data * bin_size)

print(" ---------------------------------------------------")
print(" Maximal quantum parameter from 2D spatial particle binning")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7} |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.5f} |".format(chi_max[itimestep,k])
    print(line)
    # Validation with 90% error
    # The maximal quantum parameter can vary a lot
    # for k,model in enumerate(radiation_list):
    #    Validate("Maximal quantum parameter for the {} model at iteration {}".format(model,timestep),chi_max[itimestep,k],chi_max[itimestep,k]*0.9)

print(" ---------------------------------------------------")
print(" Maximal quantum parameter from chi distribution")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7} |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.5f} |".format(chi_max_from_dists[itimestep,k])
    print(line)

print(" ---------------------------------------------------")
print(" Average quantum parameter from the 2D spatial particle binning")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7} |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.5f} |".format(chi_ave[itimestep,k])
    print(line)
    # Validation with 10% error
    for k,model in enumerate(radiation_list):
        Validate("Average quantum parameter for the {} model at iteration {}".format(model,timestep),chi_ave[itimestep,k],chi_ave[itimestep,k]*0.1)

print(" ---------------------------------------------------")
print(" Average quantum parameter from chi distribution")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7} |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.5f} |".format(chi_ave_from_dists[itimestep,k])
    print(line)
    # Validation with 10% error
    for k,model in enumerate(radiation_list):
        Validate("Average quantum parameter for the {} model at iteration {}".format(model,timestep),chi_ave_from_dists[itimestep,k],chi_ave_from_dists[itimestep,k]*0.1)

print(" ---------------------------------------------------")
print(" Total kinetic energy from energy distribution")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7}  |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.5f} |".format(ekin_from_dists[itimestep,k])
    print(line)
    # Validation with 10% error
    for k,model in enumerate(radiation_list):
        Validate("Total kinetic energy for the {} model at iteration {}".format(model,timestep),ekin_from_dists[itimestep,k],ekin_ave_from_dists[itimestep,k]*0.1)

print(" ---------------------------------------------------")
print(" Average kinetic energy from energy distribution")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<7}  |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.5f} |".format(ekin_ave_from_dists[itimestep,k])
    print(line)
    # Validation with 10% error
    for k,model in enumerate(radiation_list):
        Validate("Average kinetic energy for the {} model at iteration {}".format(model,timestep),ekin_ave_from_dists[itimestep,k],ekin_ave_from_dists[itimestep,k]*0.1)
