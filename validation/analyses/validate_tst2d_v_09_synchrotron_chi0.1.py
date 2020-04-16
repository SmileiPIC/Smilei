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
import matplotlib.pyplot as plt
import matplotlib as mpl
import happi

# ______________________________________________________________________________
# Useful functions

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
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]

# List of relativistic pushers
radiation_list = ["CLL","Niel","MC"]

# ______________________________________________________________________________
# Read scalar diagnostics

print("")
print(" 1) Analyze of scalar diags")
print("")

ukin = {}
urad = {}
utot = {}
ntot = {}

# We load successively the particle track associated
# to each radiation algorithm
for radiation in radiation_list:

    # Read scalar diagnostics
    ScalarUkinDiag = S.Scalar("Ukin_electron_"+radiation).get()
    ukin[radiation] = np.array(ScalarUkinDiag["data"])
    times = np.array(ScalarUkinDiag["times"])

    ScalarUradDiag = S.Scalar("Urad_electron_"+radiation).get()
    urad[radiation] = np.array(ScalarUradDiag["data"])

    utot[radiation] = ukin[radiation]+urad[radiation]

print(" -------------------------------------------------|")
print(" Diag scalars (Kinetic energy)                    |")
print(" iteration | CLL        | Niel       | MC         |")
print(" -------------------------------------------------|")

for it,time in enumerate(times[::5]):
    print(" {0:5d}     | {1:.4e} | {2:.4e} | {3:.4e} | ".format(it*500,ukin["CLL"][it*5],ukin["Niel"][it*5],ukin["MC"][it*5]))

print("")

for radiation in radiation_list:
    
    print(' Final kinetic energy for {}: {}'.format(radiation,ukin[radiation][-1]))
    print(' Final radiated energy for {}: {}'.format(radiation,urad[radiation][-1]))

threshols = {}
threshols["CLL"] = 0.05
threshols["Niel"] = 0.1
threshols["MC"] = 0.15

for radiation in radiation_list:

    # Validation of the kinetic energy
    for it,val in enumerate(ukin[radiation]):
        Validate("Kinetic energy evolution for {} at {}".format(radiation,it), val/utot[radiation][0], val/utot[radiation][0]*threshols[radiation])

    # Validation of the radiated energy
    for it,val in enumerate(urad[radiation]):
        Validate("Radiated energy evolution for {} at {}".format(radiation,it) , val/utot[radiation][0], val/utot[radiation][0]*threshols[radiation] )

    # Validation of the total energy
    Validate("Total energy error (max - min)/uref for {}".format(radiation),(utot[radiation].max() - utot[radiation].min())/utot[radiation][0], 1e-2)

print("")
print(" 2) Relative errors")
print("")

# ______________________________________________________________________________
# Comparison corrected Landau-Lifshitz and Niel model

urad_rel_err = abs(urad["Niel"] - urad["CLL"]) / urad["CLL"].max()
ukin_rel_err = abs(ukin["Niel"] - ukin["CLL"]) / ukin["CLL"][0]

print(" Comparison Laudau-Lifshitz/Niel methods")
print(" Maximum relative error kinetic energy: {}".format(ukin_rel_err.max()))
print(" Maximum relative error radiative energy: {}".format(urad_rel_err.max()))

# Validation difference between continuous and discontinuous methods
Validate("Relative error on the kinetic energy / ukin at t=0 (Niel/CLL) " , ukin_rel_err.max(), 0.011 )
Validate("Relative error on the radiative energy / urad max (Niel/CLL) " , urad_rel_err.max(), 0.011 )

# ______________________________________________________________________________
# Comparison corrected Landau-Lifshitz and MC model

urad_mc_rel_err = abs(urad["MC"] - urad["CLL"]) / urad["CLL"].max()
ukin_mc_rel_err = abs(ukin["MC"] - ukin["CLL"]) / ukin["CLL"][0]

print("")
print(' Comparison Laudau-Lifshitz/Monte Carlo methods')
print(' Maximum relative error kinetic energy: {}'.format(ukin_mc_rel_err.max()))
print(' Maximum relative error radiative energy: {}'.format(urad_mc_rel_err.max()))

# Validation difference between continuous and discontinuous methods
Validate("Relative error on the kinetic energy / ukin at t=0 (MC/CLL) " , ukin_mc_rel_err.max(), 0.013 )
Validate("Relative error on the radiative energy / urad max (MC/CLL) " , urad_mc_rel_err.max(), 0.013 )

# ______________________________________________________________________________
# Checking of the particle binning

maximal_iteration = 5500
period = 500
number_of_files = maximal_iteration/period

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
        
        # Local average chi from spatial distribution
        chi = weight_chi[weight>0] / weight[weight>0]
        # Maximal chi value from spatial distribution
        chi_max[itimestep,i] = chi.max()
        chi_ave[itimestep,i] = np.sum(chi)/sum(list(np.shape(chi)))

        # Maximal and average chi value from chi distribution
        if (itimestep > 0):
            chi_ave_from_dists[itimestep,i] = np.sum(chi_distr_axis * chi_distr_data * bin_size) / np.sum(chi_distr_data * bin_size)
        imax = np.argmax(chi_distr_data)
        chi_max_from_dists[itimestep,i] = chi_distr_axis[imax]
        
        # Average kinetic energy from energy distribution
        ekin_from_dists[itimestep,i] = dx*dy*np.sum(ekin_distr_axis * ekin_distr_data * bin_size)
        if (itimestep > 0):
            ekin_ave_from_dists[itimestep,i] = ekin_from_dists[itimestep,i] / np.sum(ekin_distr_data * bin_size)

print("")
print(" 3) Analyze of chi x-y particle binning diagnostics")
print("")

print(" ---------------------------------------------------")
print(" Maximal quantum parameter")
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
    
# Validation with 50% error
# The maximal quantum parameter can vary importantly
# for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    # for k,model in enumerate(radiation_list):
    #    Validate("Maximal quantum parameter for the {} model at iteration {}".format(model,timestep),chi_max[itimestep,k],chi_max[itimestep,k]*0.5)

print(" ---------------------------------------------------")
print(" Average quantum parameter")
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

print("")
print(" 4) Analyze of chi using the chi distribution")
print("")

print(" ---------------------------------------------------------")
print(" Average quantum parameter from chi distribution         |")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<9}  |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.4e} |".format(chi_ave_from_dists[itimestep,k])
    print(line)
    # Validation with 10% error
    for k,model in enumerate(radiation_list):
        Validate("Average quantum parameter for the {} model at iteration {}".format(model,timestep),chi_ave_from_dists[itimestep,k],chi_ave_from_dists[itimestep,k]*0.1)

print("")
print(" 5) Analyze of gamma using the gamma distribution")
print("")

print(" ---------------------------------------------------------")
print(" Total kinetic energy from energy distribution           |")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<9}  |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.4e} |".format(ekin_from_dists[itimestep,k])
    print(line)
    
# Validation
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    for k,model in enumerate(radiation_list):
        Validate("Total kinetic energy for the {} model at iteration {}".format(model,timestep),ekin_from_dists[itimestep,k],ekin_from_dists[itimestep,k]*0.1)

print(" ---------------------------------------------------------")
print(" Average kinetic energy from energy distribution         |")
line = "                  |"
for k,model in enumerate(radiation_list):
    line += " {0:<9}  |".format(radiation_list[k])
print(line)
print(" ---------------------------------------------------------")
# Loop over the timesteps
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    line = " Iteration {0:5d}  |".format(timestep)
    for k,model in enumerate(radiation_list):
        line += " {0:.4e} |".format(ekin_ave_from_dists[itimestep,k])
    print(line)
    
# Validation
for itimestep,timestep in enumerate(range(0,maximal_iteration,period)):
    for k,model in enumerate(radiation_list):
        Validate("Average kinetic energy for the {} model at iteration {}".format(model,timestep),ekin_ave_from_dists[itimestep,k],ekin_ave_from_dists[itimestep,k]*0.1)

# ______________________________________________________________________________
# Figures

# fig, ax = plt.subplots(figsize = (8, 6))
#
# ax.plot(times,ukin["CLL"],color='C0',label='CLL',ls='-')
# ax.plot(times,urad["CLL"],color='C0',label='',ls='--')
#
# ax.plot(times,ukin["Niel"],color='C1',label='Niel',ls='-')
# ax.plot(times,urad["Niel"],color='C1',label='',ls='--')
#
# ax.plot(times,ukin["MC"],color='C2',label='MC',ls='-')
# ax.plot(times,urad["MC"],color='C2',label='',ls='--')
#
# ax.legend(loc='best')
#
# fig.savefig('tst2d_09_scalar_energy_balance.png')
