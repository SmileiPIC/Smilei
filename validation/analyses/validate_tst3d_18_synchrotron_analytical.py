# ______________________________________________________________________________
#
# Validation script for the magnetic shower with comparison to theory
#
# In this tests case, the MC model is compared to analytics after 1 iteration.
#
# Validation:
# - Monte-Carlo radiation loss
# - Species scalar diagnostics
# - External fields
# - Particle binning with the quantum parameter
# ______________________________________________________________________________

import os, re, numpy as np, h5py
# import matplotlib.pyplot as plt
# import matplotlib as mpl
import happi
from scipy.special import kv
from scipy.integrate import quad
import scipy.integrate as integrate

# ______________________________________________________________________________
# Useful functions

def K1By3(y):
  return kv(1./3., y)

def d2N_dtdgamma(photon_chi, particle_chi, particle_gamma):
  y = (photon_chi)/(3.*particle_chi*(particle_chi - photon_chi))
  term1 = 1./particle_gamma**2
  term2 = ( 2 + 3*photon_chi*y) * kv(2./3., 2*y)
  term3 = - integrate.quad(K1By3, 2*y, float('inf'))[0]
  result = term1 * (term2 + term3)
  return result

# ______________________________________________________________________________
# opening of the case

S = happi.Open(["./restart*"], verbose=False)
dx    = S.namelist.Main.cell_length[0]
dy    = S.namelist.Main.cell_length[1]
dz    = S.namelist.Main.cell_length[2]
dt    = S.namelist.Main.timestep

electron_chi   = S.namelist.chi  # the value of chi and gamma are taken from the simualtion
electron_gamma = S.namelist.gamma  # same as electron energy in mc2

alpha          = S.namelist.alpha
electron_mass  = S.namelist.electron_mass
c              = S.namelist.c
HBar           = S.namelist.HBar
wr             = S.namelist.wr

coefficient = 2.*alpha*electron_mass*c*c*np.sqrt(3) / (6.*np.pi * HBar)

print(" Coefficient: {}".format(coefficient))

Lx = S.namelist.Main.grid_length[0]
Ly = S.namelist.Main.grid_length[1]
Lz = S.namelist.Main.grid_length[2]

volume = Lx*Ly*Lz

simulation_time = S.namelist.Main.simulation_time / wr
iteration = int(S.namelist.Main.simulation_time / dt)

print("")
print(" 1) Analyze of scalar diags")
print("")

ukin_photon = np.array(S.Scalar("Ukin_photon").getData())
dens_photon = np.array(S.Scalar("Dens_photon").getData())
ukin_electron = np.array(S.Scalar("Ukin_electron").getData())
dens_electron = np.array(S.Scalar("Dens_electron").getData())
timesteps = np.array(S.Scalar("Ukin_photon").get()["times"])

print(" Final electron density: {}".format(dens_electron[-1]))
print(" Final photon density: {}".format(dens_photon[-1]))

Validate("Final electron density",dens_electron[-1], dens_electron[-1]*1e-2)
Validate("Final photon density",dens_photon[-1], dens_photon[-1]*3e-2)

print(" Final electron kinetic energy: {}".format(ukin_electron[-1]))
print(" Final photon kinetic energy: {}".format(ukin_photon[-1]))

Validate("Final electron kinetic energy",ukin_electron[-1], ukin_electron[-1]*1e-2)
Validate("Final photon kinetic energy",ukin_photon[-1], ukin_photon[-1]*3e-2)

print("")
print(" 2) Particle binning: photon spectrum (linear scale)")
print("")

photon_spectrum = S.ParticleBinning(0).get( iteration )
gamma_axis =  photon_spectrum['gamma']
data = photon_spectrum['data'][0]
delta = gamma_axis[1] - gamma_axis[0]
total_weight = np.sum(data)*volume*delta
total_energy = np.sum(data*gamma_axis)*volume*delta

rel_total_weight = np.abs(total_weight - dens_photon[-1])/dens_photon[-1]
rel_total_energy = np.abs(total_energy - ukin_photon[-1])/ukin_photon[-1]

print(" Photon density: {} ({})".format(total_weight,rel_total_weight))
print(" Photon kinetic energy: {} ({})".format(total_energy,rel_total_energy))

Validate("Diff photon density (lin spectrum vs scalar)",rel_total_weight, 1e-5)
Validate("Diff photon ukin (lin spectrum vs scalar)",rel_total_energy, 1e-2)


print("")
print(" 3) Particle binning: photon spectrum (logscale)")
print("")

photon_spectrum = S.ParticleBinning(1).get( iteration )
gamma_axis =  photon_spectrum['gamma']
data = photon_spectrum['data'][0]*volume
log10_gamma_axis = np.log10(gamma_axis)
delta = log10_gamma_axis[1] - log10_gamma_axis[0]
bin_size =  np.power(10.,log10_gamma_axis + 0.5*delta) - np.power(10.,log10_gamma_axis - 0.5*delta)
total_weight = np.sum(bin_size * data)
total_energy = np.sum(bin_size * data * gamma_axis)

rel_total_weight = np.abs(total_weight - dens_photon[-1])/dens_photon[-1]
rel_total_energy = np.abs(total_energy - ukin_photon[-1])/ukin_photon[-1]

print(" Photon density: {} ({})".format(total_weight,rel_total_weight))
print(" Photon kinetic energy: {} ({})".format(total_energy,rel_total_energy))

Validate("Diff photon density (log spectrum vs scalar)",rel_total_weight, 1e-2)
Validate("Diff photon ukin (log spectrum vs scalar)",rel_total_energy, 1e-2)

print("")
print(" 4) Analytical estimate")
print("")

analytical_spectrum = np.zeros(len(gamma_axis))

for k,photon_gamma in enumerate(gamma_axis):
    photon_chi = electron_chi * photon_gamma / electron_gamma
    analytical_spectrum[k] = (d2N_dtdgamma(photon_chi, electron_chi, electron_gamma)) # gamma is the normalised electron energy

analytical_spectrum = analytical_spectrum * simulation_time * dens_electron[iteration] * coefficient

ukin_theory = np.sum(analytical_spectrum*gamma_axis*bin_size)

print(" Kinetic energy from theroy: {}".format(ukin_theory))

Validate("Kinetic energy from theroy",ukin_theory, ukin_theory*0.1)

diff = np.abs((data - analytical_spectrum) / analytical_spectrum)

diff = diff[0:253]

sum_diff = np.sum(diff)
mean_diff = np.mean(diff)
max_diff = np.max(diff)
min_diff = np.min(diff)

print(" Sum diff: {}".format(sum_diff))
print(" Mean diff: {}".format(mean_diff))
print(" Max diff: {} ({})".format(max_diff,np.argmax(diff)))
print(" Min diff: {} ({})".format(min_diff,np.argmin(diff)))

Validate("Sum diff",sum_diff, sum_diff*0.2)
Validate("Mean diff",mean_diff, mean_diff*0.2)

# ______________________________________________________________________________
# Figures

# fig, ax = plt.subplots(figsize = (8, 6))
# ax.plot(gamma_axis, data * gamma_axis, '-', color = "C0", label='spectrum', linewidth=2)
# ax.plot(gamma_axis, analytical_spectrum * gamma_axis, '-', color = "C1", label='analytical', linewidth=2)
# ax.set_xscale('log')
# ax.legend()
# fig.tight_layout()
