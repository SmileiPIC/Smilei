########################################################################################################################
# Setup:                                                                                                               #
# - Incoherent (synchrotron) radiation emission by ultra-relativistic electron                                         #
# - relativistic electrons radiating in a constant magnetic field                                                      #
# - EM fields are frozen: 4 electron species are considered: not radiation, LL, cLL, MC                                #
# - reference frequency is the gyrofrequency                                                                           #
# Main test:                                                                                                           #
# - RadiationSpectrum diagnostics                                                                                      #
# Author:                                                                                                              #
# mickael.grech@polytechnique.edu                                                                                      #
########################################################################################################################

import numpy as np

# quick access databases
chi0    = 0.10
g0      = 1.e3

# Physical constants in SI units
c_SI    = 299792458.
me_SI   = 9.10938356e-31
e_SI    = 1.60217662e-19
hbar_SI = 1.054571800e-34

# Simulation box properties
dx      = 1./128.
Lx      = 1.
dt      = 0.95*dx
Tsim    = 2.*np.pi

# Electron bunch properties
v0      = np.sqrt(1.-1./g0**2)
n0      = 1.
l0      = 1.
nppc    = 32

def n_(x):
    if (x<l0):
        return n0
    else:
        return 0.

# External magnetic DiagField
B0      = g0

# Estimate maximum photon energy
photon_energy_max = 10.*g0*chi0

# Compute reference angular frequency
Er_ov_Es = chi0 / g0**2
wr       = me_SI*c_SI**2/hbar_SI * Er_ov_Es

# SMILEI PARAMETERS ####################################################################################################

### MAIN
Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    cell_length = [dx],
    grid_length  = [Lx],
    number_of_patches = [16],
    timestep = dt,
    simulation_time = Tsim,
    EM_boundary_conditions = [['periodic']],
    random_seed = smilei_mpi_rank,
    reference_angular_frequency_SI = wr,
    solve_poisson = False,
    time_fields_frozen = 2.*Tsim,
    print_every = int(Tsim/dt/20.)
)

# PrescribedField(
#     field   = 'Bz',
#     profile = lambda x,t: 1000.
# )

ExternalField(
    field   = 'Bz',
    profile = B0
)


Species(
    name = "electron_noRR",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = n_,
    mean_velocity = [v0,0.,0.],
    boundary_conditions = [["periodic"]],
    radiation_model = "diagradiationspectrum"
)

Species(
    name = "electron_LL",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = n_,
    mean_velocity = [v0,0.,0.],
    boundary_conditions = [["periodic"]],
    radiation_model = "LL"
)

Species(
    name = "electron_cLL",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = n_,
    mean_velocity = [v0,0.,0.],
    boundary_conditions = [["periodic"]],
    radiation_model = "cLL"
)

Species(
    name = "electron_FP",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = n_,
    mean_velocity = [v0,0.,0.],
    boundary_conditions = [["periodic"]],
    radiation_model = "Niel"
)

Species(
    name = "electron_MC",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = n_,
    mean_velocity = [v0,0.,0.],
    boundary_conditions = [["periodic"]],
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 0
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    mass = 0.,
    charge = 0.,
    number_density = 0,
    mean_velocity = [0.],
    boundary_conditions = [["periodic"]],
)

RadiationReaction(
   Niel_computation_method = "fit5",
   # Radiation parameters
   minimum_chi_continuous = 1e-6,
   minimum_chi_discontinuous = 1e-4,
)

### Diagnostics
globalEvery = int(Tsim/dt)

DiagFields(
    every = 50
)

### (vx-vy) phase-space
DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_noRR"],
    axes = [
        ["vx", -1, 1, 50],["vy", -1, 1, 50]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_LL"],
    axes = [
        ["vx", -1, 1, 50],["vy", -1, 1, 50]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_cLL"],
    axes = [
        ["vx", -1, 1, 50],["vy", -1, 1, 50]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_FP"],
    axes = [
        ["vx", -1, 1, 50],["vy", -1, 1, 50]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_MC"],
    axes = [
        ["vx", -1, 1, 50],["vy", -1, 1, 50]
    ]
)

### distribution in chi
DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_noRR"],
    axes = [
        ["chi", 0, 4*chi0, 100]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_LL"],
    axes = [
        ["chi", 0, 4*chi0, 100]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_cLL"],
    axes = [
        ["chi", 0, 4*chi0, 100]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_FP"],
    axes = [
        ["chi", 0, 4*chi0, 100]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["electron_MC"],
    axes = [
        ["chi", 0, 4*chi0, 100]
    ]
)

### photon energy
DiagRadiationSpectrum(
    every = 1,
    species = ["electron_noRR"],
    photon_energy_axis = [photon_energy_max/1.e6,photon_energy_max, 400, 'logscale'],
    axes = []
)

DiagRadiationSpectrum(
    every = 1,
    species = ["electron_LL"],
    photon_energy_axis = [photon_energy_max/1.e6,photon_energy_max, 400, 'logscale'],
    axes = []
)

DiagRadiationSpectrum(
    every = 1,
    species = ["electron_cLL"],
    photon_energy_axis = [photon_energy_max/1.e6,photon_energy_max, 400, 'logscale'],
    axes = []
)

DiagRadiationSpectrum(
    every = 1,
    species = ["electron_FP"],
    photon_energy_axis = [photon_energy_max/1.e6,photon_energy_max, 400, 'logscale'],
    axes = []
)

def depose(particles):
    return particles.weight*np.sqrt(particles.px**2+particles.py**2+particles.pz**2)

DiagParticleBinning(
    deposited_quantity = depose,
    every = 1,
    species = ["photon"],
    axes = [
        ["gamma", photon_energy_max/1.e6, photon_energy_max, 400, 'logscale']
    ]
)

DiagScalar(every=1)
