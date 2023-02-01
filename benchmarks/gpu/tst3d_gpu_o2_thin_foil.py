# ----------------------------------------------------------------------------------------
#                                       SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Laser-thin foil interaction at extremen intensity with radiations
# - plane wave
# From: https://<>/smilei/benchmarking/-/raw/master/benchmarks/thin_foil_3d_o2/template.py
# ----------------------------------------------------------------------------------------

import math
import datetime
#import numpy as np
import os

# GENERAL PHYSICAL PARAMETERS
# ----------------------------------------

alpha = 1./137.
HBar = 1.0545718e-34                              # (*m^2kg/s*)
c = 299792458
electron_mass = 9.10938356e-31                    # (*kg*)
electron_charge = 1.60217662e-19                  # (*coulombs*)
lambdar = 1e-6                                    # Normalization wavelength
wr = 2*math.pi*c/lambdar                          # Normalization frequency
Schwinger_E_field= 1.3E18                         # Schwinger electric field
Enorm = electron_mass*wr*c/electron_charge        # Normalization electric field at lambda = 1e-6

# TARGET PROPERTIES
# ----------------------------------------

charge_density = 711.90  # normalized target charge density

Te = 0.0001 # Electron temperature
Ti = 0.00001 # Electron temperature

particles_per_cell = 8

# RESOLUTION, BOX SIZE
# ----------------------------------------

l0 = 2.0*math.pi                                  # laser wavelength

cells_per_wavelength = 32.                      # nb of cells per wavelength

grid_length = [10 * l0, 2*l0, 2*l0]             # size of the domain

cell_length = [l0/cells_per_wavelength, l0/cells_per_wavelength, l0/cells_per_wavelength]

x_start_preplasma = 3 * l0     # x position where the preplasma starts
x_start_foil = 4 * l0          # x position where the foil starts
foil_width = 3 * l0     # width of the foil in the x direction

number_of_patches = [4, 1, 1]

boundary_conditions = [["reflective", "reflective"],["periodic", "periodic"],["periodic", "periodic"]]

print(number_of_patches)

# LASER
# ----------------------------------------

t0 = l0

a0 = 10.
start = 0                               # Laser start
fwhm = 10*t0                            # Gaussian time fwhm
duration = 30*t0                        # Laser duration
center = duration*0.5                   # Laser profile center
order = 4                               # Gaussian order

# TIME
# ---------------------------------------

dt  = 0.9/math.sqrt(1./(cell_length[0]**2) + 1./(cell_length[1]**2) + 1./(cell_length[2]**2))       # timestep given by the CFL
#simulation_time = duration + x_start_foil                       # simulation duration
simulation_time = duration + x_start_foil                       # simulation duration

# OTHER PARAMETERS
# ----------------------------------------

pusher = "boris"                                    # type of pusher
radiation_model = 'cll'

# DIAGS
# ----------------------------------------

print_every = 50
diag_every = 2*int(t0 / dt)

print("diag_every: {}".format(diag_every))

random_seed = 0xDEADBEEF

gpu_computing = True
vectorization = "off"

# FUNCTIONS
# ----------------------------------------

# Density profile for inital location of the particles
def charge_density_profile(x,y,z):
    # preplasma ramp
    if (x_start_preplasma < x <= x_start_foil):
        return charge_density*(x - x_start_preplasma) / (x_start_foil - x_start_preplasma)
    # main target
    elif (x_start_foil < x <= x_start_foil + foil_width):
        return charge_density
    else:
        return 0

# NAMELIST
# ----------------------------------------

Main(
    geometry = "3Dcartesian",
    interpolation_order = 2,
    cell_length = cell_length,
    grid_length  = grid_length,
    number_of_patches = number_of_patches,
    #time_fields_frozen = simulation_time,
    timestep = dt,
    simulation_time = simulation_time,
    EM_boundary_conditions =
       [["silver-muller", "silver-muller"],
        ["periodic", "periodic"],
        ["periodic", "periodic"]],
    random_seed = random_seed,
    solve_poisson = False,
    gpu_computing=gpu_computing,
    reference_angular_frequency_SI = wr,
    print_every = print_every,

)

LaserGaussian3D(
    box_side         = "xmin",
    a0              = a0,
    omega           = 1.,
    focus           = [x_start_foil, grid_length[1]*0.5, grid_length[2]*0.5],
    waist           = 1e9,
    incidence_angle = [0., 0.],
    polarization_phi = 0.,
    ellipticity     = 0.,
    time_envelope   = tgaussian(start=start,duration=duration,fwhm=fwhm,center=center,order=order)
)

Vectorization(
    mode=vectorization,
)

# ----------------------------------------------------------------------------------------

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = charge_density_profile,
    mean_velocity = [0.],
    temperature = [Te],
    pusher = pusher,
    boundary_conditions = boundary_conditions,
#    radiation_model = radiation_model,
    time_frozen = x_start_preplasma-cell_length[0],
#    radiation_photon_species = "photon",
#    radiation_photon_gamma_threshold = 0,
)

Species(
    name = "ion",
    position_initialization = "electron",
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 49187,
    charge = 13.0,
    charge_density = charge_density_profile,
    mean_velocity = [0.],
    temperature = [Ti],
    pusher = pusher,
    boundary_conditions = boundary_conditions,
#    radiation_model = radiation_model,
    time_frozen = x_start_preplasma-cell_length[0],
#    radiation_photon_species = "photon",
#    radiation_photon_gamma_threshold = 0,
)

# ----------------------------------------------------------------------------------------
# Global parameters for the radiation reaction models
# RadiationReaction(
#     minimum_chi_discontinuous = 1e-3,
#     minimum_chi_continuous = 1e-4,
# )

# ----------------------------------------------------------------------------------------
# Scalar diagnostics
DiagScalar(
    every = 5,
)

# ----------------------------------------------------------------------------------------
# Particle Binning

# # Energy distributions of synchrotron photons
#DiagParticleBinning(
#    deposited_quantity = "weight",
#    every = 1,
#    species =["photon"],
#    axes = [["gamma", 0., gamma, 256]]
#)

#DiagParticleBinning(
#    deposited_quantity = "weight",
#    every = 1,
#    species =["photon"],
#    axes = [["gamma", 1e-3, gamma, 256, "logscale"]]
#)

#DiagParticleBinning(
#    deposited_quantity = "weight",
#    every = 1,
#    species =["electron"],
#    axes = [["chi", 1e-3, 1., 256, "logscale"]]
#)

DiagFields(
    every = diag_every,
    fields = ["Ey", "Ez", "By", "Bz"],
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., grid_length[0], 256],
        ["y", 0., grid_length[1], 64],
        ["z", 0., grid_length[2], 64]
    ]
)

DiagProbe(
    #name = "my_probe",
    every    = diag_every,
    origin   = [0, 0, 0],
    corners  = [
        [grid_length[0],0., 0.],
        [0.,grid_length[1], 0.],
        [0.,0., grid_length[2]],
    ],
    number   = [256, 64, 64],
    fields   = ["Ex", "Ey", "Ez"]
)
