# ----------------------------------------------------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
#
# Application of the particle merging on a thermal plasma
#
# ----------------------------------------------------------------------------------------

import math as m
import numpy as np

# Physical parameters
c = 299792458
lambdar = 1e-6
wr = 2*m.pi*c/lambdar

# electron & ion temperature in me c^2
Te   = 100./511.
Ti   = 10./511.

# Density
n0  = 1.
# density profile
def n0_(x):
   return n0

# Plasma length
lambdap = 2*m.pi/n0
# Debye length in units of c/\omega_{pe}
Lde = m.sqrt(Te)

# Cell length
cell_length = [Lde / 16.]
# Number of patches
number_of_patches = [16]
# grid length
grid_length = [16*Lde]
# time step
dt  = 0.95 * cell_length[0]
# Simulation time
simulation_time  = 50*dt
# Initial number of particles sper cell
particles_per_cell = 20480
# How the particles are initialized
position_initialization = 'random'

# Discretization in the momentum space
merge_momentum_cell_size = [8,8,8]

Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    timestep = dt,
    simulation_time = simulation_time,
    cell_length  = cell_length,
    grid_length = grid_length,
    number_of_patches = number_of_patches,
    EM_boundary_conditions = [ ["periodic","periodic"] ],
    random_seed = 0,
)

Species(
    name = "proton_cartesian",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    charge_density = n0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [Ti],
    pusher = "boris",
    boundary_conditions = [["periodic", "periodic"],],
    # Merging parameters
    merging_method = "vranic_cartesian",
    merge_every = 1,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = merge_momentum_cell_size,
    merge_discretization_scale = "linear",
    merge_accumulation_correction = True,
)
Species(
    name = "electron_cartesian",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [Te],
    pusher = "boris",
    boundary_conditions = [["periodic", "periodic"],],
    # Merging parameters
    merging_method = "vranic_cartesian",
    merge_every = 1,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = merge_momentum_cell_size,
    merge_discretization_scale = "linear",
    merge_accumulation_correction = True,
)

Species(
    name = "proton_spherical_lin",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    charge_density = n0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [Ti],
    pusher = "boris",
    boundary_conditions = [["periodic", "periodic"],],
    # Merging parameters
    merging_method = "vranic_spherical",
    merge_every = 1,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = merge_momentum_cell_size,
    merge_discretization_scale = "linear",
)
Species(
    name = "electron_spherical_lin",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [Te],
    pusher = "boris",
    boundary_conditions = [["periodic", "periodic"],],
    # Merging parameters
    merging_method = "vranic_spherical",
    merge_every = 1,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = merge_momentum_cell_size,
    merge_discretization_scale = "linear",
)

Species(
    name = "proton_spherical_log",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    charge_density = n0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [Ti],
    pusher = "boris",
    boundary_conditions = [["periodic", "periodic"],],
    # Merging parameters
    merging_method = "vranic_spherical",
    merge_every = 1,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = merge_momentum_cell_size,
    merge_discretization_scale = "log",
)
Species(
    name = "electron_spherical_log",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [Te],
    pusher = "boris",
    boundary_conditions = [["periodic", "periodic"],],
    # Merging parameters
    merging_method = "vranic_spherical",
    merge_every = 1,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = merge_momentum_cell_size,
    merge_discretization_scale = "log",
)

DiagScalar(every = 1,
           vars=['Uelm','Ukin','Utot','Uexp','Ubal',
          'Urad',
          'UmBWpairs',
          'Ukin_electron_cartesian',
          'Ukin_proton_cartesian',
          'Ntot_electron_cartesian',
          'Ntot_proton_cartesian',
          'Dens_electron_cartesian',
          'Dens_proton_cartesian',
          'Ukin_electron_spherical_lin',
          'Ukin_proton_spherical_lin',
          'Ntot_electron_spherical_lin',
          'Ntot_proton_spherical_lin',
          'Dens_electron_spherical_lin',
          'Dens_proton_spherical_lin',
          'Ukin_electron_spherical_log',
          'Ukin_proton_spherical_log',
          'Ntot_electron_spherical_log',
          'Ntot_proton_spherical_log',
          'Dens_electron_spherical_log',
          'Dens_proton_spherical_log'])

for species in ["proton_cartesian",
            "proton_spherical_lin",
            "proton_spherical_log"]:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 50,
		time_average = 1,
		species = [species],
		axes = [
		    ["gamma", 1., 1.0003, 128],
		]
	)
    
for species in ["electron_cartesian",
            "electron_spherical_lin",
            "electron_spherical_log"]:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 50,
		time_average = 1,
		species = [species],
		axes = [
		    ["gamma", 1., 5., 128],
		]
	)
