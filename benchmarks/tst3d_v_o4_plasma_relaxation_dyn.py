# ------------------------------------------------------------------------------
# Relaxation of a thermal electron-positron plasma in 3D
# ------------------------------------------------------------------------------

import math as m
import numpy as np

c = 299792458
lambdar = 1e-6
wr = 2*m.pi*c/lambdar

Te   = 100./511.   				# electron & ion temperature in me c^2
Ti   = 100./511.   				# electron & ion temperature in me c^2

n0  = 1.

lambdap = 2*m.pi/n0

# Debye length in units of c/\omega_{pe}
Lde = m.sqrt(Te)

cell_length = [0.5*Lde,0.5*Lde,0.5*Lde]

# timestep (0.95 x CFL)
dt  = 0.95 * cell_length[0]/m.sqrt(3.)

# Small patches that can fit in cache
cells_per_patch = [10,10,10]

# Number of patches per node
patches_per_node = [4,4,4]

# Number of nodes
nodes = [1,1,1]

# Number of patches
patches = [nodes[i]*patches_per_node[i] for i in range(3)]

# grid length
grid_length = [cells_per_patch[i]*patches[i]*cell_length[i] for i in range(3)]

# Simulation time
simulation_time  = 40*dt

particles_per_cell = 32

position_initialization = 'random'

EM_boundary_conditions = [["periodic"]]

species_boundary_conditions = [["periodic"]]

def n0_(x,y,z):
  if ((grid_length[0]*0.5 < x and x < grid_length[0]*1) and
      (grid_length[1]*0.5 < y and y < grid_length[1]*1) and
      (grid_length[2]*0.5 < z and z < grid_length[2]*1)):
    return n0
  else:
    return 0.

Main(
    geometry = "3Dcartesian",

    interpolation_order = 4,

    timestep = dt,
    simulation_time = simulation_time,

    cell_length  = cell_length,
    grid_length = grid_length,

    number_of_patches = patches,

    EM_boundary_conditions = EM_boundary_conditions,

    every_clean_particles_overhead = 1,

    patch_arrangement = "linearized_XYZ",

)

Vectorization(
    mode = "adaptive",
)

Species(
    name = "positron",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = n0_,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [Ti],
    pusher = "boris",
    boundary_conditions = species_boundary_conditions,
)
Species(
    name = "electron",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [Te],
    pusher = "boris",
    boundary_conditions = species_boundary_conditions,
)

DiagFields(
    every = 5
)

DiagPerformances(
    every = 5,
)

DiagScalar(every = 5)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., grid_length[0], 32],
        ["y", 0., grid_length[1], 32],
        ["z", 0., grid_length[2], 32]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    time_average = 1,
    species = ["positron"],
    axes = [
        ["x", 0., grid_length[0], 32],
        ["y", 0., grid_length[1], 32],
        ["z", 0., grid_length[2], 32]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 1,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., grid_length[0], 32],
        ["y", 0., grid_length[1], 32],
        ["z", 0., grid_length[2], 32]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 1,
    time_average = 1,
    species = ["positron"],
    axes = [
        ["x", 0., grid_length[0], 32],
        ["y", 0., grid_length[1], 32],
        ["z", 0., grid_length[2], 32]
    ]
)
