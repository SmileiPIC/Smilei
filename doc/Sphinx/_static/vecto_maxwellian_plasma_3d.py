# __________________________________________________________
#
# Maxwellian homogeneous plasma namelist for SMILEI
# This script is used for the vectorization study
# __________________________________________________________

import os
import numpy as np

c = 299792458
lambdar = 1e-6
wr = 2*np.pi*c/lambdar

Te   = 100./511.   				# electron & ion temperature in me c^2
Ti   = 10./511.   				# electron & ion temperature in me c^2

n0  = 1.

lambdap = 2*np.pi/n0

Lde = np.sqrt(Te)					# Debye length in units of c/\omega_{pe}

cell_length = [0.5*Lde,0.5*Lde,0.5*Lde]

# timestep (0.95 x CFL)
dt  = 0.95 * cell_length[0]/np.sqrt(3.)

# Number of particles per cells
particles_per_cell = 10

# Vectorization
vectorization = "adaptive"

# Small patches that can fit in cache
cells_per_patch = [8,8,8]

# 24 cores per sockets so at least 48 patches per node to fill all cores
# Depends on the architecture
patches_per_node = [5,8,8]

# Simulation time
simulation_time  = 100*dt

# Particle initialization
position_initialization = 'random'

# Density profile
def n0_(x,y,z):
   return n0

# Number of nodes
nodes = [1,1,1]

# Number of patches
patches = [nodes[i]*patches_per_node[i] for i in range(3)]

# grid length
grid_length = [cells_per_patch[i]*patches[i]*cell_length[i] for i in range(3)]


Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2,
    
    timestep = dt,
    simulation_time = simulation_time,
    
    cell_length  = cell_length,
    grid_length = grid_length,
    
    number_of_patches = patches,
    
    EM_boundary_conditions = [ ["periodic"] ],

    patch_arrangement = "linearized_XYZ",

    random_seed = smilei_mpi_rank
)

Vectorization(
    mode = vectorization,
)

Species(
    name = "proton",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    charge_density = n0,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [Ti],
    pusher = "boris",
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
)
Species(
    name = "electron",
    position_initialization = position_initialization,
    momentum_initialization = "mj",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [Te],
    pusher = "boris",
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
)
