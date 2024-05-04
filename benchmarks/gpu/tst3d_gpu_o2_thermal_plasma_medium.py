# ----------------------------------------------------------------------------------------
# 					Thermal plasma 3D (in a cube)
# ----------------------------------------------------------------------------------------

import math as m
import numpy as np

###############################
# HPC settings
###############################

# Patches that can fit in cache/device.
# On CPU, it may not be a bad idea to slightly oversubscribe the cores with more patches.
# On GPU, large patches work better.
kGridDimensionInCell   = [128, 128, 128]
kPatchPerGridDimension = [2, 2, 1]

# Smilei conf 2022 recommends cell length < 4 * Debye length
kCellLengthFactor = 1.0 / 2.0 # / 4.0

###############################
# Physics settings
###############################

kTimeStepCount = 500
c              = 299792458
lambdar        = 1e-6

Te = 100.0 / 511.0 # electron & ion temperature in me c^2
Ti = 10.0 / 511.0  # electron & ion temperature in me c^2

# Debye length in units of c/\omega_{pe}
Lde = m.sqrt(Te) 

# Simulation time
dt = 0.95 * ((Lde * kCellLengthFactor) / m.sqrt(3.0)) # timestep (0.95 x CFL)

kParticlesPerCell = 64

kPositionInitialization = 'regular' 

def InitialChargeDensity(x, y, z):
    # Homogeneous
    n0  = 1.
    return n0

###############################
# Smilei settings
###############################

kSimulationTime = kTimeStepCount * dt	

kCellLength = [Lde * kCellLengthFactor,
               Lde * kCellLengthFactor,
               Lde * kCellLengthFactor]

kGridLength = [kGridDimensionInCell[i] * kCellLength[i] for i in range(3)]

Main(gpu_computing = True,
     geometry = "3Dcartesian",
     interpolation_order = 2,
     timestep = dt,
     simulation_time = kSimulationTime,
     cell_length  = kCellLength,
     grid_length = kGridLength,
     number_of_patches = kPatchPerGridDimension,
     EM_boundary_conditions = [ ["periodic"] ],
     print_every = 10,
     )

Vectorization(mode = "off")

# LoadBalancing(every = 20,
#               initial_balance = False,
#               cell_load = 1.0,
#               frozen_particle_load = 0.1)

Species(name = "proton",
        position_initialization = kPositionInitialization,
        momentum_initialization = "mj",
        particles_per_cell = kParticlesPerCell, 
        c_part_max = 1.0,
        mass = 1836.0,
        charge = 1.0,
        charge_density = InitialChargeDensity,
        mean_velocity = [0.0, 0.0, 0.0],
        temperature = [Te],
        pusher = "boris",
        boundary_conditions = [
        	["periodic", "periodic"],
        	["periodic", "periodic"],
        	["periodic", "periodic"],
        ])

Species(name = "electron",
        position_initialization = kPositionInitialization,
        momentum_initialization = "mj",
        particles_per_cell = kParticlesPerCell, 
        c_part_max = 1.0,
        mass = 1.0,
        charge = -1.0,
        charge_density = InitialChargeDensity,
        mean_velocity = [0.0, 0.0, 0.0],
        temperature = [Ti],
        pusher = "boris",
        boundary_conditions = [
        	["periodic", "periodic"],
        	["periodic", "periodic"],
        	["periodic", "periodic"],
        ])

#Checkpoints(
#    dump_step = 0,
#    dump_minutes = 0,
#    exit_after_dump = False,
#)

#DiagFields(
#    every = 10
#)

DiagScalar(every = 10)

#DiagPerformances(
#    every = 10,
#)

