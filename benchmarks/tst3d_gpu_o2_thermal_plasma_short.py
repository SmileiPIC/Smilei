# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------


import math as m
import numpy as np
import os

c = 299792458
lambdar = 1e-6                  # reference wavelength
wr = 2*m.pi*c/lambdar

Te   = 100./511.   				# electron & ion temperature in me c^2
Ti   = 10./511.   				# electron & ion temperature in me c^2

n0  = 1.

lambdap = 2*m.pi/n0             # plasma wavelength

Lde = m.sqrt(Te)                # Debye length in units of c/\omega_{pe}

dx = 0.5*Lde
dy = dx
dz = dx
dt  = 0.5 * dx /m.sqrt(3.)		# timestep (0.95 x CFL)

Lx = 32*dx
Ly = 32*dy
Lz = 32*dz

# Simulation time
simulation_time  = 2001*dt

particles_per_cell = 8

number_of_patches = [4,4,4]

position_initialization = 'random'

vectorization = "off"

Main(
    geometry = "3Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = simulation_time,

    cell_length  = [dx,dy,dz],
    grid_length = [Lx,Ly,Lz],

    number_of_patches = number_of_patches,

    EM_boundary_conditions = [ ["periodic"] ],

    print_every = 100,

    gpu_computing = True,

    random_seed = smilei_mpi_rank,
)

Vectorization(
   mode=vectorization,
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
    position_initialization = "proton",
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


DiagFields(
    every = 500
)

DiagScalar(every = 10)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 50,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., Lx, 128],
        ["y", 0., Ly, 128],
        ["z", 0., Lz, 128]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 50,
    time_average = 1,
    species = ["proton"],
    axes = [
        ["x", 0., Lx, 128],
        ["y", 0., Ly, 128],
        ["z", 0., Lz, 128]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 50,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., Lx, 128],
        ["y", 0., Ly, 128]
    ]
)
