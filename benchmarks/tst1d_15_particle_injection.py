# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
#
# Particle injection from the Xmin and Xmax boundaries
#
# ----------------------------------------------------------------------------------------

import math
import numpy as np

# Mean velocity
mean_velocity = 0.999
# Electron temperature
Te = 0.01
# Ion temperature
Ti = 0.001
# Ion charge
Zi = 1
# Density
n0 = 1
# Debye length
Debye_length = 1. / np.sqrt( n0 / Te + Zi * n0 / Ti )
# Cell length
cell_length = [Debye_length*0.5]
# Number of patches
number_of_patches =[16]
# Cells per patches (patch shape)
cells_per_patch = [32.]
# Grid length
grid_length = [0.]
for i in range(1):
    grid_length[i] = number_of_patches[i] * cell_length[i] * cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 32
# Position init
position_initialization = 'random'
# Time step
timestep = 0.95 * cell_length[0]
# Total simulation time
simulation_time = ((0.5 - 0.125)*grid_length[0])/mean_velocity          # duration of the simulation
# Period of output for the diags
diag_every = int(simulation_time / timestep)

Main(
    geometry = "1Dcartesian",
    interpolation_order = 2 ,
    cell_length = cell_length,
    grid_length  = grid_length,
    number_of_patches = number_of_patches,
    #cell_sorting = True,
    timestep = timestep,
    simulation_time = simulation_time,
    EM_boundary_conditions = [
        ['silver-muller'],
    ],
    random_seed = smilei_mpi_rank,
)

# Initial plasma shape
fp = trapezoidal(1., xvacuum=0.        ,xplateau=grid_length[0]/8.)
fm = trapezoidal(1., xvacuum=7*grid_length[0]/8.,xplateau=grid_length[0])

Species(
	name = 'pon1',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1836.0,
	charge = 1.0,
	number_density = fp,
	mean_velocity = [mean_velocity,0.,0.],
	temperature = [Ti],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
	],
)

ParticleInjector(
    species = 'pon1',
    box_side = 'xmin',
)

Species(
	name = 'eon1',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fp,
	mean_velocity = [mean_velocity,0.,0.],
	temperature = [Te],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
	],
)
ParticleInjector(
    species = 'eon1',
    box_side = 'xmin',
)

Species(
	name = 'pon2',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1836.0,
	charge = 1.0,
	number_density = fm,
	mean_velocity = [-mean_velocity,0.,0.],
	temperature = [Ti],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
	],
)
ParticleInjector(
    species = 'pon2',
    box_side = 'xmax',
)

Species(
	name = 'eon2',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fm,
	mean_velocity = [-mean_velocity,0.,0.],
	temperature = [Te],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
	],
)
ParticleInjector(
    species = 'eon2',
    box_side = 'xmax',
)

DiagScalar(every=1)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["gamma", 15., 30., 128],
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["gamma", 22., 23., 128],
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon2"],
    axes = [
        ["gamma", 15., 30., 128],
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon2"],
    axes = [
        ["gamma", 22., 23., 128],
    ]
)

# DiagFields(
#     every = globalEvery,
#     fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_pon1','Rho_eon1','Rho_pon2','Rho_eon2',"Jx","Jy","Jz"]
# )

# DiagProbe(
#     every = globalEvery,
#     origin = [0., Main.grid_length[1]/2.],
#     corners = [[Main.grid_length[0],Main.grid_length[1]/2.]],
#     number = [256],
#     fields = ['Ex','Ey','Ez','Rho_pon1','Jx_pon1','Jy_pon1','Jz_pon1']
# )
