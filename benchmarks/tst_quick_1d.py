# ----------------------------------------------------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
#
# Quick test 1d
# This benchmarks is used to test a lot of features in a minimum of compute time
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
simulation_time = 10 * timestep
# Period of output for the diags
diag_every = int(simulation_time / timestep)

Main(
    geometry = "1Dcartesian",
    interpolation_order = 2 ,
    cell_length = cell_length,
    grid_length  = grid_length,
    number_of_patches = number_of_patches,
    timestep = timestep,
    simulation_time = simulation_time,
    EM_boundary_conditions = [
        ['silver-muller'],
    ],
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
	pusher = "boris",
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
	pusher = "vay",
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
	pusher = "higueracary",
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

DiagScalar(
	every=1,
)

for species_name in ["eon1","eon2"]:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 10,
		time_average = 1,
		species = [species_name],
		axes = [
			["gamma", 15., 30., 128],
		]
	)

for species_name in ["pon1","pon2"]:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 10,
		time_average = 1,
		species = [species_name],
		axes = [
			["gamma", 22.2, 22.8, 128],
		]
	)

species_list = ["eon1","eon2","pon1","pon2"]

fields = ['Ex','Ey','Ez',
			'Bx','By','Bz',
			'Rho','Jx','Jy','Jz']

for species in species_list:
    fields.append("Jx_{}".format(species))	
    fields.append("Jy_{}".format(species))	
    fields.append("Jz_{}".format(species))	
    fields.append("Rho_{}".format(species))	

DiagFields(
    every = 10,
    fields = fields,
)

fields = ['Ex','Ey','Ez',
			'Bx','By','Bz',
			'Rho','Jx','Jy','Jz']

for species in species_list:
    fields.append("Jx_{}".format(species))	
    fields.append("Jy_{}".format(species))	
    fields.append("Jz_{}".format(species))	
    fields.append("Rho_{}".format(species))	

DiagProbe(
    every = 10,
    origin = [0.],
    corners = [[Main.grid_length[0]]],
    number = [256],
    fields = fields,
)
