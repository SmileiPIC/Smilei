# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math

L0 = 2.*math.pi

Main(
    geometry = "1d3v",
    
    interpolation_order = 2,
    
    timestep = 0.005 * L0,
    sim_time  = 0.5 * L0,
    
    time_fields_frozen = 100000000000.,
    
    
    cell_length = [0.01 * L0],
    sim_length  = [1. * L0],
    
    number_of_patches = [ 4 ],
    
    bc_em_type_x  = ["periodic"],
    
    referenceAngularFrequency_SI = L0 * 3e8 /1.e-6,
    
    random_seed = 0,
    
    print_every = 10
)


Species(
	species_type = "ion1",
	initPosition_type = "random",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 2000,
	mass = 1836.0,
	charge = 1.0,
	nb_density = 10.,
	temperature = [0.00002],
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "electron1",
	initPosition_type = "random",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 2000,
	mass = 1.0,
	charge = -1.0,
	nb_density = 10.,
	mean_velocity = [0.05, 0., 0.],
	temperature = [0.00002],
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)


DiagFields(
    every = 5,
)

DiagScalar(every = 1)


DiagParticles(
	output = "density",
	every = 4,
	time_average = 2,
	species = ["electron1"],
	axes = [
		["x", 0.*L0, 1.*L0, 100],
		["vx", -0.1, 0.1, 100]
	]
)

DiagParticles(
	output = "density",
	every = 4,
	time_average = 1,
	species = ["ion1"],
	axes = [
		("x", 0.*L0, 1.*L0, 100),
		("vx", -0.001, 0.001, 100)
	]
)

DiagParticles(
	output = "px_density",
	every = 4,
	time_average = 2,
	species = ["electron1"],
	axes = [
		["x", 0.*L0, 1.*L0, 100],
		["vx", -0.1, 0.1, 100]
	]
)

DiagParticles(
	output = "density",
	every = 1,
	time_average = 1,
	species = ["electron1"],
	axes = [
		["ekin", 0.0001, 0.1, 100, "logscale", "edge_inclusive"]
	]
)
