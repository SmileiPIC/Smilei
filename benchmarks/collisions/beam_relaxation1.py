# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1d3v",

    number_of_patches = [ 8 ],

    interpolation_order = 2,

    timestep = 0.2 * L0,
    sim_time = 40 * L0,

    time_fields_frozen = 100000000000.,

    cell_length = [2.*L0],
    sim_length = [112.*L0],

    bc_em_type_x = ["periodic"],

    random_seed = 0,

	referenceAngularFrequency_SI = L0 * 3e8 /1.e-6,
    print_every = 10,
)


Species(
	species_type = "ion1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 1000,
	mass = 10., #1836.0,
	charge = 1.0,
	nb_density = 10.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.00002],
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "electron1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 1000,
	mass = 1.0,
	charge = -1.0,
	nb_density = 10.,
	mean_velocity = [0.05, 0., 0.],
	temperature = [0.0000002],
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Collisions(
	species1 = ["electron1"],
	species2 = ["ion1"],
	coulomb_log = 3
)




DiagFields(
	every = 1
)

DiagScalar(
	every = 1
)



DiagParticles(
	output = "density",
	every = 2,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.sim_length[0],   10],
		 ["vx",  -0.1,  0.1,    1000]
	]
)

DiagParticles(
	output = "density",
	every = 2,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.sim_length[0],   10],
		 ["vperp2",  0,  0.01,    1000]
	]
)

DiagParticles(
	output = "density",
	every = 2,
	time_average = 1,
	species = ["ion1"],
	axes = [
		 ["x",    0*L0,    Main.sim_length[0],   10],
		 ["vx",  -0.1,  0.1,  1000]
	]
)

DiagParticles(
	output = "density",
	every = 10,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["ekin",  0.0001,  0.1, 100, "logscale"]
	]
)
