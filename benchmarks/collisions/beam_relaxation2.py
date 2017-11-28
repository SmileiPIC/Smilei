# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1Dcartesian",

    number_of_patches = [ 8 ],

    interpolation_order = 2,

    timestep = 0.2 * L0,
    simulation_time = 40 * L0,


    time_fields_frozen = 100000000000.,

    cell_length = [2.*L0],
    grid_length = [112.*L0],

    EM_boundary_conditions = [ ["periodic"] ],


    random_seed = 0,

	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,
    print_every = 10,
)


Species(
	name = "ion1",
	position_initialization = "regular",
	momentum_initialization = "maxwell-juettner",
	particles_per_cell = 100,
	mass = 10., #1836.0,
	charge = 1.0,
	number_density = 10.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.00002],
	time_frozen = 100000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
)

Species(
	name = "electron1",
	position_initialization = "regular",
	momentum_initialization = "maxwell-juettner",
	particles_per_cell= 1000,
	mass = 1.0,
	charge = -1.0,
	number_density = 10.,
	mean_velocity = [0.05, 0., 0.],
	temperature = [0.0000002],
	time_frozen = 100000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
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



DiagParticleBinning(
	deposited_quantity = "weight",
	every = 2,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.grid_length[0],   10],
		 ["vx",  -0.1,  0.1,    1000]
	]
)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = 2,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.grid_length[0],   10],
		 ["vperp2",  0,  0.01,    1000]
	]
)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = 2,
	time_average = 1,
	species = ["ion1"],
	axes = [
		 ["x",    0*L0,    Main.grid_length[0],   10],
		 ["vx",  -0.1,  0.1,  1000]
	]
)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = 10,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["ekin",  0.0001,  0.1, 100, "logscale"]
	]
)

