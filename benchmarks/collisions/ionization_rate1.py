# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1Dcartesian",

    number_of_patches = [ 4 ],

    interpolation_order = 2,

    timestep = 0.05 * L0,
    simulation_time = 100 * L0,


    time_fields_frozen = 100000000000.,

    cell_length = [20.*L0],
    grid_length = [1600.*L0],

    EM_boundary_conditions = [ ["periodic"] ],


    random_seed = 0,

	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,
    print_every = 100,
)


Species(
	name = "electron1",
	position_initialization = "regular",
	momentum_initialization = "maxwell-juettner",
	particles_per_cell= 100,
	mass = 1.0,
	charge = -1.0,
	charge_density = 1.,
	mean_velocity = [0.03, 0., 0.],
	temperature = [0.0000001]*3,
	time_frozen = 100000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
	c_part_max = 10.
)

Species(
	name = "ion1",
	position_initialization = "regular",
	momentum_initialization = "maxwell-juettner",
	particles_per_cell= 100,
	mass = 1836.0*13.,
	charge = 3.0,
	charge_density = 1.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.00000001]*3,
	time_frozen = 100000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
	atomic_number = 13
)


Collisions(
	species1 = ["electron1"],
	species2 = ["ion1"],
	coulomb_log = 3,
	ionizing = True
)




DiagFields(
	every = 1
)


DiagScalar(
	every = 1
)



DiagParticleBinning(
	deposited_quantity = "weight_charge",
	every = 20,
	species = ["ion1"],
	axes = [
		 ["x",    0,    Main.grid_length[0],   1]
	]
)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = 20,
	species = ["ion1"],
	axes = [
		 ["x",    0,    Main.grid_length[0],   1]
	]
)


DiagParticleBinning(
	deposited_quantity = "weight",
	every = 20,
	species = ["electron1"],
	axes = [
		 ["px",    -0.3,    0.3,   1000]
	]
)

