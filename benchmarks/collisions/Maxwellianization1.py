# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
	geometry = "1Dcartesian",
	
	number_of_patches = [ 8 ],
	
	interpolation_order = 2,
	
	timestep = 0.002 * L0,
	simulation_time = 0.5 * L0,
	
	time_fields_frozen = 100000000000.,
	
	cell_length = [20.*L0],
	grid_length = [1600.*L0],
	
	EM_boundary_conditions = [ ["periodic"] ],
	
	random_seed = 0,
	
	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,
	print_every = 10,
)


Species(
	name = "electron1",
	position_initialization = "regular",
	momentum_initialization = "rectangular",
	particles_per_cell= 20000,
	mass = 1.0,
	charge = -1.0,
	charge_density = 10.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.0002, 0.0, 0.0],
	time_frozen = 100000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
)

Collisions(
	species1 = ["electron1"],
	species2 = ["electron1"],
	coulomb_log = 3,
	debug_every = 10
)




DiagFields(
	every = 1
)


DiagScalar(
	every = 1
)


DiagParticleBinning(
	deposited_quantity = "weight",
	every = 5,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.grid_length[0],   10],
		 ["vx",  -0.02,  0.02,    1000]
	]
)
DiagParticleBinning(
	deposited_quantity = "weight",
	every = 5,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.grid_length[0],   10],
		 ["vy",  -0.02,  0.02,    1000]
	]
)
DiagParticleBinning(
	deposited_quantity = "weight",
	every = 5,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.grid_length[0],   10],
		 ["vz",  -0.02,  0.02,    1000]
	]
)
