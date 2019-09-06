# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
	geometry = "2Dcartesian",
	
	number_of_patches = [ 4, 4 ],
	
	interpolation_order = 2,
	
	timestep = 1 * L0,
	simulation_time = 10 * L0,
	
	time_fields_frozen = 100000000000.,
	
	cell_length = [L0, L0],
	grid_length = [64.*L0, 64*L0],
	
	EM_boundary_conditions = [ ["periodic"] ],
	
	solve_poisson = False,
	
	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,

	random_seed = smilei_mpi_rank	
)


ion = "ion"
eon = "eon"

Species(
	name = eon,
	position_initialization = "random",
	momentum_initialization = "cold",
	particles_per_cell= 100,
	mass = 1.0,
	charge = -1.0,
	charge_density = 1.,
	mean_velocity = [0.8, 0., 0.],
	temperature = [0.]*3,
	time_frozen = 100000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
)

Species(
	name = ion,
	position_initialization = "random",
	momentum_initialization = "cold",
	particles_per_cell= 100,
	mass = 1836.0*13.,
	charge = 3.0,
	charge_density = 1.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.]*3,
	time_frozen = 100000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
	atomic_number = 13
)

Collisions(
	species1 = [eon],
	species2 = [ion],
	coulomb_log = 3,
)

DiagParticleBinning(
	deposited_quantity = "weight_vy_py",
	every = 4,
	species = [ion],
	axes = [
		 ["x",    0,    Main.grid_length[0],   64],
		 ["y",    0,    Main.grid_length[1],   64]
	]
)
