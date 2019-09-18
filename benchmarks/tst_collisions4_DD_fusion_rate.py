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
	simulation_time = 10 * L0,
	
	
	time_fields_frozen = 100000000000.,
	
	cell_length = [20.*L0],
	grid_length = [1600.*L0],
	
	EM_boundary_conditions = [ ["periodic"] ],
	
	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,

	random_seed = smilei_mpi_rank	
)


i = 0
for Da_nppc, Db_nppc, v in [
	[100, 100, 0.1],
	[100,  10, 0.1],
	[10 , 100, 0.1],
	[100, 100, 0.00082],
	[100, 100, 0.00201],
	[100, 100, 0.00493],
	[100, 100, 0.01209],
	[100, 100, 0.02960],
	[100, 100, 0.07236],
	[100, 100, 0.17545],
	[100, 100, 0.40641],
	[100, 100, 0.76966],
	[100, 100, 0.97376],
	]:
	
	Da = "Da_"+str(i)
	Db = "Db_"+str(i)
	He = "He_"+str(i)
	
	Species(
		name = Da,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= Da_nppc,
		mass = 3870.5,
		charge = 0.,
		number_density = 100.,
		mean_velocity = [v, 0., 0.],
		temperature = [0.0000001]*3,
		time_frozen = 100000000.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
	)
	
	Species(
		name = Db,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= Db_nppc,
		mass = 3870.5,
		charge = 0.,
		number_density = 100.,
		mean_velocity = [-v, 0., 0.],
		temperature = [0.00000001]*3,
		time_frozen = 100000000.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
	)
	
	Species(
		name = He,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= 0,
		mass = 5497.9,
		charge = 0.,
		number_density = 0.,
		time_frozen = 100000000.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
	)
	
	Collisions(
		species1 = [Da],
		species2 = [Db],
		coulomb_log = 0.001,
		nuclear_reaction = ["D-D", He],
	)
	
	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 20,
		species = [He],
		axes = [
			 ["x",    0,    Main.grid_length[0],   1]
		]
	)
	
	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 20,
		species = [Da, Db],
		axes = [
			 ["x",    0,    Main.grid_length[0],   1]
		]
	)
	
	i+=1

DiagScalar(
	every = 50,
)