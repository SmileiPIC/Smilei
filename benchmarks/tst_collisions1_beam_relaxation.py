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
	sim_time = 15 * L0,
	
	
	time_fields_frozen = 100000000000.,
	
	cell_length = [2.*L0],
	sim_length = [112.*L0],
	
	bc_em_type_x = ["periodic"],
	
	
	random_seed = 0,
	
	referenceAngularFrequency_SI = L0 * 3e8 /1.e-6,
	print_every = 10,
)

i = 0
for ion_nppc, eon_nppc in [[1000, 1000], [1000, 100], [100, 1000]]:
	
	ion = "ion"+str(i)
	eon = "eon"+str(i)
	
	Species(
		species_type = ion,
		initPosition_type = "regular",
		initMomentum_type = "maxwell-juettner",
		n_part_per_cell = ion_nppc,
		mass = 10., #1836.0,
		charge = 1.0,
		nb_density = 10.,
		mean_velocity = [0., 0., 0.],
		temperature = [0.00002],
		time_frozen = 100000000.0,
		bc_part_type_xmin = "none",
		bc_part_type_xmax = "none"
	)
	
	Species(
		species_type = eon,
		initPosition_type = "regular",
		initMomentum_type = "maxwell-juettner",
		n_part_per_cell= eon_nppc,
		mass = 1.0,
		charge = -1.0,
		nb_density = 10.,
		mean_velocity = [0.05, 0., 0.],
		temperature = [0.0000002],
		time_frozen = 100000000.0,
		bc_part_type_xmin = "none",
		bc_part_type_xmax = "none"
	)
	
	Collisions(
		species1 = [eon],
		species2 = [ion],
		coulomb_log = 3
	)
	
	DiagParticleBinning(
		output = "density",
		every = 4,
		time_average = 1,
		species = [eon],
		axes = [
			 ["x",    0*L0,    Main.sim_length[0],   10],
			 ["vx",  -0.1,  0.1,    1000]
		]
	)
	
	DiagParticleBinning(
		output = "density",
		every = 4,
		time_average = 1,
		species = [eon],
		axes = [
			 ["x",    0*L0,    Main.sim_length[0],   10],
			 ["vperp2",  0,  0.01,    1000]
		]
	)
	
	DiagParticleBinning(
		output = "density",
		every = 4,
		time_average = 1,
		species = [ion],
		axes = [
			 ["x",    0*L0,    Main.sim_length[0],   10],
			 ["vx",  -0.1,  0.1,  1000]
		]
	)
	
	i += 1



