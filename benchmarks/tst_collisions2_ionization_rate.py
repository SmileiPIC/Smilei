# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
	geometry = "1d3v",
	
	number_of_patches = [ 4 ],
	
	interpolation_order = 2,
	
	timestep = 0.05 * L0,
	sim_time = 100 * L0,
	
	
	time_fields_frozen = 100000000000.,
	
	cell_length = [20.*L0],
	sim_length = [1600.*L0],
	
	bc_em_type_x = ["periodic"],
	
	
	random_seed = 0,
	
	referenceAngularFrequency_SI = L0 * 3e8 /1.e-6,
)


i = 0
for ion_nppc, eon_nppc in [[100, 100], [100, 10], [10, 100]]:
	
	ion = "ion"+str(i)
	eon = "eon"+str(i)
	
	Species(
		species_type = eon,
		initPosition_type = "regular",
		initMomentum_type = "maxwell-juettner",
		n_part_per_cell= eon_nppc,
		mass = 1.0,
		charge = -1.0,
		charge_density = 1.,
		mean_velocity = [0.03, 0., 0.],
		temperature = [0.0000001]*3,
		time_frozen = 100000000.0,
		bc_part_type_xmin = "none",
		bc_part_type_xmax = "none",
		c_part_max = 10.
	)
	
	Species(
		species_type = ion,
		initPosition_type = "regular",
		initMomentum_type = "maxwell-juettner",
		n_part_per_cell= ion_nppc,
		mass = 1836.0*13.,
		charge = 3.0,
		charge_density = 1.,
		mean_velocity = [0., 0., 0.],
		temperature = [0.00000001]*3,
		time_frozen = 100000000.0,
		bc_part_type_xmin = "none",
		bc_part_type_xmax = "none",
		atomic_number = 13
	)
	
	Collisions(
		species1 = [eon],
		species2 = [ion],
		coulomb_log = 3,
		ionizing = True
	)
	
	DiagParticles(
		output = "charge_density",
		every = 20,
		species = [ion],
		axes = [
			 ["x",    0,    Main.sim_length[0],   1]
		]
	)
	
	DiagParticles(
		output = "density",
		every = 20,
		species = [ion],
		axes = [
			 ["x",    0,    Main.sim_length[0],   1]
		]
	)
	
	i+=1

