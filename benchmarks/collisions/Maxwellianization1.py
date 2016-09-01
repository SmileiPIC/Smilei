# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1d3v",

    number_of_patches = [ 16 ],

    interpolation_order = 2,

    timestep = 0.002 * L0,
    sim_time = 0.5 * L0,

    time_fields_frozen = 100000000000.,

    cell_length = [20.*L0],
    sim_length = [1600.*L0],

    bc_em_type_x = ["periodic"],


    random_seed = 0,

	referenceAngularFrequency_SI = L0 * 3e8 /1.e-6,
    print_every = 10,
)


Species(
	species_type = "electron1",
	initPosition_type = "regular",
	initMomentum_type = "rectangular",
	n_part_per_cell= 20000,
	mass = 1.0,
	charge = -1.0,
	charge_density = 10.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.0002, 0.0, 0.0],
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
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


DiagParticles(
	output = "density",
	every = 5,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.sim_length[0],   10],
		 ["vx",  -0.02,  0.02,    1000]
	]
)
DiagParticles(
	output = "density",
	every = 5,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.sim_length[0],   10],
		 ["vy",  -0.02,  0.02,    1000]
	]
)
DiagParticles(
	output = "density",
	every = 5,
	species = ["electron1"],
	axes = [
		 ["x",    0*L0,    Main.sim_length[0],   10],
		 ["vz",  -0.02,  0.02,    1000]
	]
)
