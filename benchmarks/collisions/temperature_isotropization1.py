# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1Dcartesian",

    number_of_patches = [ 32 ],

    interpolation_order = 2,

    timestep = 0.01 * L0,
    simulation_time = .4 * L0,


    time_fields_frozen = 100000000000.,

    cell_length = [4.*L0],
    grid_length = [4*1024.*L0],

    EM_boundary_conditions = [ ["periodic"] ],


    random_seed = smilei_mpi_rank+1,

	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,
    print_every = 10,
)


Species(
	name = "electron1",
	position_initialization = "regular",
	momentum_initialization = "maxwell-juettner",
	particles_per_cell= 1000,
	mass = 1.0,
	charge = -1.0,
	charge_density = 10.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.000011, 0.00001, 0.00001],
	time_frozen = 100000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
)

Collisions(
	species1 = ["electron1"],
	species2 = ["electron1"],
	coulomb_log = 2
)


DiagParticleBinning(
	deposited_quantity = "weight",
	every = 1,
	species = ["electron1"],
	axes = [
		 ["vx",  -0.02,  0.02,    1000]
	]
)
DiagParticleBinning(
	deposited_quantity = "weight",
	every = 1,
	species = ["electron1"],
	axes = [
		 ["vy",  -0.02,  0.02,    1000]
	]
)
DiagParticleBinning(
	deposited_quantity = "weight",
	every = 1,
	species = ["electron1"],
	axes = [
		 ["vz",  -0.02,  0.02,    1000]
	]
)
