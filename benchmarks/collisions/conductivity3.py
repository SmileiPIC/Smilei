# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1d3v",

    number_of_patches = [ 4 ],

    interpolation_order = 2,

    timestep = 0.002 * L0,
    sim_time = 0.8 * L0,


    time_fields_frozen = 100000000000.,

    cell_length = [1*L0],
    sim_length = [20*L0],

    bc_em_type_x = ["periodic"],


    random_seed = 0,

	referenceAngularFrequency_SI = L0 * 3e8 /1.e-6,
    print_every = 10,
)


# EXTERNAL FIELDS
ExtField(
	fields = ["Ex"],
	profile = 0.001
)


Species(
	species_type = "copper1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 50000,
	mass = 115845.,      # =  mass of Cu atom
	charge = 17.,
	charge_density = constant(1268.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.0006], # 300 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)
Species(
	species_type = "electron1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 100000,
	mass = 1.0,
	charge = -1.0,
	charge_density = constant(1268.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.0006], # 300 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "copper2",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 50000,
	mass = 115845.,      # =  mass of Cu atom
	charge = 25.,
	charge_density = constant(1869.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.002], # 1000 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)
Species(
	species_type = "electron2",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 100000,
	mass = 1.0,
	charge = -1.0,
	charge_density = constant(1869.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.002], # 1000 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Collisions(
	species1 = ["copper1"],
	species2 = ["electron1"],
	coulomb_log = 2.
)
Collisions(
	species1 = ["copper2"],
	species2 = ["electron2"],
	coulomb_log = 2.
)




DiagFields(
	every = 5
)


DiagScalar(
	every = 1
)



DiagParticles(
	output = "jx_density",
	every = 5,
	time_average = 4,
	species = ["electron1"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)
DiagParticles(
	output = "jx_density",
	every = 5,
	time_average = 4,
	species = ["electron2"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)


DiagParticles(
	output = "density",
	every = 5,
	time_average = 4,
	species = ["electron1"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)
DiagParticles(
	output = "density",
	every = 5,
	time_average = 4,
	species = ["electron2"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)

