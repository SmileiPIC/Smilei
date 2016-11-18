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
    sim_time = 0.4 * L0,


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
	field = "Ex",
	profile = 0.001
)

Species(
	species_type = "copper1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 1000,
	mass = 115845.,      # =  mass of Cu atom
	charge = 4.4,
	charge_density = constant(330.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.000002], # 1eV
	time_frozen = 0.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
)
Species(
	species_type = "electron1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 10000,
	mass = 1.0,
	charge = -1.0,
	charge_density = constant(330.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.000002], # 1 eV
	time_frozen = 0.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
)

Species(
	species_type = "copper2",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 1000,
	mass = 115845.,      # =  mass of Cu atom
	charge = 4.4,
	charge_density = constant(333.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.000006], # 3eV
	time_frozen = 0.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
)
Species(
	species_type = "electron2",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 10000,
	mass = 1.0,
	charge = -1.0,
	charge_density = constant(333.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.000006], # 3 eV
	time_frozen = 0.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
)

Species(
	species_type = "copper3",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 1000,
	mass = 115845,      # =  mass of Cu atom
	charge = 5.,
	charge_density = constant(368.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.00002], # 10eV
	time_frozen = 0.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
)
Species(
	species_type = "electron3",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 10000,
	mass = 1.0,
	charge = -1.0,
	charge_density = constant(368.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.00002], # 10 eV
	time_frozen = 0.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
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
Collisions(
	species1 = ["copper3"],
	species2 = ["electron3"],
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
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)
DiagParticles(
	output = "jx_density",
	every = 5,
	time_average = 1,
	species = ["electron2"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)
DiagParticles(
	output = "jx_density",
	every = 5,
	time_average = 1,
	species = ["electron3"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)

DiagParticles(
	output = "density",
	every = 5,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)
DiagParticles(
	output = "density",
	every = 5,
	time_average = 1,
	species = ["electron2"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)
DiagParticles(
	output = "density",
	every = 5,
	time_average = 1,
	species = ["electron3"],
	axes = [
		 ["x",  0, 20*L0, 1]
	]
)
