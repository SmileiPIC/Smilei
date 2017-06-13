# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1d3v",

    number_of_patches = [ 4 ],

    interpolation_order = 2,

    timestep = 0.01 * L0,
    sim_time = 10 * L0,


    time_fields_frozen = 100000000000.,

    cell_length = [20.*L0],
    sim_length = [1600.*L0],

    bc_em_type_x = ["periodic"],


    random_seed = 0,

	referenceAngularFrequency_SI = L0 * 3e8 /1.e-6,
    print_every = 100,
)



electrons = []
energy = []
emin = 1e-2 #keV
emax = 1e10
npoints = 20
for i in range(npoints):
	el = "electron"+str(i)
	electrons .append(el)
	E = math.exp(math.log(emin) + float(i)/npoints*math.log(emax/emin)) #logscale
	E /= 511.
	energy.append(E)
	vel = math.sqrt(1.-1./(1.+E)**2)
	Species(
		species_type = el,
		initPosition_type = "regular",
		initMomentum_type = "maxwell-juettner",
		n_part_per_cell= 10,
		mass = 1.0,
		charge = -1.0,
		charge_density = 1e-9,
		mean_velocity = [vel, 0., 0.],
		temperature = [0.0000001]*3,
		time_frozen = 100000000.0,
		bc_part_type_xmin = "none",
		bc_part_type_xmax = "none",
		bc_part_type_ymin = "none",
		bc_part_type_ymax = "none",
		c_part_max = 10.
	)

Species(
	species_type = "ion1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 100,
	mass = 1836.0*27.,
	charge = 0,
	nb_density = 1.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.00000001]*3,
	time_frozen = 100000000.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none",
	bc_part_type_ymin = "none",
	bc_part_type_ymax = "none",
	atomic_number = 13
)


Collisions(
	species1 = electrons,
	species2 = ["ion1"],
	coulomb_log = 0.,
	ionizing = True
)




DiagFields(
	every = 1000000
)


DiagScalar(
	every = 1000000000
)



for i in range(npoints):
	DiagParticleBinning(
		output = "p_density",
		every = 20,
		species = [electrons[i]],
		axes = [
			 ["x",    0.,    Main.sim_length[0],   1]
		]
	)
	DiagParticleBinning(
		output = "density",
		every = 20,
		species = [electrons[i]],
		axes = [
			 ["x",    0.,    Main.sim_length[0],   1]
		]
	)

