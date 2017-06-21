# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1d3v",

    number_of_patches = [ 4 ],

    interpolation_order = 2,

    timestep = 1 * L0,
    sim_time = 400 * L0,


    time_fields_frozen = 100000000000.,

    cell_length = [5.*L0],
    sim_length = [160.*L0],

    EM_boundary_conditions = [ ["periodic"] ],


    random_seed = 0,

	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,
    print_every = 100,
)



Z = 1
A = 1
density = 10.

electrons = []
ions = []
temperature = []
Tmin = 0.003 # keV
Tmax = 0.4
npoints = 10

for i in range(npoints):
	eon = "electron"+str(i)
	ion = "ion"+str(i)
	electrons.append(eon)
	ions.append(ion)
	
	T = math.exp(math.log(Tmin) + float(i)/(npoints-1)*math.log(Tmax/Tmin)) #logscale
	T /= 511.
	temperature.append(T)
	
	Zstar = (Z)*(1.-math.exp(-T/(0.03/511.)))**2
	Zstar = round(Zstar, 2)
	
	Species(
		species_type = eon,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		n_part_per_cell= 1000,
		mass = 1.0,
		charge = -1.0,
		charge_density = Zstar*density,
		mean_velocity = [0., 0., 0.],
		temperature = [T]*3,
		time_frozen = 100000000.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
		boundary_conditions = [
			["periodic", "periodic"],
		],
		c_part_max = 10.
	)
	
	Species(
		species_type = ion,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		n_part_per_cell= 1000,
		mass = 1836.0*A,
		charge = Zstar,
		nb_density = 10.,
		mean_velocity = [0., 0., 0.],
		temperature = [T]*3,
		time_frozen = 100000000.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
		boundary_conditions = [
			["periodic", "periodic"],
		],
		atomic_number = Z
	)
	
	Collisions(
		species1 = [eon],
		species2 = [ion],
		coulomb_log = 0.,
		ionizing = True
	)
	
	DiagParticleBinning(
		output = "ekin_density",
		every = 10,
		species = [eon],
		axes = [ ["x", 0, Main.sim_length[0], 1] ]
	)
	DiagParticleBinning(
		output = "density",
		every = 10,
		species = [eon],
		axes = [ ["x", 0, Main.sim_length[0], 1] ]
	)
	DiagParticleBinning(
		output = "charge_density",
		every = 10,
		species = [ion],
		axes = [ ["x", 0, Main.sim_length[0], 1] ]
	)
	DiagParticleBinning(
		output = "density",
		every = 10,
		species = [ion],
		axes = [ ["x", 0, Main.sim_length[0], 1] ]
	)
	#DiagParticleBinning(
	#	output = "density",
	#	every = 50,
	#	species = [ion],
	#	axes = [ ["charge", -0.5, Z+0.5, Z+1] ]
	#)




DiagFields(
	every = 1000000
)


DiagScalar(
	every = 10000000
)





