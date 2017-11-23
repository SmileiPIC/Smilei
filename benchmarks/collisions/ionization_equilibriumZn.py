# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
    geometry = "1Dcartesian",

    number_of_patches = [ 4 ],

    interpolation_order = 2,

    timestep = 1 * L0,
    simulation_time = 400 * L0,


    time_fields_frozen = 100000000000.,

    cell_length = [5.*L0],
    grid_length = [160.*L0],

    EM_boundary_conditions = [ ["periodic"] ],


    random_seed = 0,

	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,
    print_every = 100,
)



Z = 30
A = 65.
density = 10.

electrons = []
ions = []
temperature = []
Tmin = 0.001 # keV
Tmax = 2.
npoints = 10

for i in range(npoints):
	eon = "electron"+str(i)
	ion = "ion"+str(i)
	electrons.append(eon)
	ions.append(ion)
	
	T = math.exp(math.log(Tmin) + float(i)/(npoints-1)*math.log(Tmax/Tmin)) #logscale
	T /= 511.
	temperature.append(T)
	
	Zstar = (Z)*(1.-math.exp(-T/(1./511.)))
	Zstar = round(Zstar, 2)
	
	Species(
		name = eon,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= 100,
		mass = 1.0,
		charge = -1.0,
		charge_density = Zstar*density,
		mean_velocity = [0., 0., 0.],
		temperature = [T]*3,
		time_frozen = 100000000.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
		c_part_max = 10.
	)
	
	Species(
		name = ion,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= 100,
		mass = 1836.0*A,
		charge = Zstar,
		number_density = 10.,
		mean_velocity = [0., 0., 0.],
		temperature = [T]*3,
		time_frozen = 100000000.0,
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
		deposited_quantity = "weight_ekin",
		every = 10,
		species = [eon],
		axes = [ ["x", 0, Main.grid_length[0], 1] ]
	)
	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 10,
		species = [eon],
		axes = [ ["x", 0, Main.grid_length[0], 1] ]
	)
	DiagParticleBinning(
		deposited_quantity = "weight_charge",
		every = 10,
		species = [ion],
		axes = [ ["x", 0, Main.grid_length[0], 1] ]
	)
	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 10,
		species = [ion],
		axes = [ ["x", 0, Main.grid_length[0], 1] ]
	)
	#DiagParticleBinning(
	#	deposited_quantity = "weight",
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





