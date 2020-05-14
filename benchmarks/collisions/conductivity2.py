# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

from math import pi, sqrt
L0 = 2.*pi # conversion from normalization length to wavelength


Main(
    geometry = "1Dcartesian",

    number_of_patches = [ 4 ],

    interpolation_order = 2,

    timestep = 0.0004 * L0,
    simulation_time = 0.4 * L0,


    time_fields_frozen = 100000000000.,

    cell_length = [0.002*L0],
    grid_length = [10*L0],

    EM_boundary_conditions = [ ["periodic"] ],


    random_seed = 0,

    reference_angular_frequency_SI = L0 * 3e8 /1.e-6,
    print_every = 10,

    solve_poisson = False
)


# EXTERNAL FIELDS
ExternalField(
	field = "Ex",
	profile = 0.01
)


ion_nppc = 10
eon_nppc = 10

charge_density = [415., 554., 757]
charge = [5.6, 7.4, 10.]
temperature = [4e-5, 1e-4, 2e-4] # 20eV, 50eV, 100eV

for i in range(3):
	Species(
		name = "copper"+str(i+1),
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell = ion_nppc,
		mass = 115845.,      # =  mass of Cu atom
		charge = charge[i],
		charge_density = constant(charge_density[i]),
		mean_velocity = [0., 0., 0.],
		temperature = [temperature[i]],
		time_frozen = 0.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
	)
	Species(
		name = "electron"+str(i+1),
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= eon_nppc,
		mass = 1.0,
		charge = -1.0,
		charge_density = constant(charge_density[i]),
		mean_velocity = [0., 0., 0.],
		temperature = [temperature[i]],
		time_frozen = 0.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
	)


	Collisions(
		species1 = ["copper"+str(i+1)],
		species2 = ["electron"+str(i+1)],
		coulomb_log = 2,
		debug_every = 10
	)

	DiagParticleBinning(
		deposited_quantity = "weight_charge_vx",
		every = 10,
		time_average = 10,
		species = ["electron"+str(i+1)],
		axes = [
			 ["x",  0, Main.grid_length[0], 1]
		]
	)
	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 10,
		time_average = 10,
		species = ["electron"+str(i+1)],
		axes = [
			 ["x",  0, Main.grid_length[0], 1]
		]
	)
	
	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 5,
		time_average = 1,
		species = ["electron"+str(i+1)],
		axes = [
			 ["vx",  -5.*sqrt(temperature[i]), 5.*sqrt(temperature[i]), 100]
		]
	)



DiagFields(
	every = 50
)


DiagScalar(
	every = 1
)
