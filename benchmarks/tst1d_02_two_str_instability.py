# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
L  = 1.12			# wavelength=simulation box length
dn = 0.001			# amplitude of the perturbation


Main(
    geometry = "1Dcartesian",
     
    interpolation_order = 2,
     
    cell_length = [0.01],
    grid_length  = [L],
    
    number_of_patches = [ 16 ],
    
    timestep = 0.0095,
    simulation_time = 50.,
     
    EM_boundary_conditions = [ ['periodic'] ],
     
    random_seed = smilei_mpi_rank
)


Species(
	name = "ion",
	position_initialization = "regular",
	momentum_initialization = "cold",
	particles_per_cell = 10,
	mass = 1836.0,
	charge = 1.0,
	number_density = 1.,
	#time_frozen = 10000.0,
	boundary_conditions = [
		["periodic", "periodic"],
	],
	time_frozen = 0.1
)
Species(
	name = "eon1",
	position_initialization = "regular",
	momentum_initialization = "cold",
	particles_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	number_density = cosine(0.5,xamplitude=dn,xlength=L, xnumber=1),
	mean_velocity = [-0.1,0.0,0.0],
	boundary_conditions = [
		["periodic", "periodic"],
	],
)
Species(
	name = "eon2",
	position_initialization = "regular",
	momentum_initialization = "cold",
	particles_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	number_density = cosine(0.5,xamplitude=dn,xlength=L, xnumber=1),
	mean_velocity = [0.1,0.0,0.0],
	boundary_conditions = [
		["periodic", "periodic"],
	],
)


every = 100

DiagScalar(
    every = every,
)

DiagFields(
    every = every,
    fields = ['Ex','Ey','Ez','By_m','Bz_m','Rho']
)

DiagProbe(
    every = every,
    origin = [0],
    corners = [[0.99999*L]],
    number = [100],
    fields = ['Ex', 'Rho_eon1', 'Rho_eon2']
)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = every,
	species = ["eon1","eon2"],
	axes = [
		["x", 0., L, 50],
		["px", -0.4, 0.4, 100]
	]
)

DiagTrackParticles(
	species = "ion",
	every = 1000,
	filter = lambda particles: (particles.x<0.02),
	attributes = ["x", "px", "Ex"]
)

DiagPerformances(
    every = 300,
)
