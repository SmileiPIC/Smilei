# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cicle
Lsim = [40.*l0,5.*l0]   # length of the simulation
Tsim = 30.*t0           # duration of the simulation
resx = 50.              # nb of cells in on laser wavelength
rest = 75.              # time of timestep in one optical cycle 

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 16, 2 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['periodic'],
    ],
    
    random_seed = smilei_mpi_rank
)


fp = trapezoidal(1., xvacuum=l0, xplateau=10.*l0)
fm = trapezoidal(1., xvacuum=11.*l0, xplateau=28.5*l0)

Species(
	name = 'pon1',
	position_initialization = 'regular',
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = 4,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	number_density = fp,
	mean_velocity = [0.,0.,0.],
	temperature = [0.001],
	time_frozen = 100000000.0,
	boundary_conditions = [
		["stop", "reflective"],
		["periodic", "periodic"],
	],
)
Species(
	name = 'eon1',
	position_initialization = 'regular',
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = 4,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fp,
	mean_velocity = [0.,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	boundary_conditions = [
		["stop", "reflective"],
		["periodic", "periodic"],
	],
)
Species(
	name = 'pon2',
	position_initialization = 'regular',
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = 4,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	number_density = fm,
	mean_velocity = [-0.5,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	boundary_conditions = [
		["stop", "reflective"],
		["periodic", "periodic"],
	],
)
Species(
	name = 'eon2',
	position_initialization = 'regular',
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = 4,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fm,
	mean_velocity = [-0.5,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	boundary_conditions = [
		["stop", "reflective"],
		["periodic", "periodic"],
	],
)



every = 50
globalEvery = int(rest/2.)

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_pon1','Rho_eon1','Rho_pon2','Rho_eon2']
)

