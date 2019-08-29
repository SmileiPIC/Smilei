# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math


Lsim = [32.,8.]    # length of the simulation
Tsim = 64          # duration of the simulation
resx = 8.          # nb of cells in on laser wavelength
rest = 16.         # time of timestep in one optical cycle
nppc = 16

mean_V = 0.4
pos_all = 'random'

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2 ,
    cell_length = [1./resx,1./resx],
    grid_length  = Lsim,
    number_of_patches = [ 4, 4 ],
    timestep = 1./rest,
    simulation_time = Tsim,
    EM_boundary_conditions = [
        ['silver-muller'],
        ['periodic'],
    ],
    random_seed = smilei_mpi_rank
)

fp = trapezoidal(1., xvacuum=0.        ,xplateau=Lsim[0]/2.)
fm = trapezoidal(1., xvacuum=Lsim[0]/2.,xplateau=Lsim[0]/2.)

Species(
	name = 'pon1',
	position_initialization = pos_all,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = nppc,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	number_density = fp,
	mean_velocity = [mean_V,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

ParticleInjector(
    species = 'pon1',
    box_side = 'xmin',
)

Species(
	name = 'eon1',
	position_initialization = pos_all,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = nppc,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fp,
	mean_velocity = [mean_V,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)
ParticleInjector(
    species = 'eon1',
    box_side = 'xmin',
)

Species(
	name = 'pon2',
	position_initialization = pos_all,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = nppc,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	number_density = fm,
	mean_velocity = [-mean_V,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)
ParticleInjector(
    species = 'pon2',
    box_side = 'xmax',
)

Species(
	name = 'eon2',
	position_initialization = pos_all,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = nppc,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fm,
	mean_velocity = [-mean_V,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)
ParticleInjector(
    species = 'eon2',
    box_side = 'xmax',
)


globalEvery = rest/2

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_pon1','Rho_eon1','Rho_pon2','Rho_eon2',"Jx","Jy","Jz"]
)

DiagProbe(
    every = globalEvery,
    origin = [0., Main.grid_length[1]/2.],
    corners = [[Main.grid_length[0],Main.grid_length[1]/2.]],
    number = [256],
    fields = ['Ex','Ey','Ez','Rho_pon1','Jx_pon1','Jy_pon1','Jz_pon1']
)

