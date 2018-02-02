# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi	# laser wavelength
t0 = l0			# optical cicle
Lsim = [20.*l0,20.*l0]	# length of the simulation
Tsim = 30.*t0		# duration of the simulation
resx = 64.		# nb of cells in on laser wavelength
rest = 100.		# time of timestep in one optical cycle 

Main(
    geometry = "2d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    sim_length  = Lsim,
    
    number_of_patches = [ 64, 64 ],
    clrw = 20,
    
    timestep = t0/rest,
    sim_time = Tsim/3.,
     
    bc_em_type_x = ['silver-muller'],
    bc_em_type_y = ['periodic'],
    
    random_seed = 0
)

LoadBalancing(
    every = 20
)

fp = trapezoidal(1., xvacuum=l0, xplateau=5.*l0)
fm = trapezoidal(1., xvacuum=5.5*l0, xplateau=14.25*l0)

Species(
	species_type = 'pon1',
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	ionization_model = 'none',
	n_part_per_cell = 36,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	nb_density = fp,
	mean_velocity = [0.,0.,0.],
	temperature = [0.001],
	time_frozen = 100000000.0,
	bc_part_type_xmin  = 'stop',
	bc_part_type_xmax  = 'refl',
	bc_part_type_ymin = 'none',
	bc_part_type_ymax = 'none'
)
Species(
	species_type = 'eon1',
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	ionization_model = 'none',
	n_part_per_cell = 36,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	nb_density = fp,
	mean_velocity = [0.,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	bc_part_type_xmin  = 'stop',
	bc_part_type_xmax  = 'refl',
	bc_part_type_ymin = 'none',
	bc_part_type_ymax = 'none'
)
Species(
	species_type = 'pon2',
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	ionization_model = 'none',
	n_part_per_cell = 36,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	nb_density = fm,
	mean_velocity = [-0.5,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	bc_part_type_xmin  = 'stop',
	bc_part_type_xmax  = 'refl',
	bc_part_type_ymin = 'none',
	bc_part_type_ymax = 'none'
)
Species(
	species_type = 'eon2',
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	ionization_model = 'none',
	n_part_per_cell = 36,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	nb_density = fm,
	mean_velocity = [-0.5,0.,0.],
	temperature = [0.001],
	time_frozen = 0.0,
	bc_part_type_xmin  = 'stop',
	bc_part_type_xmax  = 'refl',
	bc_part_type_ymin = 'none',
	bc_part_type_ymax = 'none'
)



every = 50
globalEvery = int(rest/2.)

DiagScalar(every=globalEvery)

DiagFields(
    every = 1000,
    fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_pon1','Rho_eon1','Rho_pon2','Rho_eon2']
)

