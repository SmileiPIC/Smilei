# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi		# laser wavelength
t0 = l0					# optical cicle
Lsim = [20.*l0,20.*l0]	# length of the simulation
Tsim = 20.*t0			# duration of the simulation
resx = 5.				# nb of cells in on laser wavelength
rest = 8.				# time of timestep in one optical cycle 

Main(
    geometry = "2d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    sim_length  = Lsim,
    
    number_of_patches = [ 4, 4  ],
    
    timestep = t0/rest,
    sim_time = Tsim,
     
    bc_em_type_x = ['silver-muller'],
    bc_em_type_y = ['silver-muller'],
    
    random_seed = 0
)

globalEvery = 5

Species(
    species_type = "electron",
	initPosition_type = 'random',
	initMomentum_type = 'cold',
	n_part_per_cell = 2,
	mass = 1.0,
	charge = -1.0,
	nb_density = trapezoidal(1.0,xvacuum=1.*l0,xplateau=4.*l0,yvacuum=5.*l0,yplateau=10.*l0),
	bc_part_type_xmin  = 'refl',
	bc_part_type_xmax  = 'refl',
	bc_part_type_ymin = 'none',
	bc_part_type_ymax = 'none',
	mean_velocity=[0.9,0.01,0]
)

Species(
	initPosition_type = 'random',
	initMomentum_type = 'cold',
	n_part_per_cell = 2,
	mass = 1.0,
	charge = 1.0,
	nb_density = trapezoidal(1.0,xvacuum=1.*l0,xplateau=4.*l0,yvacuum=5.*l0,yplateau=10.*l0),
	bc_part_type_xmin  = 'refl',
	bc_part_type_xmax  = 'refl',
	bc_part_type_ymin = 'none',
	bc_part_type_ymax = 'none',
	mean_velocity=[0.9,0.01,0]
)

Antenna(
    field='Jz',
    time_profile= lambda t: math.sin(t/t0),
    space_profile=gaussian(0.2, xfwhm=l0, yfwhm=l0, xcenter=Main.sim_length[0]*0.6, ycenter=Main.sim_length[1]*0.5)
)

PartWall (
    x= 15.0*l0,
    kind="refl"
)

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ez','Jz','Rho_electron','Rho_species1']
)

DiagScreen(
    shape = "plane",
    point = [10.*l0, 10.*l0],
    vector = [1., 0.],
    output = "density",
    species = ["electron"],
    axes = [["a", -10.*l0, 10.*l0, 40],
            ["px", 0., 3., 30]],
    every = 10
)
