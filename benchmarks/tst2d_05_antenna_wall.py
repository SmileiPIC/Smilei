# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi		# laser wavelength
t0 = l0					# optical cicle
Lsim = [20.*l0,20.*l0]	# length of the simulation
Tsim = 20.*t0			# duration of the simulation
resx = 10.				# nb of cells in on laser wavelength
rest = 16.				# time of timestep in one optical cycle 

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4, 4  ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
    ],
    
)

globalEvery = 5

Species(
    name = "electron",
	position_initialization = 'random',
	momentum_initialization = 'cold',
	particles_per_cell = 2,
	mass = 1.0,
	charge = -1.0,
	number_density = trapezoidal(1.0,xvacuum=1.*l0,xplateau=4.*l0,yvacuum=5.*l0,yplateau=10.*l0),
	boundary_conditions = [
		["reflective", "reflective"],
		["periodic", "periodic"],
	],
	mean_velocity=[0.9,0.01,0]
)

Species(
	position_initialization = 'random',
	momentum_initialization = 'cold',
	particles_per_cell = 2,
	mass = 1.0,
	charge = 1.0,
	number_density = trapezoidal(1.0,xvacuum=1.*l0,xplateau=4.*l0,yvacuum=5.*l0,yplateau=10.*l0),
	boundary_conditions = [
		["reflective", "reflective"],
		["periodic", "periodic"],
	],
	mean_velocity=[0.9,0.01,0]
)

Antenna(
    field='Jz',
    time_profile= lambda t: math.sin(2.*t/t0),
    space_profile=gaussian(0.2, xfwhm=l0, yfwhm=l0, xcenter=Main.grid_length[0]*0.6, ycenter=Main.grid_length[1]*0.5)
)

PartWall (
    x= 15.0*l0,
    kind="reflective"
)

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ez','Jz','Rho_electron','Rho_species1']
)

DiagScreen(
    shape = "plane",
    point = [13.*l0, 10.*l0],
    vector = [1., 0.],
    direction = "canceling",
    deposited_quantity = "weight",
    species = ["electron"],
    axes = [["a", -10.*l0, 10.*l0, 40],
            ["p", 0., 3., 30]],
    every = 10
)

DiagScreen(
    shape = "sphere",
    point = [5.*l0, 10.*l0],
    vector = [5.*l0, 0.],
    direction = "both",
    deposited_quantity = "weight",
    species = ["electron"],
    axes = [["theta", -math.pi, math.pi, 40],
            ["p", 0., 3., 30]],
    every = 10
)
