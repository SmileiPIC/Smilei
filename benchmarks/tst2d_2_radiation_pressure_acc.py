# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi       # laser wavelength
t0 = l0                # optical cycle
Lsim = [6.*l0,10.*l0]  # length of the simulation
Tsim = 10.*t0          # duration of the simulation
resx = 100.            # nb of cells in on laser wavelength
rest = 150.            # time of timestep in one optical cycle 

Main(
    geometry = "2d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    sim_length  = Lsim,
    
    number_of_patches = [ 8, 8 ],
    
    timestep = t0/rest,
    sim_time = Tsim,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['periodic'],
    ],
    
    random_seed = 0
)

LaserGaussian2D(
    box_side         = "xmin",
    a0              = 150.,
    focus           = [10.*l0, 5.*l0],
    waist           = 2.0*l0,
    ellipticity     = 1.,
    time_envelope   = ttrapezoidal(slope1=t0)
)


Species(
    species_type = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    n_part_per_cell = 4,
    mass = 1836.0,
    charge = 1.0,
    nb_density = trapezoidal(100.0,xvacuum=l0,xplateau=0.44*l0),
    boundary_conditions = [
        ["refl", "refl"],
        ["periodic", "periodic"],
    ],
)
Species(
    species_type = 'eon',
    position_initialization = 'regular',
    momentum_initialization = 'mj',
    n_part_per_cell = 4,
    mass = 1.0,
    charge = -1.0,
    nb_density = trapezoidal(100.0,xvacuum=l0,xplateau=0.44*l0),
    temperature = [0.001],
    boundary_conditions = [
        ["refl", "refl"],
        ["periodic", "periodic"],
    ], 
    time_frozen = 0.1
)



globalEvery = int(rest/2.)

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_ion','Rho_eon']
)


for direction in ["forward", "backward", "both", "canceling"]:
	DiagScreen(
	    shape = "sphere",
	    point = [0., Lsim[1]/2.],
	    vector = [Lsim[0]/2., 0.1],
	    direction = direction,
	    output = "density",
	    species = ["eon"],
	    axes = [["theta", -math.pi, math.pi, 10],],
	    every = 350
	)
	DiagScreen(
	    shape = "plane",
	    point = [Lsim[0]/2., Lsim[1]/2.],
	    vector = [1., 0.1],
	    direction = direction,
	    output = "density",
	    species = ["eon"],
	    axes = [["a", -Lsim[1]/2., Lsim[1]/2., 10],],
	    every = 350
	)

DiagTrackParticles(
    species = "eon",
    every = 500,
)

