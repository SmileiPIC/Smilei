# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi       # laser wavelength
t0 = l0                # optical cycle

resx = 100.            # nb of cells in on laser wavelength
rest = 150.            # time of timestep in one optical cycle 

dx = l0/resx
dy = dx
dt = t0/rest

npatch_x = 8
npatch_y = 8

Lx = npatch_x*64*dx 
Ly = npatch_y*128*dy 


Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [dx,dy],
    grid_length  = [Lx,Ly],
    
    number_of_patches = [ npatch_x, npatch_y ],
    
    timestep = dt,
    simulation_time = 1500*dt,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['periodic'],
    ],
    cluster_width = 16,
)


# We build a gaussian laser from scratch instead of using LaserGaussian2D
# The goal is to test the space_time_profile attribute
omega = 1.
a0 = 150.
waist = 2.0*l0
Zr = omega * waist**2/2.
focus = [10.*l0, 5.*l0]

amplitude = a0 / math.sqrt(2)
w  = math.sqrt(1./(1.+(focus[0]/Zr)**2))
invWaist2 = (w/waist)**2
coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
def By(y,t):
    B = amplitude * w * math.exp( -invWaist2*(y-focus[1])**2 ) \
        * math.sin(omega*t - coeff*(y-focus[1])**2 + math.pi*0.5)
    if t < t0: return t/t0 * B
    else     : return B
def Bz(y,t):
    B = amplitude * w * math.exp( -invWaist2*(y-focus[1])**2 ) \
        * math.sin(omega*t - coeff*(y-focus[1])**2)
    if t < t0: return t/t0 * B
    else     : return B
Laser(
    box_side           = "xmin",
    space_time_profile = [By, Bz]
)


Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = 4,
    mass = 1836.0,
    charge = 1.0,
    number_density = trapezoidal(100.0,xvacuum=l0,xplateau=0.44*l0),
    boundary_conditions = [
        ["reflective", "reflective"],
        ["periodic", "periodic"],
    ],
)
Species(
    name = 'eon',
    position_initialization = 'regular',
    momentum_initialization = 'mj',
    particles_per_cell = 4,
    mass = 1.0,
    charge = -1.0,
    number_density = trapezoidal(100.0,xvacuum=l0,xplateau=0.44*l0),
    temperature = [0.001],
    boundary_conditions = [
        ["reflective", "reflective"],
        ["periodic", "periodic"],
    ], 
    time_frozen = 0.1
)

LoadBalancing(
    initial_balance = True,
    every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
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
	    point = [0., Ly/2.],
	    vector = [Lx/2., 0.1],
	    direction = direction,
	    deposited_quantity = "weight",
	    species = ["eon"],
	    axes = [["theta", -math.pi, math.pi, 10],],
	    every = 350
	)
	DiagScreen(
	    shape = "plane",
	    point = [Lx/2., Ly/2.],
	    vector = [1., 0.1],
	    direction = direction,
	    deposited_quantity = "weight",
	    species = ["eon"],
	    axes = [["a", -Ly/2., Ly/2., 10],],
	    every = 350
	)

DiagTrackParticles(
    species = "eon",
    every = 500,
)
