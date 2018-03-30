# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
from numpy import exp, sin

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Lsim = [20.*l0,50.*l0]  # length of the simulation
Tsim = 50.*t0           # duration of the simulation
resx = 28.              # nb of cells in on laser wavelength
rest = 40.              # time of timestep in one optical cycle 

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 8, 4 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
    ],
    
    random_seed = smilei_mpi_rank
)

#LaserGaussian2D(
#    a0              = 1.,
#    omega           = 1.,
#    focus           = [Lsim[0], Lsim[1]/2.],
#    waist           = 8.,
#    incidence_angle = 0.5,
#    time_envelope   = tgaussian()
#)


LaserOffset(
    space_time_profile = [None, lambda y,t: exp( -(400.*(y/Lsim[1]-0.5))**2 - (10.*(t-40.-Tsim/5.)/Tsim)**2 ) * sin(t)],
    offset = 40.,
    time_envelope = tpolygonal(points=[t0, Tsim/2.], values=[1., 1.]),
    keep_n_best_frequencies = 800
)

Species(
	name = "eon",
	position_initialization = "regular",
	momentum_initialization = "cold",
	particles_per_cell = 0.001,
	mass = 1.0,
	charge = -1.0,
	number_density = 1.,
	mean_velocity = [0.,0.,0.],
	boundary_conditions = [
		["periodic", "periodic"],
	],
	is_test = True
)


globalEvery = int(rest)

DiagScalar(every=globalEvery)

DiagFields(
    every = 25,
    fields = ['Ex','Ey','Ez','Bx','By','Bz']
)
from numpy import s_
DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez'],
    subgrid = s_[4:100:3, 5:400:10]
)

DiagProbe(
    every = 100,
    number = [100, 100],
    origin = [0., 10.*l0],
    corners = [
        [20.*l0, 0.*l0],
        [3.*l0 , 40.*l0],
    ],
    fields = []
)

DiagProbe(
    every = 10,
    origin = [0.1*Lsim[0], 0.5*Lsim[1]],
    fields = []
)


DiagTrackParticles(
    species = "eon",
    every = 4*rest,
	attributes = ["x", "y", "Ex", "Ey", "Ez", "Bx", "By", "Bz"]
)

