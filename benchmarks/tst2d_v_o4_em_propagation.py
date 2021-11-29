# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

from math import pi, cos, sin

l0 = 2.0*pi        # laser wavelength
t0 = l0                 # optical cycle
Lsim = [20.*l0,50.*l0]  # length of the simulation
Tsim = 50.*t0           # duration of the simulation
resx = 28.              # nb of cells in on laser wavelength
rest = 40.              # time of timestep in one optical cycle

ang = 0.5

Main(
    geometry = "2Dcartesian",

    interpolation_order = 4 ,

    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,

    number_of_patches = [ 8, 4 ],

    timestep = t0/rest,
    simulation_time = Tsim,

    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
    ],

    EM_boundary_conditions_k = [[cos(ang), sin(ang)],[-1.,0.],[0.,1.],[0.,-1.]],

    random_seed = smilei_mpi_rank
)

Vectorization(
    mode="on",
)

LaserGaussian2D(
    a0              = 1.,
    omega           = 1.,
    focus           = [Lsim[0], Lsim[1]/2.],
    waist           = 8.,
    incidence_angle = ang,
    polarization_phi = 20./180.*pi,
    time_envelope   = tgaussian()
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
    every = globalEvery,
    fields = ['Ex','Ey','Ez']
)
from numpy import s_
DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez'],
    subgrid = s_[4:100:3, 5:400:10]
)
DiagFields(
    every = 2*globalEvery,
    fields = ['Ex','Ey','Ez'],
    time_average = globalEvery
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

DiagProbe(
    every = 100,
    number = [10, 1000],
    origin = [0., 0.],
    corners = [
        [Lsim[0], 0.],
        [0., Lsim[1]],
    ],
    fields = ["PoyX", "PoyY", "PoyZ"],
    time_integral = True
)


DiagTrackParticles(
    species = "eon",
    every = 4*rest,
	attributes = ["x", "y", "px", "py", "pz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"]
)
