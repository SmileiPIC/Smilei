import math

l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [7.*l0,10.*l0,10.*l0]  # length of the simulation
Tsim = 8.*t0                 # duration of the simulation
resx = 16.                    # nb of cells in one laser wavelength
rest = 30.                    # nb of timesteps in one optical cycle 

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4,4,4 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    
)

LaserGaussian3D(
    a0              = 1.,
    omega           = 1.,
    focus           = [0.9*Lsim[0], 0.6*Lsim[1], 0.3*Lsim[2]],
    waist           = l0,
    incidence_angle = [0.2, 0.1],
#    time_envelope   = tgaussian()
)


globalEvery = int(rest)

DiagScalar(
    every=globalEvery
)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez']
)
from numpy import s_
DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez'],
    subgrid = s_[4:100:3, 5:400:10, 6:300:80]
)
DiagFields(
    every = 2*globalEvery,
    fields = ['Ex','Ey','Ez'],
    time_average = globalEvery
)

DiagProbe(
    every = 10,
    origin = [0.1*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [30],
    origin = [0.1*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    corners = [[0.9*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [10, 10],
    origin = [0.1*Lsim[0], 0.*Lsim[1], 0.5*Lsim[2]],
    corners = [
        [0.9*Lsim[0], 0. *Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0.9*Lsim[1], 0.5*Lsim[2]],
    ],
    fields = []
)

DiagProbe(
    every = 100,
    number = [4, 4, 4],
    origin = [0.1*Lsim[0], 0.*Lsim[1], 0.5*Lsim[2]],
    corners = [
        [0.9*Lsim[0], 0. *Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0.9*Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0. *Lsim[1], 0.9*Lsim[2]],
    ],
    fields = []
)
