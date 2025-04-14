from math import pi, cos, sin

l0 = 2.0*pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [7.*l0,10.*l0,20.*l0]  # length of the simulation
Tsim = 12.*t0                 # duration of the simulation
resx = 16.                    # nb of cells in one laser wavelength
rest = 30.                    # nb of timesteps in one optical cycle 

angle1 = 0.2*pi
angle2 = 0.07*pi

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4,4,4 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    EM_boundary_conditions_k = [
        [cos(angle1)*cos(angle2), sin(angle2), -sin(angle1)*cos(angle2)],
        [-cos(angle1)*cos(angle2), -sin(angle2), sin(angle1)*cos(angle2)],
        [0., 1., 0.],[0., -1., 0.],
        [0., 0., 1.],[0., 0., -1.],],
)

LaserGaussian3D(
    a0              = 1.,
    omega           = 1.,
    focus           = [0.5*Lsim[0], 0.6*Lsim[1], 0.3*Lsim[2]],
    waist           = 2*l0,
    incidence_angle = [angle1, angle2],
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
