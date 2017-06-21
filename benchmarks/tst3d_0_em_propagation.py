import math

l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [7.*l0,10.*l0,10.*l0]  # length of the simulation
Tsim = 8.*t0                 # duration of the simulation
resx = 16.                    # nb of cells in one laser wavelength
rest = 30.                    # nb of timesteps in one optical cycle 

Main(
    geometry = "3d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx,l0/resx],
    sim_length  = Lsim,
    
    number_of_patches = [ 4,4,4 ],
    
    timestep = t0/rest,
    sim_time = Tsim,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    
    random_seed = 0
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

DiagProbe(
    every = 10,
    pos = [0.1*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [30],
    pos = [0.1*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    pos_first = [0.9*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [10, 10],
    pos = [0.1*Lsim[0], 0.*Lsim[1], 0.5*Lsim[2]],
    pos_first  = [0.9*Lsim[0], 0. *Lsim[1], 0.5*Lsim[2]],
    pos_second = [0.1*Lsim[0], 0.9*Lsim[1], 0.5*Lsim[2]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [4, 4, 4],
    pos = [0.1*Lsim[0], 0.*Lsim[1], 0.5*Lsim[2]],
    pos_first  = [0.9*Lsim[0], 0. *Lsim[1], 0.5*Lsim[2]],
    pos_second = [0.1*Lsim[0], 0.9*Lsim[1], 0.5*Lsim[2]],
    pos_third  = [0.1*Lsim[0], 0. *Lsim[1], 0.9*Lsim[2]],
    fields = []
)
