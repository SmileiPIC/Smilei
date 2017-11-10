# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi                 # laser wavelength [in code units]
t0 = l0                          # optical cycle
Lsim = [32.*l0,64.*l0]           # length of the simulation
Tsim = 98.*t0                    # duration of the simulation
resx = 16.                       # nb of cells in on laser wavelength
rest = resx*mat.sqrt(2.)/0.95   # time of timestep in one optical cycle 

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 8, 8 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
    ],
    
    random_seed = smilei_mpi_rank
)

LaserGaussian2D(
    a0              = 1.,
    omega           = 1.,
    focus           = [Lsim[0]/2., Lsim[1]/2.],
    waist           = 5.*l0,
    incidence_angle = 0.,
    time_envelope   = tgaussian(fwhm=4*t0)
)


globalEvery = int(rest)

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez','Bx','By','Bz']
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


