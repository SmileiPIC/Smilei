# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
from numpy import exp, sin, cos, tan, sqrt, pi

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Lsim = [10.*l0,40.*l0]  # length of the simulation
Tsim = 25.*t0           # duration of the simulation
resx = 28.              # nb of cells in on laser wavelength
rest = 40.              # time of timestep in one optical cycle 

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4, 4 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
    ],
    
    random_seed = smilei_mpi_rank
)

angle = 25./180.*pi
offset_x = 15.
offset_y = Lsim[1]*0.7
FWHM_y = l0
FWHM_t = 50.

LaserOffset(
    space_time_profile = [None, lambda y,t: exp( -2.77259*((y-offset_y)/FWHM_y)**2 - (2.77259*(t-80.)/FWHM_t)**2 ) * sin(t)],
    offset = offset_x - (Lsim[1]-offset_y)*sin(angle),
    extra_envelope = lambda y,t: 0. if t < 20. else 1.,
    keep_n_strongest_modes = 10000,
    angle = angle
)

DiagFields(
    every = 25,
    fields = ['Ex','Ey','Ez','Bx','By','Bz']
)

DiagScalar(
    every = 10,
)
