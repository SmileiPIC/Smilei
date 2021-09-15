# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
from numpy import exp, sin, cos, tan, sqrt, pi

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Lsim = [40.*l0,10.*l0]  # length of the simulation
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
offset_y = 15.
offset_x = Lsim[0]*0.7
FWHM_x = l0
FWHM_t = 50.

LaserOffset(
    box_side = "ymin",
    space_time_profile = [None, lambda x,t: exp( -2.77259*((x-offset_x)/FWHM_x)**2 - (2.77259*(t-80.)/FWHM_t)**2 ) * sin(t)],
    offset = offset_y - (Lsim[0]-offset_x)*sin(angle),
    extra_envelope = lambda x,t: 0. if t < 20. else 1.,
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
