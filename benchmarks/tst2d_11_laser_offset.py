# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
from numpy import exp, sin

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Lsim = [10.*l0,50.*l0]  # length of the simulation
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

LaserOffset(
    space_time_profile = [None, lambda y,t: exp( -(600.*(y/Lsim[1]-0.5))**2 - (10.*(t-40.-Tsim/5.)/Tsim)**2 ) * sin(t)],
    offset = 40.,
    time_envelope = tpolygonal(points=[t0, Tsim/2.], values=[1., 1.]),
    keep_n_strongest_modes = 800
)

DiagFields(
    every = 25,
    fields = ['Ex','Ey','Ez','Bx','By','Bz']
)

