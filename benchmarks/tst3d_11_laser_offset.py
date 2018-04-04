import math

from numpy import exp, sin
l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [4.*l0,10.*l0,10.*l0]  # length of the simulation
Tsim = 5.*t0                 # duration of the simulation
resx = 16.                    # nb of cells in one laser wavelength
rest = 30.                    # nb of timesteps in one optical cycle 

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 2,4,4 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    
    random_seed = smilei_mpi_rank
)


LaserOffset(
    space_time_profile = [None, lambda y,z,t: exp( -(40.*(y/Lsim[1]-0.5))**2-(40.*(z/Lsim[2]-0.5))**2 - (1.*(t-Tsim/5.)/Tsim)**2 ) * sin(t)],
    offset = 10.,
    time_envelope = 1.,
    keep_n_strongest_modes = 60
)

DiagFields(
    every = 10,
    fields = ['Ex','Ey','Ez','Bx','By','Bz']
)
