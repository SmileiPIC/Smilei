import math

from numpy import exp, sin, pi
l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [4.*l0,10.*l0,10.*l0]  # length of the simulation
Tsim = 8.*t0                 # duration of the simulation
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
    
)


angle = 20./180.*pi
offset_x = 4.
offset_y = Lsim[1]*0.6
offset_z = Lsim[2]*0.5

LaserOffset(
    space_time_profile = [None, lambda y,z,t: exp( -(20.*(y-offset_y)/Lsim[1])**2-(20.*(z-offset_z)/Lsim[2])**2 - (3.*(t-30.)/Tsim)**2 ) * sin(t)],
    offset = 4.,
    keep_n_strongest_modes = 100,
    angle = angle
)

DiagFields(
    every = 10,
    fields = ['Ex','Ey','Ez','Bx','By','Bz']
)
