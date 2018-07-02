import math
l0 = 2.0*math.pi  # wavelength in normalized units
t0 = l0           # optical cycle in normalized units
rest = 102.0      # nb of timestep in 1 optical cycle
resx = 100.0      # nb cells in 1 wavelength

Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    
    cell_length = [l0/resx],
    grid_length  = [6.0*l0],
    
    number_of_patches = [ 8 ],
    
    timestep = t0/rest,
    simulation_time = 16.0*t0,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    
    random_seed = smilei_mpi_rank,
    
    print_every = int(rest/2.0)
)

MovingWindow(
    time_start = 12.*t0,
    velocity_x = 0.9997
)

Laser(
    omega          = 1.,
    chirp_profile  = tpolynomial(order3=5e-6),
    time_envelope  = tgaussian(fwhm=1.*t0),
    space_envelope = [1., 0.],
)

DiagScalar(
    every = 5
)

DiagFields(
    every = int(rest/2.0),
    fields = ['Ex','Ey','Ez','By_m','Bz_m']
)
from numpy import s_
DiagFields(
    every = int(rest/2.0),
    fields = ['Ex','Ey','Ez','By_m','Bz_m'],
    subgrid = s_[4:100:3]
)
DiagFields(
    every = int(rest*2.),
    fields = ['Ex','Ey','Ez','By_m','Bz_m'],
    time_average = int(rest)
)

DiagProbe(
    every = 5, 
    origin = [Main.grid_length[0]*0.2]   
)

DiagProbe(
    every = 5,
    origin = [0.0],
    corners = [Main.grid_length],
    number = [1000]
)
