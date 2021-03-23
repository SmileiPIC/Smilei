import math

l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [10.*l0,7.*l0,10.*l0]  # length of the simulation
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
    
    random_seed = smilei_mpi_rank
)

LaserGaussian3D(
    box_side        = "ymin",
    a0              = 1.,
    omega           = 1.,
    focus           = [0.6*Lsim[0], 0.9*Lsim[1], 0.3*Lsim[2]],
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
