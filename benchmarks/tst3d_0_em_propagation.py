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
    
    bc_em_type_x = ['silver-muller'],
    bc_em_type_y = ['silver-muller'],
    bc_em_type_z = ['silver-muller'],
    
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
