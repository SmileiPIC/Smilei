# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

from math import pi, cos, sin

l0 = 2.0*pi        # laser wavelength
t0 = l0                 # optical cycle
Lsim = [40.*l0,20.*l0]  # length of the simulation
Tsim = 50.*t0           # duration of the simulation
resx = 28.              # nb of cells in on laser wavelength
rest = 40.              # time of timestep in one optical cycle 

ang = 0.5

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4, 8 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
    ],
    
    EM_boundary_conditions_k = [[1.,0.],[-1.,0.],[cos(ang), sin(ang)],[0.,-1.]],
    
    random_seed = smilei_mpi_rank
)

LaserGaussian2D(
    box_side        = "ymin",
    a0              = 1.,
    omega           = 1.,
    focus           = [Lsim[0]*0.8, Lsim[1]],
    waist           = 8.,
    incidence_angle = ang,
    polarization_phi = 20./180.*pi,
    time_envelope   = tgaussian()
)


globalEvery = int(rest)

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez']
)


