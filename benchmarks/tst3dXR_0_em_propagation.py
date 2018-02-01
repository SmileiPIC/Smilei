# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Lsim = [20.,50.]  # length of the simulation
Tsim = t0           # duration of the simulation
resx = 28.              # nb of cells in on laser wavelength
rest = 60.              # time of timestep in one optical cycle 
laser_fwhm = 19.80
Main(
    geometry = "3drz",
    Nmode = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [l0/resx,20/resx],
    grid_length  = Lsim,
    number_of_patches = [ 1, 1 ],
    timestep = t0/rest,
    simulation_time = 10*t0/rest ,
     
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    
    random_seed = smilei_mpi_rank
)

LaserGaussian2D(
    a0              = 1.,
    omega           = 1.,
    focus           = [Lsim[0]/2., 0.],
    waist           = 25.,
    #time_envelope   = tconstant()
)
#LaserGaussian2D(
#    box_side         = "xmin",
#    a0              = 2.,
#    focus           = [0., 0.],
#    waist           = 26.16,
#    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
#)

globalEvery = int(1)


#DiagScalar(every=globalEvery)

DiagFields(
    every = 3,
    fields = []
)

#DiagProbe(
#    every = 100,
#    number = [100, 100],
#    pos = [0., 10.*l0],
#    pos_first = [20.*l0, 0.*l0],
#    pos_second = [3.*l0 , 40.*l0],
#    fields = []
#)
#
#DiagProbe(
#    every = 10,
#    pos = [0.1*Lsim[0], 0.5*Lsim[1]],
#    fields = []
#)

