# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Tsim = t0           # duration of the simulation
resx = 28.              # nb of cells in on laser wavelength
rest = 60.              # time of timestep in one optical cycle 
laser_fwhm = 19.80
dx=0.2
nx=240
dt=0.18
Lsim=[dx*nx,50.] # length of the simulation

Main(
    geometry = "3drz",
    nmodes = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [dx,20/resx],
    grid_length  = Lsim,
    number_of_patches = [ 8, 1 ],
    timestep = dt,
    simulation_time = 1000*t0/rest ,
     
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    
    random_seed = smilei_mpi_rank
)

MovingWindow(
    time_start = Main.grid_length[0]-40*dx,
    velocity_x = 0.9997
)

LaserGaussian2D(
    a0              = 1.,
    omega           = 1.,
    focus           = [Lsim[0]/2., 0.],
    waist           = 8.,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
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
    every = 10,
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

