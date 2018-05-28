# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Tsim = t0           # duration of the simulation
laser_fwhm = 19.80
dx = 0.2
dr = 1.5
nx = 800
nr = 100
Lsim = [dx*nx,nr*dr] # length of the simulation
dt = 0.18

Main(
    geometry = "3drz",
    nmodes = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [dx, dr],
    grid_length  = Lsim,
    number_of_patches = [ 32, 4 ],
    timestep = dt,
    simulation_time = 2*nx*dt,
     
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    
    random_seed = smilei_mpi_rank,
    print_every = 1,
)

MovingWindow(
    time_start = Main.grid_length[0]-50*dx,
    velocity_x = 0.9997
)

LaserGaussian2D(
    a0              = 1.,
    omega           = 1.,
    #focus           = [Lsim[0]/2., 0.],
    focus           = [0. , 0.],
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
    every = 100,
    fields = ["Br_m_mode_0", "Br_m_mode_1","Bx_m_mode_0","Bx_m_mode_1","Bt_m_mode_0","Bt_m_mode_1","Bt_mode_0","Bt_mode_1","Bx_mode_0","Bx_mode_1","Br_mode_0","Br_mode_1","Er_mode_0","Er_mode_1","Et_mode_0","Et_mode_1","Ex_mode_0","Ex_mode_1" ]
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

