# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

dx = 0.2
dr = 1.5
nx = 800
nr = 100
Lsim = [dx*nx,nr*dr] # length of the simulation
dt = 0.18

Main(
    geometry = "AMcylindrical",
    nmodes = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [dx, dr],
    grid_length  = Lsim,
    number_of_patches = [ 32, 4 ],
    timestep = dt,
    simulation_time = 2.5*nx*dt,
     
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    
    random_seed = smilei_mpi_rank,
    print_every = 100,
)

MovingWindow(
    time_start = Main.grid_length[0]-50*dx,
    velocity_x = 0.9997
)

laser_fwhm = 19.80
LaserGaussian2D(
    box_side         = "xmin",
    a0              = 2.,
    focus           = [0., 0.],
    waist           = 25.,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)

DiagFields(
    every = 100,
    fields = ["Br_m_mode_0", "Br_m_mode_1","Bx_m_mode_0","Bx_m_mode_1","Bt_m_mode_0","Bt_m_mode_1","Bt_mode_0","Bt_mode_1","Bx_mode_0","Bx_mode_1","Br_mode_0","Br_mode_1","Er_mode_0","Er_mode_1","Et_mode_0","Et_mode_1","Ex_mode_0","Ex_mode_1" ]
)

#DiagProbe(
#    every = 10,
#    origin = [1., 10., 0.],
#    fields = []
#)
#DiagProbe(
#    every = 10,
#    origin = [0., 10., 0.],
#    corners = [[Lsim[0], 10., 0.]],
#    number=[100],
#    fields = []
#)
#DiagProbe(
#    every = 10,
#    origin = [0., -10., 0.],
#    corners = [[Lsim[0], -10., 0.]],
#    number=[100],
#    fields = []
#)


