# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
dx = 0.125
nx = 896
Lx = nx * dx
npatch_x = 128
laser_fwhm = 19.80

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2,
    
    timestep = 0.124,
    simulation_time = 100,
    
    cell_length  = [dx, 3.],
    grid_length = [ Lx,  120.],
    
    number_of_patches = [npatch_x, 4],
    
    clrw = nx/npatch_x,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    
    random_seed = 0,
    solve_poisson = False,
    print_every = 100
)

MovingWindow(
    time_start = Main.grid_length[0],
    velocity_x = 0.9997
)

LoadBalancing(
    initial_balance = False,
    every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

Species( 
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 50,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    number_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    time_frozen = 0.0,
    boundary_conditions = [
        ["remove", "remove"],
        ["stop", "stop"],
        ]
)

LaserGaussian2D(
    box_side         = "xmin",
    a0              = 2.,
    focus           = [0., Main.grid_length[1]/2.],
    waist           = 26.16/2,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)

DiagFields(
    every = 100000000,
    fields = ['Ex','Ey','Rho_electron','Rho_proton','Jx_electron']
)

DiagScalar(every = 100)

