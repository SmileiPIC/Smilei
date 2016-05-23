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
    geometry = "2d3v",
    
    interpolation_order = 2,
    
    timestep = 0.124,
    sim_time = 200,
    
    cell_length  = [dx, 3.],
    sim_length = [ Lx,  120.],
    
    number_of_patches = [npatch_x, 8],
    
    clrw = nx/npatch_x,
    
    bc_em_type_x = ["silver-muller","silver-muller"],
    bc_em_type_y = ["silver-muller","silver-muller"],
    
    random_seed = 0,
    
    print_every = 100
)

MovingWindow(
    time_start = Main.sim_length[0],
    velocity_x = 0.9997
)

LoadBalancing(
    every = 20,
    coef_cell = 1.,
    coef_frozen = 0.1
)


Species(
    species_type = "proton",
    initPosition_type = "regular",
    initMomentum_type = "cold",
    ionization_model = "none",
    n_part_per_cell = 10, 
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    dynamics_type = "norm",
    time_frozen = 100000.,
    radiating = False,
    bc_part_type_west  = "supp",
    bc_part_type_east  = "supp",
    bc_part_type_south = "supp",
    bc_part_type_north = "supp"
)

Species( 
    species_type = "electron",
    initPosition_type = "regular",
    initMomentum_type = "cold",
    n_part_per_cell = 10,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    dynamics_type = "norm",    
    time_frozen = 0.0,
    radiating = False,
    bc_part_type_west = "supp",
    bc_part_type_east = "supp",
    bc_part_type_south ="stop",
    bc_part_type_north ="stop"
)

LaserGaussian2D(
    boxSide         = "west",
    a0              = 2.,
    focus           = [0., Main.sim_length[1]/2.],
    waist           = 26.16,
    time_envelope   = tgaussian(center=2*laser_fwhm, fwhm=laser_fwhm)
)

DumpRestart(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

DiagFields(
    every = 100,
    fields = ['Ex','Ey','Rho_electron','Rho_proton','Jx_electron']
)

DiagScalar(every = 100)

