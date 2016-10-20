
dx = 0.125
dt = 0.124
nx = 896
Lx = nx * dx
npatch_x = 128
laser_fwhm = 19.80

Main(
    geometry = "2d3v",
    
    interpolation_order = 2,
    
    timestep = dt,
    sim_time = int(2*Lx/dt)*dt,
    
    cell_length  = [dx, 3.],
    sim_length = [ Lx,  120.],
    
    number_of_patches = [npatch_x, 8],
    
    clrw = nx/npatch_x,
    
    bc_em_type_x = ["silver-muller","silver-muller"],
    bc_em_type_y = ["silver-muller","silver-muller"],
    
    random_seed = 0,
    solve_poisson = False,
    print_every = 100
)

MovingWindow(
    time_start = Main.sim_length[0],
    velocity_x = 0.9997
)

LoadBalancing(
    initial_balance = False,
    every = 20,
    coef_cell = 1.,
    coef_frozen = 0.1
)

Species( 
    species_type = "electron",
    initPosition_type = "regular",
    initMomentum_type = "cold",
    n_part_per_cell = 16,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    dynamics_type = "norm",    
    time_frozen = 0.0,
    radiating = False,
    bc_part_type_xmin = "supp",
    bc_part_type_xmax = "supp",
    bc_part_type_ymin ="supp",
    bc_part_type_ymax ="supp"
)

LaserGaussian2D(
    boxSide         = "xmin",
    a0              = 2.,
    focus           = [0., Main.sim_length[1]/2.],
    waist           = 26.16,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)

DumpRestart(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Ey','Rho_electron','Rho_proton','Jx_electron']

DiagFields(
    every = 100,
    fields = list_fields
)

DiagProbe(
    every = 10,
    pos = [0., Main.sim_length[1]/2.],
    pos_first = [Main.sim_length[0], Main.sim_length[1]/2.],
    number = [nx],
    fields = ['Ex','Ey','Rho','Jx']
)

DiagScalar(
    every = 10, 
    vars=['Uelm','Ukin_electron','ExMax','ExMaxCell','EyMax','EyMaxCell', 'RhoMax', 'RhoMaxCell']
)

