
dx = 0.125
dt = 0.124
nx = 128*20
Lx = nx * dx
npatch_x = 256
laser_fwhm = 40.
laser_waist = 200.

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    timestep = dt,
    simulation_time = int(Lx/dt)*dt*10, #Laser propagates over 10 grid length
    cell_length  = [dx, 3.],
    grid_length = [ Lx,  1536],
    number_of_patches = [npatch_x, 64],
    clrw = nx/npatch_x,
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["silver-muller","silver-muller"],
    ],
    solve_poisson = False,
    print_every = 100,
    random_seed = smilei_mpi_rank
)

MovingWindow(
    time_start = Main.grid_length[0]*0.98,
    velocity_x = 0.9997
)

LoadBalancing(
    initial_balance = False,
    every = 20,
    cell_load = 1.,
)

Species( 
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 25,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = trapezoidal(0.003952, xvacuum=0., xplateau=1000000000, xslope1=785.),
    mean_velocity = [0.0, 0.0, 0.0],
    pusher = "boris",    
    time_frozen = 0.0,
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
    ],
)

LaserGaussian2D(
    box_side         = "xmin",
    a0              = 6.,
    focus           = [0., Main.grid_length[1]/2.],
    waist           = laser_waist,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)

Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Ey','Rho','Jx']

DiagFields(
    every = 500,
    fields = list_fields
)

DiagProbe(
    every = 100,
    origin = [0., Main.grid_length[1]/2.],
    corners = [
        [Main.grid_length[0], Main.grid_length[1]/2.],
    ],
    number = [nx],
    fields = ['Ex','Ey','Rho','Jx']
)

DiagScalar(
    every = 100,
    vars=[
        'Uelm','Ukin_electron',
        'ExMax','ExMaxCell','EyMax','EyMaxCell','RhoMin','RhoMinCell',
        'Ukin_bnd','Uelm_bnd','Ukin_out_mvw','Ukin_inj_mvw','Uelm_out_mvw','Uelm_inj_mvw'
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    species = ["electron"],
    axes = [
        ["moving_x", 0, Lx, 300],
        ["ekin", 1, 400, 100]
    ]
)
