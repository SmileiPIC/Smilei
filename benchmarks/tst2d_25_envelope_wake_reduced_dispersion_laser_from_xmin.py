################### 2D Laser Wakefield with envelope
dx = 1. 
dtrans = 3.
dt = 0.8*dx
nx = 192
ntrans = 32 
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 32
laser_fwhm = 20. 
center_laser = 2.*laser_fwhm 
time_start_moving_window = Lx


Main(
    geometry = "2Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 1501.*dt,

    cell_length  = [dx, dtrans],
    grid_length = [ Lx,  Ltrans],

    number_of_patches =[npatch_x, 4],
    
    cluster_width = nx/npatch_x,

    EM_boundary_conditions = [ ["silver-muller","silver-muller"],["PML","PML"] ],
    number_of_pml_cells = [[0,0],[15,15]],
    

    solve_poisson = False,
    print_every = 100,
    use_BTIS3_interpolation = True,

)

MovingWindow(
    time_start = time_start_moving_window,
    velocity_x = 1. 
)

LoadBalancing(
    initial_balance = False,
        every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

start_plasma = 5*dx
Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 1,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = polygonal(xpoints=[start_plasma+2.*laser_fwhm,start_plasma+2.1*laser_fwhm,15000,20000],xvalues=[0.,0.0045,0.0045,0.]),
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    pusher = "ponderomotive_borisBTIS3", # pusher to interact with envelope
    #pusher = "boris",
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
       ["remove", "remove"],
    ],
)

LaserEnvelopeGaussian2D(
    box_side        = "xmin",
    a0              = 0.1,  # linear regime of LWFA  
    focus           = [center_laser, Main.grid_length[1]/2.],
    waist           = 30.,
    time_envelope   = tgaussian(center=center_laser, fwhm=laser_fwhm),
    envelope_solver = 'explicit_reduced_dispersion',
    Envelope_boundary_conditions = [ ["reflective", "reflective"], ["PML", "PML"], ],
)

Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Ey','Rho','Jx','Env_A_abs','Env_Chi','Bz_mBTIS3']

DiagFields(
   every = 50,
        fields = list_fields
)

DiagProbe(
        every = 50,
        origin = [0., Main.grid_length[1]/2.],
        corners = [
            [Main.grid_length[0], Main.grid_length[1]/2.]
        ],
        number = [nx],
        fields = ['Ex','Ey','Rho','Env_A_abs','Env_Chi','Env_E_abs','BzBTIS3']
)


DiagScalar(every = 10, vars=['Env_A_absMax','Env_E_absMax'])

