################### 1D Laser Wakefield with envelope
dx = 1. 
dtrans = 3.
dt = 0.8*dx
nx = 192
Lx = nx * dx
npatch_x = 32
laser_fwhm = 20. 
center_laser = Lx-2.*laser_fwhm # the temporal center here is the same as waist position, but in principle they can differ
time_start_moving_window =  0.


Main(
    geometry = "1Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 350.*dt,

    cell_length  = [dx],
    grid_length = [ Lx],

    number_of_patches =[npatch_x],
    
    clrw = nx/npatch_x,

    EM_boundary_conditions = [ ["silver-muller"] ],
   

    solve_poisson = False,
    print_every = 100,

    random_seed = smilei_mpi_rank
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

Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 1,
    c_part_max = 1.0,
    ponderomotive_dynamics = True, # = this species interacts with laser envelope
    mass = 1.0,
    charge = -1.0,
    charge_density = polygonal(xpoints=[center_laser+2.*laser_fwhm,center_laser+2.1*laser_fwhm,15000,20000],xvalues=[0.,0.0045,0.0045,0.]),
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    pusher = "ponderomotive_boris", # pusher to interact with envelope
    #pusher = "boris",
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
    ],
)

LaserEnvelopePlanar1D( # linear regime of LWFA
    a0              = 0.1,     
    time_envelope   = tgaussian(center=center_laser, fwhm=laser_fwhm),
    envelope_solver = 'explicit',
     Envelope_boundary_conditions = [ ["reflective", "reflective"],
     ],
)


Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Rho','Env_A_abs','Env_Chi','Env_E_abs']

DiagFields(
   every = 50,
        fields = list_fields
)


DiagScalar(every = 10, vars=['Env_A_absMax','Env_E_absMax'])

                                                                                                                                                               

