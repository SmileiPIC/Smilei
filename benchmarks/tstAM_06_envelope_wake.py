############################# Laser envelope propagation in vacuum
dx = 0.69 
dr = 5. 
dt = 0.57#0.8*dx
nx = 512
nr = 60
Lx = nx * dx
Lr = nr*dr
npatch_x=32
laser_fwhm = 70.7 
center_laser = 2*laser_fwhm # here is the same as focus position of laser but in principle they can differ
time_start_moving_window = Lx-4*laser_fwhm


Main(
    geometry = "AMcylindrical",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 1700.*dt,

    cell_length  = [dx, dr],
    grid_length = [ Lx,  Lr],

    number_of_AM = 1,

    number_of_patches = [npatch_x,4],
    clrw = nx/npatch_x,

    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],

    solve_poisson = False,
    print_every = 100,

    random_seed = smilei_mpi_rank
)

MovingWindow(
    time_start = time_start_moving_window,
    velocity_x = 1.0
)

LoadBalancing(
    initial_balance = False,
        every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

LaserEnvelopeGaussianAM(
    a0              = 0.01,
    focus           = [center_laser, Main.grid_length[1]/2.],
    waist           = 94.26,
    time_envelope   = tgaussian(center=center_laser, fwhm=laser_fwhm),
    envelope_solver = 'explicit',
    Envelope_boundary_conditions = [ ["reflective", "reflective"],
        ["reflective", "reflective"], ],
)


n0 = 0.0017
Radius_plasma = 300.
longitudinal_profile = polygonal(xpoints=[3.*laser_fwhm,3.5*laser_fwhm,24.5*laser_fwhm,25.*laser_fwhm],xvalues=[0.,n0,n0,0.])
def nplasma(x,r):
    profile_r = 0.
    if ((r)**2<Radius_plasma**2):
        profile_r = 1.
    return profile_r*longitudinal_profile(x,r)



Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 4,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = nplasma,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.,0.,0.],
    pusher = "ponderomotive_boris",
    ponderomotive_dynamics = "True",
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
       ["remove", "remove"],
    ],
)


Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)


DiagFields(
    every = 100,
    fields = ['Env_A_abs_mode_0','Env_Chi_mode_0','Env_E_abs_mode_0' ]
)

DiagFields(
    every = 100,
    fields = ['El_mode_0' ]
)

DiagProbe(
        every = 50,
        origin = [0., 2.*dr, 2.*dr],
        corners = [
            [Main.grid_length[0], 2.*dr, 2.*dr]
        ],
        number = [nx],
        fields = ['Ex','Ey','Rho','Jx','Rho_electron','Jx_electron','Env_A_abs','Env_Chi','Env_E_abs']
)


DiagProbe(
    every = 50,
    origin   = [0., -nr*dr,0.],
    corners  = [ [nx*dx,-nr*dr,0.], [0,nr*dr,0.] ],
    number   = [nx, 2*nr],
    fields = ['Ex','Ey','Rho','Jx','Rho_electron','Jx_electron','Env_A_abs','Env_Chi','Env_E_abs']
)



#DiagScalar(every = 10, vars=['Env_A_absMax','Env_E_absMax'])


                                                                                                                                                                 

