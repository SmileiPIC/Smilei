############################# Laser envelope propagation in vacuum
dx = 1. 
dr = 3. 
dt = 0.8*dx
nx = 128
nr = 64
Lx = nx * dx
Lr = nr*dr
npatch_x=16
laser_fwhm = 20. 
center_laser = 2*laser_fwhm # here is the same as waist position of laser but in principle they can differ
time_start_moving_window = 0.


Main(
    geometry = "AMcylindrical",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 300.*dt,

    cell_length  = [dx, dr],
    grid_length = [ Lx,  Lr],

    number_of_AM = 1,

    number_of_patches = [npatch_x,8],
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
    a0              = 1.,
    focus           = [center_laser, Main.grid_length[1]/2.],
    waist           = 30.,
    time_envelope   = tgaussian(center=center_laser, fwhm=laser_fwhm),
    envelope_solver = 'explicit',
    Envelope_boundary_conditions = [ ["reflective", "reflective"],
        ["reflective", "reflective"], ],
)

Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)


#DiagFields(
#    every = 100,
#)

DiagProbe(
        every = 50,
        origin = [0., 2.*dr, 2.*dr],
        corners = [
            [Main.grid_length[0], 2.*dr, 2.*dr]
        ],
        number = [nx],
        fields = ['Ex','Ey','Rho','Jx','Env_A_abs','Env_Chi','Env_E_abs']
)


DiagProbe(
    every = 200,
    origin   = [0., -nr*dr,0.],
    corners  = [ [nx*dx,-nr*dr,0.], [0,nr*dr,0.] ],
    number   = [nx, 2*nr],
    fields = ['Ex','Ey','Rho','Jx','Env_A_abs','Env_Chi','Env_E_abs']
)


#DiagScalar(every = 10, vars=['Env_A_absMax','Env_E_absMax'])


                                                                                                                                                                 

