############################# Laser envelope propagation in vacuum
dx = 1. 
dtrans = 3. 
dt = 0.8*dx
nx = 128
Lx = nx * dx
npatch_x=16
laser_fwhm = 20. 
center_laser = 2*laser_fwhm # here is the same as waist position of laser but in principle they can differ
time_start_moving_window = 0.


Main(
    geometry = "1Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 300.*dt,

    cell_length  = [dx],
    grid_length = [ Lx],

    number_of_patches = [npatch_x],
    clrw = nx/npatch_x,

    EM_boundary_conditions = [ ["silver-muller"] ],

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

LaserEnvelopePlanar1D(
    a0              = 1.,
    time_envelope   = tgaussian(center=center_laser, fwhm=laser_fwhm),
    envelope_solver = 'explicit_reduced_dispersion',
    Envelope_boundary_conditions = [ ["reflective", "reflective"],
         ],
)

Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Rho','Jx','Env_A_abs']

DiagFields(
    every = 100,
        fields = list_fields
)

# DiagProbe(
#         every = 50,
#         origin = [0., Main.grid_length[1]/2.],
#         corners = [
#             [Main.grid_length[0], Main.grid_length[1]/2.]
#         ],
#         number = [nx],
#         fields = ['Ex','Rho','Jx','Env_A_abs','Env_Chi']
# )



#DiagScalar(every = 10, vars=['Uelm','Ukin_electron','ExMax','ExMaxCell','EyMax','EyMaxCell', 'RhoMin', 'RhoMinCell'])
DiagScalar(every = 10, vars=['Env_A_absMax'])


                                                                                                                                                                 

