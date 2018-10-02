################### 3D Laser Wakefield with envelope
dx = 1.
dtrans = 3.
dt = 0.8*dx
nx = 192
ntrans = 64
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 32
laser_fwhm = 20.
center_laser = Lx-2.*laser_fwhm # the temporal center here is the same as waist position, but in principle they can differ
time_start_moving_window =  0.


Main(
    geometry = "3Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 350.*dt,

    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],

    number_of_patches =[npatch_x, 8, 8],

    clrw = nx/npatch_x,

    EM_boundary_conditions = [ ["silver-muller"] ],
    Envelope_boundary_conditions = [ ["reflective", "reflective"],
        ["reflective", "reflective"],
        ["reflective", "reflective"], ],

    solve_poisson = False,
    print_every = 100,

    random_seed = smilei_mpi_rank
)

Vectorization(
    mode = "normal",
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
       ["remove", "remove"],
       ["remove", "remove"],
    ],
)

LaserEnvelopeGaussian3D( # linear regime of LWFA
    a0              = 0.1,
    focus           = [center_laser, Main.grid_length[1]/2.,Main.grid_length[2]/2.],
    waist           = 30.,
    time_envelope   = tgaussian(center=center_laser, fwhm=laser_fwhm),
    envelope_solver = 'explicit',
)


Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Ey','Rho','Jx','Env_A_abs','Env_Chi']

DiagFields(
   every = 50,
        fields = list_fields
)

DiagProbe(
        every = 50,
        origin = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
        corners = [
            [Main.grid_length[0], Main.grid_length[1]/2., Main.grid_length[2]/2.]
        ],
        number = [nx],
        fields = ['Ex','Ey','Rho','Jx','Env_A_abs','Env_Chi']
)

#DiagProbe(
#        every = 10,
#        origin = [0., Main.grid_length[1]/4., Main.grid_length[2]/2.],
#        corners = [
#            [Main.grid_length[0], Main.grid_length[1]/4., Main.grid_length[2]/2.],
#            [0., 3*Main.grid_length[1]/4., Main.grid_length[2]/2.],
#        ],
#        number = [nx, ntrans],
#        fields = ['Ex','Ey','Rho','Jx']
#)

#DiagScalar(every = 10, vars=['Uelm','Ukin_electron','ExMax','ExMaxCell','EyMax','EyMaxCell', 'RhoMin', 'RhoMinCell'])
DiagScalar(every = 10, vars=['Env_A_absMax'])

#DiagParticleBinning(
#       deposited_quantity = "weight_charge",
#       every = 50,
#       species = ["electron"],
#       axes = [
#               ["moving_x", 0, Main.grid_length[0], nx],
#               ["px", -1, 2., 100]
#       ]
#)
