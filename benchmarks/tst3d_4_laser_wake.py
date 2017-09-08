
dx = 0.125
dtrans = 3.
dt = 0.124
nx = 896
ntrans = 40
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 128
laser_fwhm = 19.80

Main(
    geometry = "3d3v",
    
    interpolation_order = 2,
    
    timestep = dt,
    sim_time = int(2*Lx/dt)*dt,
    
    cell_length  = [dx, dtrans, dtrans],
    sim_length = [ Lx,  Ltrans, Ltrans],
    
    number_of_patches = [npatch_x, 4, 4],
    
    clrw = nx/npatch_x,
    
    EM_boundary_conditions = [ ["silver-muller"] ],
    
    solve_poisson = False,
    print_every = 100,

    random_seed = smilei_mpi_rank
)

MovingWindow(
    time_start = Main.sim_length[0],
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
    position_initialization = "regular",
    momentum_initialization = "cold",
    n_part_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    pusher = "boris",    
    time_frozen = 0.0,
    radiating = False,
    boundary_conditions = [
    	["supp", "supp"],
    	["supp", "supp"],
    	["supp", "supp"],
    ],
)

LaserGaussian3D(
    box_side         = "xmin",
    a0              = 2.,
    focus           = [0., Main.sim_length[1]/2., Main.sim_length[2]/2.],
    waist           = 26.16,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)

Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Ey','Rho','Jx']

DiagFields(
    every = 100,
    fields = list_fields
)

DiagProbe(
	every = 10,
	origin = [0., Main.sim_length[1]/2., Main.sim_length[2]/2.],
	corners = [
	    [Main.sim_length[0], Main.sim_length[1]/2., Main.sim_length[2]/2.]
	],
	number = [nx],
	fields = ['Ex','Ey','Rho','Jx']
)

DiagProbe(
	every = 10,
	origin = [0., Main.sim_length[1]/4., Main.sim_length[2]/2.],
	corners = [
	    [0., 3*Main.sim_length[1]/4., Main.sim_length[2]/2.],
	    [Main.sim_length[0], Main.sim_length[1]/4., Main.sim_length[2]/2.]
	],
	number = [nx, ntrans],
	fields = ['Ex','Ey','Rho','Jx']
)

DiagScalar(every = 10, vars=['Uelm','Ukin_electron','ExMax','ExMaxCell','EyMax','EyMaxCell', 'RhoMin', 'RhoMinCell'])

DiagParticleBinning(
	output = "charge_density",
	every = 50,
	species = ["electron"],
	axes = [
		["moving_x", 0, Main.sim_length[0], nx],
		["px", -1, 2., 100]
	]
)

