
dx = 0.2
dtrans = 3.
dt = 0.19
nx = 512
ntrans = 40
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 64
laser_fwhm = 19.80

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2,

    timestep = dt,
    simulation_time = int(2*Lx/dt)*dt,

    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],

    number_of_patches = [npatch_x, 4, 4],
    uncoupled_grids = "True",

    clrw = nx/npatch_x,
    
    EM_boundary_conditions = [ ["silver-muller"] ],
    EM_boundary_conditions_k = [ [1., 0., 0.],[-1., 0., 0.],[1., 0.005, 0.],[1., -0.005, 0.],[1., 0., 0.005],[1., 0., -0.005] ],
    
    solve_poisson = False,
    print_every = 100,

    random_seed = smilei_mpi_rank
)

MovingWindow(
    time_start = Main.grid_length[0],
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
    particles_per_cell = 1,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    pusher = "boris",
    time_frozen = 0.0,
    boundary_conditions = [
    	["remove", "remove"],
    	["remove", "remove"],
    	["remove", "remove"],
    ],
)


# We build a gaussian laser from scratch instead of using LaserGaussian3D
# The goal is to test the space_time_profile attribute
omega = 1.
a0 = 2.
focus = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.]
waist = 10.
time_envelope = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)

Zr = omega * waist**2/2.
w  = math.sqrt(1./(1.+(focus[0]/Zr)**2))
invWaist2 = (w/waist)**2
coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
def By(y,z,t):
    return 0.
def Bz(y,z,t):
    r2 = (y-focus[1])**2 + (z-focus[2])**2
    omegat = omega*t - coeff*r2
    return a0 * w * math.exp( -invWaist2*r2  ) * time_envelope( omegat/omega ) * math.sin( omegat )
Laser(
    box_side = "xmin",
    space_time_profile = [By, Bz]
)

#LaserGaussian3D(
#    box_side         = "xmin",
#    a0              = 2.,
#    focus           = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
#    waist           = 10.,
#    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
#)



Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Ey','Rho','Jx']

#DiagFields(
#    every = 100,
#    fields = list_fields
#)

DiagProbe(
	every = 10,
	origin = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
	corners = [
	    [Main.grid_length[0], Main.grid_length[1]/2., Main.grid_length[2]/2.]
	],
	number = [nx],
	fields = list_fields+["Jx_electron"]
)

DiagProbe(
	every = 40,
	origin = [0., Main.grid_length[1]/4., Main.grid_length[2]/2.],
	corners = [
	    [Main.grid_length[0], Main.grid_length[1]/4., Main.grid_length[2]/2.],
	    [0., 3*Main.grid_length[1]/4., Main.grid_length[2]/2.],
	],
	number = [nx, ntrans],
	fields = list_fields
)

DiagScalar(every = 10, vars=['Uelm','Ukin_electron','ExMax','ExMaxCell','EyMax','EyMaxCell', 'RhoMin', 'RhoMinCell'])

DiagParticleBinning(
	deposited_quantity = "weight_charge",
	every = 50,
	species = ["electron"],
	axes = [
		["moving_x", 0, Main.grid_length[0], nx],
		["px", -1, 2., 100]
	]
)
