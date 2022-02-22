import math

gam = 1000.0
Lx = 20.0
Te = 0.01
dt_multiplier = 0.2

nx = 160
nppc = 8

simtime = 100.0
diagtime1 = 1.0

dx = Lx / nx
dt = dx * dt_multiplier
beta = math.sqrt(1.0- 1.0 / gam**2)

Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    interpolator = "wt",
    timestep = dt,
    simulation_time = simtime,
    cell_length = [dx],
    grid_length  = [Lx],
    number_of_patches = [16],
    print_every = 100,
)

Species(
    name = "electon",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = nppc,
    mass = 1.0, 
    charge = -1.0,
    number_density = gam,
    mean_velocity = [beta, 0.0, 0.0],
    temperature = [Te],
)

Species(
    name = "positron",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = nppc,
    mass = 1.0, 
    charge = 1.0,
    number_density = gam,
    mean_velocity = [beta, 0.0, 0.0],
    temperature = [Te],
)

DiagScalar(
    every = [int(diagtime1/dt)],
    vars = ['Utot', 'Ukin', 'Uelm', 'Uexp', 'Ubal', 'Ubal_norm'],
)
