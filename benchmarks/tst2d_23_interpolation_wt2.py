import math

gam = 1000.0
Lx = 20.0
Ly = 8.0
Te = 0.01
dt_multiplier = 0.2

nx = 40
ny = 80
nppc = 4

simtime = 100.0
diagtime1 = 1.0

dx = Lx / nx
dy = Ly / ny
dt = min(dx, dy) * dt_multiplier
beta = math.sqrt(1.0- 1.0 / gam**2)

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    interpolation_WT = True,
    timestep = dt,
    simulation_time = simtime,
    cell_length = [dx, dy],
    grid_length  = [Lx, Ly],
    number_of_patches = [4, 8],
    print_every = 100,
)

Vectorization(
    mode = "off",
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
