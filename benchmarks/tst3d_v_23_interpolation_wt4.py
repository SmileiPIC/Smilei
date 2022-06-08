import math

gam = 1000.0
Lx = 8.0
Ly = 7.0
Lz = 6.0
Te = 0.01
dt_multiplier = 0.2

nx = 20
ny = 32
nz = 40
nppc = 8

simtime = 10.0
diagtime1 = 1.0

dx = Lx / nx
dy = Ly / ny
dz = Lz / nz
dt = min(dx, dy, dz) * dt_multiplier
beta = math.sqrt(1.0- 1.0 / gam**2)

Main(
    geometry = "3Dcartesian",
    interpolation_order = 4,
    interpolator = "wt",
    timestep = dt,
    simulation_time = simtime,
    cell_length = [dx, dy, dz],
    grid_length  = [Lx, Ly, Lz],
    number_of_patches = [2, 2, 4],
    print_every = 50,
)

Vectorization(
    mode = "on",
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
