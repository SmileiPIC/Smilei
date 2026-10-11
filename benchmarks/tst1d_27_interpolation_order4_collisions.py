# Regression for fourth-order 1D initialization with collision cell sorting.
Main(
    geometry = "1Dcartesian",
    interpolation_order = 4,
    cell_length = [0.1],
    grid_length = [6.4],
    number_of_patches = [4],
    timestep = 0.05,
    simulation_time = 0.1,
    EM_boundary_conditions = [["periodic"]],
    reference_angular_frequency_SI = 1.0e15,
    random_seed = 0,
    print_every = 1,
)

Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "mj",
    particles_per_cell = 4,
    mass = 1.0,
    charge = -1.0,
    number_density = 1.0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.001],
    pusher = "boris",
    boundary_conditions = [["periodic"]],
)

Species(
    name = "ion",
    position_initialization = "regular",
    momentum_initialization = "mj",
    particles_per_cell = 4,
    mass = 1836.0,
    charge = 1.0,
    number_density = 1.0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.001],
    pusher = "boris",
    boundary_conditions = [["periodic"]],
)

Collisions(
    species1 = ["electron"],
    species2 = ["ion"],
    coulomb_log = 3.0,
)

DiagScalar(every = 1)
