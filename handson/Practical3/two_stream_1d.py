# ----------------------------------------------------------------------------------------
#                  SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

# -------------------
# MY PYTHON VARIABLES
# -------------------

output_every=25
n_part=100
velocity=0.1
amplitude=0.001
length = 1.03

# --------------------------------------
# SMILEI's VARIABLES (DEFINED IN BLOCKS)
# --------------------------------------

# MAIN
Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    simulation_time = 10.0*2*math.pi,
    timestep = 0.01,
    grid_length  = [length],
    cell_length = [length/64],
    number_of_patches = [ 4 ],
    EM_boundary_conditions = [ ['periodic'] ]
)

Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = n_part,
    mass = 1836.0,
    charge = 1.0,
    number_density = 1.0,
    time_frozen = 10000.0,
    boundary_conditions = [
        ['periodic'],
    ]
)

Species(
    name = 'eon1',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = n_part/2,
    mass = 1.0,
    charge = -1.0,
    number_density = cosine(0.5,xamplitude=amplitude,xnumber=1),
    mean_velocity = [velocity, 0, 0],
    boundary_conditions = [
        ['periodic'],
    ]
)

Species(
    name = 'eon2',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = n_part/2,
    mass = 1.0,
    charge = -1.0,
    number_density = cosine(0.5,xamplitude=amplitude,xnumber=1),
    mean_velocity = [-velocity, 0, 0],
    boundary_conditions = [
        ['periodic'],
    ]
)


DiagScalar (
    precision = 3,
    every=output_every,
    vars = ['Utot', 'Ukin', 'Uelm', 'Ukin_eon1', 'Ukin_eon2', 'Uelm_Ex']
)
 
DiagParticleBinning(
    deposited_quantity = "weight",
    every = output_every,
    species = ['eon1', 'eon2'],
    axes = [
        ["x", 0., length, 100],
        ["px", -4*velocity, 4*velocity, 100]
    ]
)

DiagFields(
    every = output_every,
)


