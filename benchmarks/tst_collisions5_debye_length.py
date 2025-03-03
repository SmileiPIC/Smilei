# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength


Main(
	geometry = "2Dcartesian",
	
	number_of_patches = [ 16, 16 ],
	
	interpolation_order = 2,
	
	timestep = 0.2 * L0,
	number_of_timesteps = 2,
	
	time_fields_frozen = 100000000000.,
	
	number_of_cells = [16*6, 16*6],
	grid_length = [L0, L0],
	
	EM_boundary_conditions = [ ["periodic"] ],
	
	reference_angular_frequency_SI = L0 * 3e8 /1.e-6,
	print_every = 1,
    
    patch_arrangement = "linearized_XY"
)

def n1(x,y):
    return 20*math.exp(-x/L0*10.)
def T1(x,y):
    return 0.001*math.exp(-y/L0*10.)

def n2(x,y):
    return 20*math.exp(-x/L0*20.)
def T2(x,y):
    return 0.001*math.exp(-y/L0*20.)


Species(
    name = "ions1",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = 9,
    mass = 1836.0,
    charge = 4.0,
    charge_density = n1,
    mean_velocity = [0., 0., 0.],
    temperature = [T1],
    time_frozen = 100000000.0,
    boundary_conditions = [
        ["periodic", "periodic"],
    ],
)

Species(
    name = "ions2",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = 9,
    mass = 7*1836.0,
    charge = 2.0,
    charge_density = n2,
    mean_velocity = [0., 0., 0.],
    temperature = [T2],
    time_frozen = 100000000.0,
    boundary_conditions = [
        ["periodic", "periodic"],
    ],
)

Species(
    name = "eons1",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell= 9,
    mass = 1.0,
    charge = -1.0,
    charge_density = n1,
    mean_velocity = [0., 0., 0.],
    temperature = [T1],
    time_frozen = 100000000.0,
    boundary_conditions = [
        ["periodic", "periodic"],
    ],
)

Species(
    name = "eons2",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell= 9,
    mass = 1.0,
    charge = -1.0,
    charge_density = n2,
    mean_velocity = [0., 0., 0.],
    temperature = [T2],
    time_frozen = 100000000.0,
    boundary_conditions = [
        ["periodic", "periodic"],
    ],
)

Collisions(
    species1 = ["eons1"],
    species2 = ["eons2"],
    coulomb_log = 0., # force calculation of Debye length
    debug_every = 1
)

for s in ["ions1", "ions2", "eons1", "eons2"]:
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = 4,
        species = [s],
        axes = [
                ["x",    0*L0,    Main.grid_length[0],   16],
                ["y",    0*L0,    Main.grid_length[1],   16],
        ]
    )
    DiagParticleBinning(
        deposited_quantity = "weight_charge",
        every = 4,
        species = [s],
        axes = [
                ["x",    0*L0,    Main.grid_length[0],   16],
                ["y",    0*L0,    Main.grid_length[1],   16],
        ]
    )
    DiagParticleBinning(
        deposited_quantity = "weight_ekin",
        every = 4,
        species = [s],
        axes = [
                ["x",    0*L0,    Main.grid_length[0],   16],
                ["y",    0*L0,    Main.grid_length[1],   16],
        ]
    )

