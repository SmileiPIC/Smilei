# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
#
# Particle injection from the Xmin and Xmax boundaries
#
# ----------------------------------------------------------------------------------------

import math

# Mean velocity
mean_velocity = 0.999
# Electron temperature
Te = 0.01
# Ion temperature
Ti = 0.001
# Ion charge
Zi = 1
# Density
n0 = 1
# Debye length
Debye_length = 1. / np.sqrt( n0 / Te + Zi * n0 / Ti )
# Cell length
cell_length = [Debye_length*0.5, Debye_length*0.5, Debye_length*0.5]
# Number of patches
number_of_patches =[2, 4, 16]
# Cells per patches (patch shape)
cells_per_patch = [8., 8., 8.]
# Grid length
grid_length = [0.,0.,0.]
for i in range(3):
    grid_length[i] = number_of_patches[i] * cell_length[i] * cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 16
# Position init
position_initialization = 'random'
# Time step
timestep = 0.95/np.sqrt(1./ cell_length[0]**2 + 1./ cell_length[1]**2 + 1./ cell_length[2]**2)
# Total simulation time
simulation_time = ((0.5 - 0.125)*grid_length[2])/mean_velocity          # duration of the simulation
# Period of output for the diags
diag_every = int(simulation_time / timestep)

EM_boundary_conditions = [
    ['periodic'],
    ['periodic'],
    ['silver-muller']
]

boundary_conditions = [
    ["periodic", "periodic"],
    ["periodic", "periodic"],
    ["remove", "remove"],
]

mean_velocity_1 = [0.,0.,mean_velocity]
mean_velocity_2 = [0.,0.,-mean_velocity]

Main(
    geometry = "3Dcartesian",
    interpolation_order = 2 ,
    cell_length = cell_length,
    grid_length  = grid_length,
    number_of_patches = number_of_patches,
    timestep = timestep,
    simulation_time = simulation_time,
    EM_boundary_conditions = EM_boundary_conditions,
)

# Initial plasma shape
fp = trapezoidal(1., zvacuum=0.        ,zplateau=grid_length[2]/8.)
fm = trapezoidal(1., zvacuum=7*grid_length[2]/8.,zplateau=grid_length[2])

Species(
	name = 'pon1',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1836.0,
	charge = 1.0,
	number_density = fp,
	mean_velocity = mean_velocity_1,
	temperature = [Ti],
	time_frozen = 0.0,
	boundary_conditions = boundary_conditions,
)

ParticleInjector(
    species = 'pon1',
    box_side = 'zmin',
)

Species(
	name = 'eon1',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fp,
	mean_velocity = mean_velocity_1,
	temperature = [Te],
	time_frozen = 0.0,
	boundary_conditions = boundary_conditions,
)
ParticleInjector(
    species = 'eon1',
    box_side = 'zmin',
)

Species(
	name = 'pon2',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1836.0,
	charge = 1.0,
	number_density = fm,
	mean_velocity = mean_velocity_2,
	temperature = [Ti],
	time_frozen = 0.0,
	boundary_conditions = boundary_conditions,
)
ParticleInjector(
    species = 'pon2',
    box_side = 'zmax',
)

Species(
	name = 'eon2',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fm,
	mean_velocity = mean_velocity_2,
	temperature = [Te],
	time_frozen = 0.0,
	boundary_conditions = boundary_conditions,
)
ParticleInjector(
    species = 'eon2',
    box_side = 'zmax',
)

# Diags _______________________________________________________________

DiagScalar(every=1)

for species in ["eon1","pon1","eon2","pon2"]:
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = diag_every,
        time_average = 1,
        species = [species],
        axes = [
            ["x", 0, grid_length[0], int(grid_length[0]/cell_length[0])],
            ["y", 0, grid_length[1], int(grid_length[1]/cell_length[1])],
            ["z", 0, grid_length[2], int(grid_length[2]/cell_length[2])],
        ]
    )

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["gamma", 15., 30., 128],
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["gamma", 22.3, 22.45, 128],
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon2"],
    axes = [
        ["gamma", 15., 30., 128],
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon2"],
    axes = [
        ["gamma", 22.3, 22.45, 128],
    ]
)

# DiagFields(
#     every = globalEvery,
#     fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_pon1','Rho_eon1','Rho_pon2','Rho_eon2',"Jx","Jy","Jz"]
# )

# DiagProbe(
#     every = globalEvery,
#     origin = [0., Main.grid_length[1]/2.],
#     corners = [[Main.grid_length[0],Main.grid_length[1]/2.]],
#     number = [256],
#     fields = ['Ex','Ey','Ez','Rho_pon1','Jx_pon1','Jy_pon1','Jz_pon1']
# )
