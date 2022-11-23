# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

dx = 0.1                 # longitudinal resolution
dy = 0.1
dz = dy                  # transverse resolution
dt = 0.9*dx/math.sqrt(3) # timestep
Lx = 64*dx
Ly = 64*dy
Lz = 64*dz

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [dx,dy,dz],
    grid_length  =[Lx,Ly,Lz],
    
    number_of_patches = [ 8, 8, 8 ],
    
    timestep = dt,
    simulation_time = 2000*dt,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['periodic'],
        ['periodic'],
    ],

    print_every = 100,
)

l0 = 2.*math.pi
thickness_plasma  = 0.2*l0

Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = 27,
    mass = 1836.0,
    charge = 1.0,
    number_density = trapezoidal(100.0,xvacuum=Lx/2.-thickness_plasma/2.,xplateau=thickness_plasma),
    boundary_conditions = [
        ["reflective", "reflective"],
        ["periodic", "periodic"],
        ["periodic", "periodic"],
    ],
)
Species(
    name = 'eon',
    position_initialization = 'regular',
    momentum_initialization = 'mj',
    particles_per_cell = 27,
    mass = 1.0,
    charge = -1.0,
    number_density = trapezoidal(100.0,xvacuum=Lx/2.-thickness_plasma/2.,xplateau=thickness_plasma),
    temperature = [0.05],
    boundary_conditions = [
        ["reflective", "reflective"],
        ["periodic", "periodic"],
        ["periodic", "periodic"],
    ], 
    
)

LoadBalancing(
    every = 100,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

Vectorization(
    mode = "adaptive",
    reconfigure_every = 20,
    initial_mode = "on"
)

DiagScalar(every=50)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["eon"],
    axes = [
        ["x", 0., Lx, 128],
        ["y", 0., Ly, 128],
	["z", 0., Lz, 128]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["ion"],
    axes = [
        ["x", 0., Lx, 128],
        ["y", 0., Ly, 128],
        ["z", 0., Lz, 128]
    ]
)

DiagScreen(
    #name = "my screen",
    shape = "plane",
    point = [0., 0., 0.],
    vector = [1., 0., 0.],
    direction = "canceling",
    deposited_quantity = "weight",
    species = ["eon"],
    axes = [["a", 0, Ly, 100],
            ["b", 0., Lz, 100]],
    every = 500
)

DiagScreen(
    #name = "my screen",
    shape = "plane",
    point = [0., 0., 0.],
    vector = [1., 0., 0.],
    direction = "canceling",
    deposited_quantity = "weight",
    species = ["ion"],
    axes = [["a", 0, Ly, 100],
            ["b", 0., Lz, 100]],
    every = 500
)


DiagFields(
    every = 500,
    fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_ion','Rho_eon','Rho']
)

