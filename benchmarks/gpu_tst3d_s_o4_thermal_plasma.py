# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math as m


TkeV = 10.						# electron & ion temperature in keV
T   = TkeV/511.   				# electron & ion temperature in me c^2
n0  = 1.
Lde = m.sqrt(T)					# Debye length in units of c/\omega_{pe}
dx  = 0.5*Lde 					# cell length (same in x & y)
dy  = dx
dz  = dx
dt  = 0.95 * dx/m.sqrt(3.)		# timestep (0.95 x CFL)

Lx    = 40.*dx
Ly    = 40.*dy
Lz    = 40.*dz
Tsim  = 2.*m.pi

def n0_(x,y,z):
	if (0.1*Lx<x<0.9*Lx) and (0.1*Ly<y<0.9*Ly) and (0.1*Lz<z<0.9*Lz):
		return n0
	else:
		return 0.


Main(gpu_computing = True,
    geometry = "3Dcartesian",
    
    interpolation_order = 4,
    
    timestep = dt,
    simulation_time = Tsim,
    
    cell_length  = [dx,dy,dz],
    grid_length = [Lx,Ly,Lz],
    
    number_of_patches = [4,4,4],
    
    EM_boundary_conditions = [ ["periodic"] ],
    
    print_every = 1,

)


LoadBalancing(
    every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

Vectorization(
    mode = "off",
)

Species(
    name = "proton",
    position_initialization = "regular",
    momentum_initialization = "mj",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    charge_density = n0_,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [T],
    pusher = "boris",
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
)
Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "mj",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [T],
    pusher = "boris",
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
)

Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

DiagFields(
    every = 4
)

DiagScalar(every = 1)

for direction in ["forward", "backward", "both", "canceling"]:
	DiagScreen(
	    shape = "sphere",
	    point = [0., Ly/2., Lz/2.],
	    vector = [Lx*0.9, 0.1, 0.1],
	    direction = direction,
	    deposited_quantity = "weight",
	    species = ["electron"],
	    axes = [
	    	["theta", 0, math.pi, 10],
	    	["phi", -math.pi, math.pi, 10],
	    	],
	    every = 40,
	    time_average = 30
	)
	DiagScreen(
	    shape = "plane",
	    point = [Lx*0.9, Ly/2., Lz/2.],
	    vector = [1., 0.1, 0.1],
	    direction = direction,
	    deposited_quantity = "weight",
	    species = ["electron"],
	    axes = [
	    	["a", -Ly/2., Ly/2., 10],
	    	["b", -Lz/2., Lz/2., 10],
	    	],
	    every = 40,
	    time_average = 30
	)
