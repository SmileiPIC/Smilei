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
Lx    = 32.*dx
Ly    = 32.*dy
Lz    = 32.*dz
Tsim  = 2.*m.pi			

def n0_(x,y,z):
	if (0.1*Lx<x<0.9*Lx) and (0.1*Ly<y<0.9*Ly) and (0.1*Lz<z<0.9*Lz):
		return n0
	else:
		return 0.


Main(
    geometry = "3d3v",
    interpolation_order = 2,
    global_factor = [2,2,2],
    norder = [2,2,2],
    is_spectral = True,
    is_pxr = True,
    timestep = dt,
    sim_time = Tsim,
    cell_length  = [dx,dy,dz],
    sim_length = [Lx,Ly,Lz],
    number_of_patches = [4,4,4],
    bc_em_type_x = ["periodic"],
    bc_em_type_y = ["periodic"],
    bc_em_type_z = ["periodic"],
    print_every = 1,
    random_seed = smilei_mpi_rank
)


LoadBalancing(
    every = 20,
    coef_cell = 1.,
    coef_frozen = 0.1
)


Species(
    species_type = "proton",
    initPosition_type = "regular",
    initMomentum_type = "mj",
    n_part_per_cell = 8, 
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    charge_density = n0_,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [T],
    dynamics_type = "norm",
    bc_part_type_xmin  = "none",
    bc_part_type_xmax  = "none",
    bc_part_type_ymin = "none",
    bc_part_type_ymax = "none",
    bc_part_type_zmin = "none",
    bc_part_type_zmax = "none"
)
Species(
    species_type = "electron",
    initPosition_type = "regular",
    initMomentum_type = "mj",
    n_part_per_cell = 8, 
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [T],
    dynamics_type = "norm",
    bc_part_type_xmin  = "none",
    bc_part_type_xmax  = "none",
    bc_part_type_ymin = "none",
    bc_part_type_ymax = "none",
    bc_part_type_zmin = "none",
    bc_part_type_zmax = "none"
)

DumpRestart(
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
	    output = "density",
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
	    output = "density",
	    species = ["electron"],
	    axes = [
	    	["a", -Ly/2., Ly/2., 10],
	    	["b", -Lz/2., Lz/2., 10],
	    	],
	    every = 40,
	    time_average = 30
	)


