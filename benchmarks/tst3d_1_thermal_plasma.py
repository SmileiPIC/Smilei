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
    
    timestep = dt,
    sim_time = Tsim,
    
    cell_length  = [dx,dy,dz],
    sim_length = [Lx,Ly,Lz],
    
    number_of_patches = [4,4,4],
    
    bc_em_type_x = ["periodic"],
    bc_em_type_y = ["periodic"],
    bc_em_type_z = ["periodic"],
    
    random_seed = 0,
    
    print_every = 1
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
    bc_part_type_west  = "none",
    bc_part_type_east  = "none",
    bc_part_type_south = "none",
    bc_part_type_north = "none",
    bc_part_type_bottom = "none",
    bc_part_type_up = "none"
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
    bc_part_type_west  = "none",
    bc_part_type_east  = "none",
    bc_part_type_south = "none",
    bc_part_type_north = "none",
    bc_part_type_bottom = "none",
    bc_part_type_up = "none"
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

