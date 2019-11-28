###### Namelist for the field initialization of an electron sphere in AM geometry

import math


dx = 1.
dtrans = 1.
dt = 0.8*dx
nx =  128 
ntrans = 64 
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 8 
npatch_r = 8 

# Sphere density
n0 = 1.

# Sphere position and Radius
R_sphere      = 20.
center_sphere = nx*dx/2. 


# normalized density of a sphere with uniform density
def nsphere_(x,r):
    if ( (x-center_sphere)**2+r**2 <  R_sphere**2 ):
    	return n0
    else:
        return 0.

Main(
    geometry = "AMcylindrical",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 1.*dt, 

    cell_length  = [dx, dtrans],
    grid_length = [ Lx,  Ltrans],

    number_of_AM = 1,

    number_of_patches = [npatch_x,npatch_r ],
    
    clrw = nx/npatch_x,

    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],

    solve_poisson = False,
    
    solve_poisson = True,
    poisson_max_iteration = 50000,    
    print_every = 100,

    random_seed = smilei_mpi_rank
)



Species(
    name = "electronSphere",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 20,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = nsphere_,
    mean_velocity = [0.0, 0.0, 0.0], 
    pusher = "boris",
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
       ["reflective", "remove"],
    ],
)



list_fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho','Jx','Jy']

DiagFields(
    every = 200,
#        fields = list_fields
)


DiagProbe(
    every = 200,
    origin   = [0., -ntrans*dtrans,0.],
    corners  = [ [nx*dx,-ntrans*dtrans,0.], [0,ntrans*dtrans,0.] ],
    number   = [nx, 2*ntrans]
)

