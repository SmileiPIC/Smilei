###### Namelist for the field initialization of an electron sphere in AM geometry

import math


dx = 1.
dtrans = 1.
dt = 0.8*dx
nx =  128 
ntrans = 128
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 8 
npatch_trans = 8 

# Sphere density
n0 = 1.

# Sphere position and Radius
R_sphere      = 20.
center_sphere = nx*dx/2. 


# normalized density of a sphere with uniform density
def nsphere_(x,y,z):
    r_sq = (y-Ltrans/2.)**2+(z-Ltrans/2.)**2
    if ( (x-center_sphere)**2+r_sq <  R_sphere**2 ):
    	return n0
    else:
        return 0.

Main(
    geometry = "3Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 1.*dt, 

    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],


    number_of_patches = [npatch_x,npatch_trans,npatch_trans ],
    
    clrw = nx/npatch_x,

    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["silver-muller","silver-muller"],
        ["silver-muller","silver-muller"],
    ],
    
    solve_poisson = True,
    poisson_max_iteration = 50000,    
    print_every = 100,

    random_seed = smilei_mpi_rank
)



Species(
    name = "electronSphere",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = nsphere_,
    mean_velocity = [0.0, 0.0, 0.0], 
    pusher = "boris",
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
       ["remove", "remove"],
       ["remove", "remove"],
    ],
)



list_fields = ['Ex','Ey','Ez','Rho']


DiagFields(
    every = 20,
        fields = list_fields
)

DiagProbe(
        every = 10,
        origin = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
        corners = [
            [Main.grid_length[0], Main.grid_length[1]/2., Main.grid_length[2]/2.]
        ],
        number = [nx],
        fields = list_fields
)

DiagProbe(
        every = 10,
        origin = [0., 0., Main.grid_length[2]/2.],
        corners = [
	    [Main.grid_length[0], 0., Main.grid_length[2]/2.],
            [0., Main.grid_length[1], Main.grid_length[2]/2.],
        ],
        number = [nx, ntrans],
        fields = list_fields
)


