###### Namelist for the field initialization of a relativistic electron bunch, with cylindrical uniform density

import math
dx = 0.8 
dtrans = 4.
dt = 0.75
nx = 320
ntrans = 48
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 32

###################################

# Electron bunch parameters
center_bunch = nx*dx/2.           # centered at half the grid length
n0 = 0.0017                       # plasma density
alpha = 1.                        # bunch density normalized to plasma density
gamma= 100.                       # relativistic lorentz factor of the bunch, moving in the positive x direction
beta = math.sqrt(1.-1/gamma**2)   # bunch velocity, normalized to c
L_bunch = 30.                     # bunch length
R_bunch = 20.                     # bunch radius

# Normalized density of a cylindrical uniformly charged bunch
def nbunch_uniform_cylindrical(x,y,z):
        bunch_density = alpha*n0
        if ( (abs(x-center_bunch)<L_bunch/2.) and   ( (y-Main.grid_length[1]/2.)**2+(z-Main.grid_length[2]/2.)**2) < R_bunch**2 ):
                return bunch_density
        else:
                return 0.

###################################


Main(
    geometry = "3Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = dt*2., 

    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],

    number_of_patches = [npatch_x, 8, 8],
    
    clrw = nx/npatch_x,

    EM_boundary_conditions = [ ["silver-muller"] ],

    solve_poisson = False,
    
    solve_relativistic_poisson = True,
    time_fields_frozen = 100.,
   
    print_every = 100,

    random_seed = smilei_mpi_rank
)


Species(
    name = "bunch_electrons",
    position_initialization = "regular",
    momentum_initialization = "cold",
    relativistic_field_initialization = True,
    particles_per_cell = 1,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = nbunch_uniform_cylindrical,
    mean_velocity = [beta, 0.0, 0.0], 
    pusher = "boris",
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
       ["remove", "remove"],
       ["remove", "remove"],
    ],
)


list_fields = ['Ex','Ey','Ez','Bx','By','Bz']

DiagFields(
    every = 1,
        fields = list_fields
)

DiagProbe(
        every = 1,
        origin = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
        corners = [
            [Main.grid_length[0], Main.grid_length[1]/2., Main.grid_length[2]/2.]
        ],
        number = [nx],
        fields = ['Ex','Ey','Rho','Jx']
)
                                                                                                                                                                 
