###### Namelist for plasma wake excited by a relativistic electron bunch

import math
dx = 0.2 
dtrans = 1.
dt = 0.18
#nx = 896
nx = 1280
ntrans = 80
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 128
bunch_sigma_x = 15.
bunch_sigma_r = 10. 
center_bunch = 6*bunch_sigma_x
n0 = 0.0017
alpha = 1.
gamma= 100. # relativistic lorentz factor
beta = math.sqrt(1.-1/gamma**2)
relative_energy_spread = 0.01
electron_mass_eV = 0.511e6
norm_emittance_m = 1.e-6 # transverse normalized emittance
L_bunch = 2.*bunch_sigma_x   
R_bunch = bunch_sigma_r




# normalized density of a cylindrical bunch, uniformly charged
def nbunch_uniform_cylindrical(x,y,z):

        bunch_density = alpha*n0

        if ( (abs(x-center_bunch)<L_bunch/2.) and   ( (y-Main.grid_length[1]/2.)**2+(z-Main.grid_length[2]/2.)**2) < R_bunch**2 ):
                return bunch_density
        else:
                return 0.





Main(
    geometry = "3Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = dt*2., 

    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],

    number_of_patches = [npatch_x, 4, 4],
    
    cluster_width = nx/npatch_x,

    EM_boundary_conditions = [ ["silver-muller"] ],

    solve_poisson = False,
    
    solve_relativistic_poisson = True,
    
    print_every = 100,

)

MovingWindow(
    time_start = Main.grid_length[0]/2.-center_bunch,
    velocity_x = beta
)

LoadBalancing(
    initial_balance = False,
        every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
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


Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Ey','Ez','Rho','Jx']

DiagFields(
    every = 1,
        fields = list_fields
)

DiagProbe(
        every = 10,
        origin = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
        corners = [
            [Main.grid_length[0], Main.grid_length[1]/2., Main.grid_length[2]/2.]
        ],
        number = [nx],
        fields = ['Ex','Ey','Rho','Jx']
)

DiagProbe(
        every = 10,
        origin = [0., Main.grid_length[1]/4., Main.grid_length[2]/2.],
        corners = [
            [Main.grid_length[0], Main.grid_length[1]/4., Main.grid_length[2]/2.],
            [0., 3*Main.grid_length[1]/4., Main.grid_length[2]/2.],
        ],
        number = [nx, ntrans],
        fields = ['Ex','Ey','Rho','Jx']
)

#DiagScalar(every = 10, vars=['Uelm','Ukin_electron','ExMax','ExMaxCell','EyMax','EyMaxCell', 'RhoMin', 'RhoMinCell'])

#DiagParticleBinning(
#       deposited_quantity = "weight_charge",
#       every = 50,
#       species = ["electron"],
#       axes = [
#               ["moving_x", 0, Main.grid_length[0], nx],
#               ["px", -1, 2., 100]
#       ]
#)
                                                                                                                                                                 
