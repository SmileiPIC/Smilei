###### Namelist for field initialization of a relativistic electron bunch

import math

dx = 0.4 
dtrans = 1.
dt = 0.33
nx = 128 
ntrans = 128 
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 8
npatch_r = 8

# Plasma density
n0 = 0.0035

# Bunch position and rms dimensions (gaussian density distribution)
bunch_sigma_x = 4.
bunch_sigma_r = 3. 
center_bunch = dx*nx/2. 

# Bunch normalized density
alpha = 0.9

# Bunch mean energy
gamma= 200. # relativistic lorentz factor
beta = math.sqrt(1.-1/gamma**2)
relative_energy_spread = 0.01



# normalized density of a bunch with gaussian density
def nbunch_(x,y,z):
    r2 = (y-ntrans*dtrans/2)**2+(z-ntrans*dtrans/2)**2
    profile_x = math.exp(-(x-center_bunch)**2/2./bunch_sigma_x**2)
    profile_r = math.exp(-(r2/2./bunch_sigma_r**2))
    profile = alpha*n0*profile_x*profile_r
    x2 = (x-center_bunch)**2
    if x2/(5.*bunch_sigma_x)**2 + r2/(5.*bunch_sigma_r)**2 < 1.:
        return profile
    else:
        return 0.


Main(
    geometry = "3Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = dt*1, 

    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],


    number_of_patches = [npatch_x,npatch_r,npatch_r],
    
    clrw = nx/npatch_x,

    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],["buneman","buneman"],
    ],

    solve_poisson = False,
    
    solve_relativistic_poisson = True,
    relativistic_poisson_max_iteration = 50000,    
    print_every = 100,

    random_seed = smilei_mpi_rank
)

#MovingWindow(
#    time_start = 0.,
#    velocity_x = 1.
#)

#LoadBalancing(
#    initial_balance = False,
#        every = 20,
#    cell_load = 1.,
#    frozen_particle_load = 0.1
#)


Species(
    name = "bunch_electrons",
    position_initialization = "regular",
    momentum_initialization = "cold",
    relativistic_field_initialization = True,
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = nbunch_,
    mean_velocity = [beta, 0.0, 0.0], 
    pusher = "boris",
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
       ["remove", "remove"],["remove", "remove"],
    ],
)



list_fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho','Jx','Jy']

DiagFields(
    every = 200,
    fields = list_fields
)



DiagProbe(
    every = 200,
    origin   = [0., 0.,ntrans*dtrans/2.],
    corners  = [ [nx*dx,0.,ntrans*dtrans/2.], [0,ntrans*dtrans,ntrans*dtrans/2.] ],
    number   = [nx, ntrans]
)

