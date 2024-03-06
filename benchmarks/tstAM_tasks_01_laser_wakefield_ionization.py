# ----------------------------------------------------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

from numpy import pi, linspace, random, meshgrid, cos, sin, vstack
import scipy.constants

dx = 0.251327
dtrans = 1.96349
dt = 0.96 * dx
nx =  1024
ntrans = 256
Lx = nx * dx
Ltrans = ntrans * dtrans
npatch_x = 32
npatch_trans = 16
Nit = 1000

ne = 0.0045
begin_upramp = 10.  #Density is 0 before that and up ramp starts.
Lupramp = 100. #Length of the upramp 
Lplateau = 2000.  #Length of the plateau 
Ldownramp = 2356.19 #Length of the down ramp
xplateau = begin_upramp + Lupramp # Start of the plateau
begin_downramp = xplateau + Lplateau # Beginning of the output ramp. 
finish = begin_downramp + Ldownramp # End of plasma

Lmu  = 0.8 # laser wavelength

Main(
    geometry = "AMcylindrical",
    number_of_AM=2,
    interpolation_order = 2,
    timestep = dt,
    simulation_time = dt*Nit,
    cell_length  = [dx, dtrans],
    grid_length = [ Lx,  Ltrans],

    number_of_patches = [npatch_x, npatch_trans],

    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    
    solve_poisson = False,
    print_every = 100,
    reference_angular_frequency_SI = 2.*math.pi * 3.e8/(Lmu*1.e-6),
    cluster_width = 8,

)

MovingWindow(
    time_start = Main.grid_length[0] - 50*dx, #Leaves 2 patches untouched, in front of the laser.
    velocity_x = 0.996995486
)


Radius_plasma = 250.
longitudinal_profile = polygonal(xpoints=[begin_upramp,xplateau,begin_downramp,finish],xvalues=[0.,ne,ne,0.])
def my_profile(x,r):
        profile_r = 0.
        if ((r)**2<Radius_plasma**2):
                profile_r = 1.
        return profile_r*longitudinal_profile(x,r)


Species( 
    name = "electronfromion",
    position_initialization = "regular",
    momentum_initialization = "cold",
    ionization_model = "none",
    particles_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.,  
    mean_velocity = [0., 0., 0.],
    time_frozen = 0.0,
    pusher = "boris",
    boundary_conditions = [
       ["remove", "remove"],
       ["reflective", "remove"],
    ], 
  

)

Species( 
    name = "nitrogen5plus",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 8,
    regular_number = [1,1,8], 
    atomic_number = 7,
    ionization_model = "tunnel",
    #ionization_model = "none",
    ionization_electrons = "electronfromion",
    maximum_charge_state = 7,
    c_part_max = 1.0,
    mass = 7*2000+ 7*2000  + 2,
    charge = 5.0,
    charge_density = my_profile,  
    mean_velocity = [0., 0., 0.],
    time_frozen = 20000.0,
    pusher = "boris",
    boundary_conditions = [
       ["remove", "remove"],
       ["reflective", "remove"],
    ],
)

Species( 
    name = "neutralizingelectron",
    position_initialization = "nitrogen5plus",
    #position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 8, 
    ionization_model = "none",
    c_part_max = 1.0,
    mass = 1,
    charge = -1,
    charge_density = my_profile,  
    mean_velocity = [0., 0., 0.],
    time_frozen = 0.0,
    pusher = "boris",
    boundary_conditions = [
       ["remove", "remove"],
       ["reflective", "remove"],
    ],

)

laser_fwhm = 82. 
LaserGaussianAM(
    box_side         = "xmin",
    a0              = 3.,
    focus           = [10.],  
    waist           = 120.,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)


LoadBalancing(
    initial_balance = False,
        every = 20,
    cell_load = 1.,
    #frozen_particle_load = 0.1
)


DiagProbe(
        every = 100,
        origin = [0., 0, 0.],
        corners = [
              [Lx, 0, 0.],
                  ],
        number = [nx],
        fields = ["Ex","Ey","Rho_electronfromion","Rho_neutralizingelectron","Rho"],
        
)


DiagProbe(
	every = 100,
	origin = [0., -Ltrans, 0.],
	corners = [
              [Lx, -Ltrans, 0.],
              [0.,  Ltrans, 0.]
                  ],
	number = [nx, 2*ntrans],
        fields = ["Ex","Ey","Rho_electronfromion","Rho_neutralizingelectron","Rho"],
)

