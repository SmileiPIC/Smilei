# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

from numpy import pi, linspace, random, meshgrid, cos, sin, vstack
import scipy.constants

dx = 0.8
dtrans = 1.96349
dt = 0.8 * dx
nx =  512
ntrans = 128
Lx = nx * dx
Ltrans = ntrans * dtrans
npatch_x = 64
npatch_trans = 8
Nit = 1100

ne = 0.0017
begin_upramp = Lx  #Density is 0 before that and up ramp starts.
Lupramp = 10. #Length of the upramp 
Lplateau = 2000.  #Length of the plateau 
Ldownramp = 2356.19 #Length of the down ramp
xplateau = begin_upramp + Lupramp # Start of the plateau
begin_downramp = xplateau + Lplateau # Beginning of the output ramp. 
finish = begin_downramp + Ldownramp # End of plasma

mp_over_me = scipy.constants.proton_mass / scipy.constants.electron_mass
mn_over_me = scipy.constants.neutron_mass / scipy.constants.electron_mass

Lmu  = 0.8

Main(
    geometry = "AMcylindrical",
    number_of_AM=1,
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
    reference_angular_frequency_SI = 2.*math.pi * scipy.constants.c /(Lmu*1.e-6),
    cluster_width = 2,

)

MovingWindow(
    time_start = 0., #Leaves 2 patches untouched, in front of the laser.
    velocity_x = 0.996995486
)


Radius_plasma = 200.
longitudinal_profile = polygonal(xpoints=[begin_upramp,xplateau,begin_downramp,finish],xvalues=[0.,ne,ne,0.])
def my_profile(x,r):
        profile_r = 0.
        if ((r)**2<Radius_plasma**2):
                profile_r = 1.
        return profile_r*longitudinal_profile(x,r)

def my_profile_dopant(x,r):
        profile_r = 0.
        if ((r)**2<Radius_plasma**2):
                profile_r = 1.
        return profile_r*longitudinal_profile(x,r)*0.1


Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    ionization_model = "none",
    particles_per_cell = 16,
    regular_number = [1,2,8],
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = my_profile,
    mean_velocity = [0., 0., 0.],
    time_frozen = 0.0,
    pusher="ponderomotive_boris",
    boundary_conditions = [
    	["remove", "remove"],
    	["reflective", "remove"],
    ],
)

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
    pusher = "ponderomotive_boris",
    boundary_conditions = [
       ["remove", "remove"],
       ["reflective", "remove"],
    ], 
  

)

Species( 
    name = "nitrogen5plus",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 16, 
    regular_number = [1,2,8],
    atomic_number = 7,
    ionization_model = "tunnel_envelope_averaged",
    #ionization_model = "none",
    ionization_electrons = "electronfromion",
    maximum_charge_state = 7,
    c_part_max = 1.0,
    mass = 7*mp_over_me + 7*mn_over_me + 2,
    charge = 5.0,
    charge_density = my_profile_dopant,  
    mean_velocity = [0., 0., 0.],
    time_frozen = 20000.0,
    pusher = "ponderomotive_boris",
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
    particles_per_cell = 16, 
    ionization_model = "none",
    c_part_max = 1.0,
    mass = 1,
    charge = -1,
    charge_density = my_profile_dopant,  
    mean_velocity = [0., 0., 0.],
    time_frozen = 0.0,
    pusher = "ponderomotive_boris",
    boundary_conditions = [
       ["remove", "remove"],
       ["reflective", "remove"],
    ],

)

laser_fwhm = 120. 
LaserEnvelopeGaussianAM(
    a0              = 2.,
    focus           = [10.],  
    waist           = 95.,
    time_envelope   = tgaussian(center=(Lx-1.8*laser_fwhm), fwhm=laser_fwhm)
)


LoadBalancing(
    initial_balance = False,
        every = 30,
    cell_load = 1.,
    #frozen_particle_load = 0.1
)


DiagProbe(
        every = 500,
        origin = [0., 0, 0.],
        corners = [
              [Lx, 0, 0.],
                  ],
        number = [nx],
        fields = ["Ex","Ey","Rho","Env_E_abs","Rho_electronfromion"],
)


DiagProbe(
	every = 500,
	origin = [0., -Ltrans, 0.],
	corners = [
              [Lx, -Ltrans, 0.],
              [0.,  Ltrans, 0.]
                  ],
	number = [nx, 2*ntrans],
	fields = ["Ex","Ey","Rho","Env_E_abs","Rho_electronfromion"],
)



