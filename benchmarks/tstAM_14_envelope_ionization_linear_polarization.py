import math
import scipy.constants



################### AM ionization with Laser Envelope
dx = 0.8
dtrans = 2.
dt = 0.6386  
nx = 2560
ntrans = 192
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 256
laser_fwhm = 92.5


Lmu  = 0.8 # laser wavelength

mp_over_me = scipy.constants.proton_mass / scipy.constants.electron_mass
mn_over_me = scipy.constants.neutron_mass / scipy.constants.electron_mass


Main(
    geometry = "AMcylindrical",

    interpolation_order = 2,
    number_of_AM=1,
    timestep = dt,
    simulation_time = 2002*dt, #6100, #350.*dt,

    cell_length  = [dx, dtrans],
    grid_length = [ Lx,  Ltrans],

    number_of_patches =[npatch_x, 32],
    
    clrw = nx/npatch_x,

    EM_boundary_conditions = [
       ["silver-muller","silver-muller"],
       ["buneman","buneman"],
    ],

    solve_poisson = False,
    print_every = 100,
    reference_angular_frequency_SI = 2.*math.pi * 3.e8/(Lmu*1.e-6),

    random_seed = smilei_mpi_rank
)


LoadBalancing(
    initial_balance = False,
        every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)


ne = 0.0000056

begin_upramp = 2.5*350.+2.*laser_fwhm  
Lupramp = 0. #Length of the upramp 
Lplateau = 120.  #Length of the plateau 
Ldownramp = 0.  #Length of the down ramp
xplateau = begin_upramp + Lupramp # Start of the plateau
begin_downramp = xplateau + Lplateau # Beginning of the output ramp. 
finish = begin_downramp + Ldownramp # End of plasma


g = polygonal(xpoints=[begin_upramp, xplateau, begin_downramp, finish], xvalues=[0, ne, ne, 0.])
R_plasma = 100. 


def my_profile(x,y):
        radial_profile = 1.
        if (y>R_plasma):
		radial_profile = 0.
	return g(x,y)*radial_profile


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
    ponderomotive_dynamics = "True",
    boundary_conditions = [
       ["remove", "remove"],
       ["reflective", "remove"],
    ], 
  

)

Species( 
    name = "nitrogen5plus",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 15, 
    atomic_number = 7,
    ionization_model = "tunnel_envelope_averaged",
    #ionization_model = "none",
    ionization_electrons = "electronfromion",
    maximum_charge_state = 7,
    c_part_max = 1.0,
    mass = 7*mp_over_me + 7*mn_over_me + 2,
    charge = 5.0,
    charge_density = my_profile,  
    mean_velocity = [0., 0., 0.],
    time_frozen = 20000.0,
    pusher = "ponderomotive_boris",
    ponderomotive_dynamics = "True",
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
    particles_per_cell = 15, 
    ionization_model = "none",
    c_part_max = 1.0,
    mass = 1,
    charge = -1,
    charge_density = my_profile,  
    mean_velocity = [0., 0., 0.],
    time_frozen = 0.0,
    pusher = "ponderomotive_boris",
    ponderomotive_dynamics = "True",
    boundary_conditions = [
       ["remove", "remove"],
       ["reflective", "remove"],
    ],

)

LaserEnvelopeGaussianAM( 
    a0              = 2.,     
    focus           = [begin_upramp, 0.],
    waist           = 100.,
    time_envelope   = tgaussian(center=(2.0*laser_fwhm), fwhm=laser_fwhm),
    Envelope_boundary_conditions = [ ["reflective"] ],
    polarization_phi = 0.,
    ellipticity      = 0.
  
)

Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)


field_lists_forprobes=["Ex","Ey","Rho","Jx","Jy","Jz","Env_E_abs","Rho_neutralizingelectron","Rho_electronfromion","Rho_nitrogen5plus"]

DiagProbe(	
	every = 100,
	
	origin = [0., 2*dtrans, 0.],
	
	corners = [
              [Main.grid_length[0], 2*dtrans, 0.]
                  ],

	number = [nx],
	fields = field_lists_forprobes,
)
	
DiagProbe(
        every = 100,

        origin = [0., 0., 0.],

        corners = [
              [Main.grid_length[0], 0., 0.]
                  ],

        number = [nx],
        fields = field_lists_forprobes,
)
	
DiagProbe(	
	every = 100,
	
	origin = [0., -Main.grid_length[1], 0.],
	
	corners =  [	
           [Main.grid_length[0], -Main.grid_length[1], 0.],	
	   [0., Main.grid_length[1], 0.],
                   ],
	
 	number = [nx,2*ntrans],	
        fields = field_lists_forprobes,
)

DiagTrackParticles(
    species = "electronfromion",
    every = 2000,
    attributes = ["x","y", "z", "px", "py", "pz","weight"]
)


DiagParticleBinning(
    deposited_quantity = "weight",
    every = 2000,
    time_average = 1,
    species = ["electronfromion"],
    axes = [ ["px",     -0.3,  2.0,    100], ],
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 2000,
    time_average = 1,
    species = ["electronfromion"],
    axes = [ ["py",    -2.1,  2.1,    50], ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 2000,
    time_average = 1,
    species = ["electronfromion"],
    axes = [ ["pz",    -0.7,  0.7,    50] ]
)

