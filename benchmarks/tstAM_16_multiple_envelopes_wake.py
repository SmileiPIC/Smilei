# Multiple laser envelopes propagation in plasma

import cmath
from numpy import exp, sqrt, arctan, vectorize

dx = 0.69 
dr = 5. 
dt = 0.57#0.8*dx
nx = 1024
nr = 60
Lx = nx * dx
Lr = nr*dr
npatch_x=64

time_start_moving_window = 0.


Main(
    geometry = "AMcylindrical",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = 1700.*dt,

    cell_length  = [dx, dr],
    grid_length = [ Lx,  Lr],

    number_of_AM = 1,

    number_of_patches = [npatch_x,4],
    clrw = nx/npatch_x,

    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],

    solve_poisson = False,
    print_every = 100,

    random_seed = smilei_mpi_rank
)

MovingWindow(
    time_start = time_start_moving_window,
    velocity_x = 1.0
)

LoadBalancing(
    initial_balance = False,
        every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)


####### Define multiple laser envelopes
####### IMPORTANT: they represent laser pulses with the same carrier frequency omega, BC, ellipticity and polarization angle
####### Their initial vector potentials can thus be summed at t=0 and be evolved with the same envelope solver

# parameters in common for the envelopes
omega = 1.
ellipticity = 0.
polarization_phi = 0.
envelope_solver = "explicit"
Envelope_boundary_conditions = [ ["reflective", "reflective"],["reflective", "reflective"], ]

# a0 of the envelopes
a0_0 = 0.01
a0_1 = 0.02
a0_2 = 0.03

a0_envelopes = [a0_0,a0_1,a0_2] 

Number_of_envelopes = len(a0_envelopes)

# fwhm duration of the envelopes

laser_fwhm_0 = 10. #70.
laser_fwhm_1 = 10. #75.
laser_fwhm_2 = 10. #80.

laser_fwhm_envelopes = [laser_fwhm_0,laser_fwhm_1,laser_fwhm_2]

# waist of the envelopes
waist_0 = 70.
waist_1 = 80.
waist_2 = 90.

waist_envelopes = [waist_0,waist_1,waist_2]

# initial position of the envelopes

center_laser_0 = Lx-150.
center_laser_1 = center_laser_0-150.
center_laser_2 = center_laser_1-150.

center_laser_envelopes = [center_laser_0,center_laser_1,center_laser_2]

# focus of the envelopes
focus_0           = [center_laser_0, 0.]
focus_1           = [center_laser_1, 0.]
focus_2           = [center_laser_2, 0.]

focus_envelopes = [focus_0,focus_1,focus_2]

# temporal profile of the envelopes
time_envelope_0 = tgaussian(center=center_laser_envelopes[0], fwhm=laser_fwhm_envelopes[0])
time_envelope_1 = tgaussian(center=center_laser_envelopes[1], fwhm=laser_fwhm_envelopes[1])
time_envelope_2 = tgaussian(center=center_laser_envelopes[2], fwhm=laser_fwhm_envelopes[2])

time_envelope_envelopes = [time_envelope_0,time_envelope_1,time_envelope_2]

# profile of the vector potential of the single envelope

def gaussian_beam_with_temporal_profile(x,r,t,i_envelope):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        Zr = omega * waist_envelopes[i_envelope]**2/2.
        w  = sqrt(1./(1.+   ( (x-focus_envelopes[i_envelope][0])/Zr  )**2 ) )
        coeff = omega * (x-focus_envelopes[i_envelope][0]) * w**2 / (2.*Zr**2)
        phase = coeff * ( r**2 )
        exponential_with_total_phase = exp(1j*(phase-arctan( (x-focus_envelopes[i_envelope][0])/Zr )))
        invWaist2 = (w/waist_envelopes[i_envelope])**2
        spatial_amplitude = a0_envelopes[i_envelope]*omega *polarization_amplitude_factor* w * exp( -invWaist2*(  r**2  ) )
        space_time_envelope = spatial_amplitude * vectorize(time_envelope_envelopes[i_envelope])(t)
        return space_time_envelope * exponential_with_total_phase

# sum the vector potentials of the single envelopes

def summed_envelope_profiles(x,r,t):
	total_vector_potential = 0.
	for i_envelope in range(0,Number_of_envelopes):
		total_vector_potential += gaussian_beam_with_temporal_profile(x,r,t,i_envelope)
	return total_vector_potential 


LaserEnvelope(
        omega               = omega,
        envelope_profile    = summed_envelope_profiles,
        envelope_solver     = envelope_solver,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi    = polarization_phi,
        ellipticity         = ellipticity
    )


########### Define the plasma

n0 = 0.0017
Radius_plasma = 300.
longitudinal_profile = polygonal(xpoints=[center_laser_0+laser_fwhm_0,center_laser_0+1.1*laser_fwhm_0,240000,240000],xvalues=[0.,n0,n0,0.])
def nplasma(x,r):
    profile_r = 0.
    if ((r)**2<Radius_plasma**2):
        profile_r = 1.
    return profile_r*longitudinal_profile(x,r)



Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 4,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = nplasma,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.,0.,0.],
    pusher = "ponderomotive_boris",
    ponderomotive_dynamics = "True",
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
       ["remove", "remove"],
    ],
)


Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)


DiagFields(
    every = 100,
    fields = ['Env_A_abs_mode_0','Env_Chi_mode_0','Env_E_abs_mode_0' ]
)

DiagFields(
    every = 100,
    fields = ['El_mode_0' ]
)

DiagProbe(
        every = 50,
        origin = [0., 2.*dr, 2.*dr],
        corners = [
            [Main.grid_length[0], 2.*dr, 2.*dr]
        ],
        number = [nx],
        fields = ['Ex','Ey','Rho','Jx','Rho_electron','Jx_electron','Env_A_abs','Env_Chi','Env_E_abs']
)


DiagProbe(
    every = 50,
    origin   = [0., -nr*dr,0.],
    corners  = [ [nx*dx,-nr*dr,0.], [0,nr*dr,0.] ],
    number   = [nx, 2*nr],
    fields = ['Ex','Ey','Rho','Jx','Rho_electron','Jx_electron','Env_A_abs','Env_Chi','Env_E_abs']
)



#DiagScalar(every = 10, vars=['Env_A_absMax','Env_E_absMax'])


                                                                                                                                                                 


