# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import cmath, math
import scipy.constants
from numpy import exp, sqrt, arctan, vectorize, real
from math import log
import numpy as np

# Reference wavelength, the ionizing wavelength
lambda_0  = 2.e-6    # [m]  #

omega0 = 2*math.pi*scipy.constants.c/lambda_0 # [1/s]


dx = 0.8             # longitudinal resolution
dy = 4.
dz = dy              # transverse resolution
dt = 0.8*dx         # timestep
Lx = 256*dx
Ly = 32*dy
Lz = 32*dz

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [dx,dy,dz],
    grid_length  =[Lx,Ly,Lz],
    
    number_of_patches = [ 32, 4, 4 ],
    
    timestep = dt,
    simulation_time = 1000*dt,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
        ['silver-muller'],
    ],
    
    solve_poisson = False,
    
    print_every = 100,
    cluster_width=2,
)


Radius_plasma  = 20e-6/(scipy.constants.c/omega0)


n0 = 0.0013

begin_upramp = Lx
Lupramp  = 10*dx
Lplateau = 2e15
Ldownramp = 10*dx
xplateau = begin_upramp + Lupramp # Start of the plateau
begin_downramp = xplateau + Lplateau # Beginning of the output ramp.
finish = begin_downramp + Ldownramp + Ldownramp # End of plasma
 
longitudinal_profile = polygonal(xpoints=[begin_upramp, xplateau, begin_downramp, finish], xvalues=[0, n0, n0, 0.])


def my_density_profile(x,y,z):
	radial_factor = 0
	if (((y-Ly/2.)**2+(z-Lz/2.)**2)<Radius_plasma**2):
		radial_factor = 1.
	return radial_factor*longitudinal_profile(x,y,z)


focus_x_laser = begin_upramp
fwhm_duration_laser_pulse = 20e-15 * omega0 *2**0.5
time_laser_pulse_peak_enters = fwhm_duration_laser_pulse

a0 = 1.
waist = 8.e-6/(scipy.constants.c/omega0)

distance_laser_peak_window_border = 1.7*fwhm_duration_laser_pulse

LaserEnvelopeGaussian3D(
    omega            = 2.5,
    a0               = a0,
    focus            = [focus_x_laser, Ly/2., Lz/2.],
    waist            = waist, # from um to normalized units
    time_envelope    = tgaussian(center=(Lx-distance_laser_peak_window_border), fwhm=fwhm_duration_laser_pulse),
)

distance_laser_peak_window_border = 1.7*fwhm_duration_laser_pulse
MovingWindow(
    time_start = 0.,
    velocity_x = 0.99, #np.sqrt(1.-n0),
)


Species(
    name = 'eon',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = 8,
    mass = 1.0,
    charge = -1.0,
    number_density = my_density_profile,
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
        ["remove", "remove"],
    ],
    pusher = "ponderomotive_boris",
)

LoadBalancing(
    every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)


DiagScalar(every=25)

DiagFields(
    every = 200,
    fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_eon','Rho','Env_E_abs']
)

