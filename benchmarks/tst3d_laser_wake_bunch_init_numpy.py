#### Namelist for laser wakefield acceleration, external injection of an electron bunch made of numpy particles

import numpy as np
import math
dx = 0.2
dtrans = 1.
dt = 0.19
nx = 1280
ntrans = 80
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 128

lambda0_laser     = 0.8e-6              # laser central frequency, m
c_over_omega0     = lambda_0/2./math.pi # converts from c/omega0 to um

# Plasma plateau density
n0 = 0.0017

# Laser parameters
laser_fwhm = 19.8
omega = 1.
a0 = 2.5
focus_laser = 3.*laser_fwhm
center_laser=3.*laser_fwhm
waist = 15. #94.26


# Bunch parameters
npart = 2000000
Q_bunch = -2.e-12                                        # C (negative for electrons as usual)
sigma_x = 3.                                             # bunch longitudinal rms size
sigma_r = 4.                                             # bunch transverse rms size (cylindrical symmetry)
bunch_energy_spread = 0.01                               # sigma_gamma/gamma (not in percent)
bunch_normalized_emittance = 1.e-6/c_over_omega0         # for both planes, m-rad /c_over_omega0
center_bunch = 40.                                       # bunch at 0.5*lambdap behind the laser (ok a little more because of relativistic elongation of plasma wavelength)  
gamma_bunch = 200.                                       # relativistic lorentz factor


Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2,

    timestep = dt,
    simulation_time = int(5*Lx/dt)*dt,

    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],

    number_of_patches = [npatch_x, 4, 4],

    clrw = nx/npatch_x,
    #clrw = 1,

    solve_relativistic_poisson = True,
    
    EM_boundary_conditions = [ ["silver-muller"] ],
    
    solve_poisson = False,
    print_every = 100,

)

MovingWindow(
    time_start = center_laser+focus_laser+135.,
    velocity_x = 0.9997
)


LoadBalancing(
    initial_balance = False,
    every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

# Auxiliary physical quantities
lambda0 = 0.8e-6               # laser wavelength, m
c=299792458                    # lightspeed, m/s
omega0 = 2*math.pi*c/lambda0   # laser frequency
eps0=8.854187e-12              # Vacuum permittivity, F/m
e=1.60217646e-19               # Elementary charge, C
me=9.10938215e-31              # Electron mass, kg
ncrit = eps0*omega0**2*me/e**2 # Plasma critical number density [m-3]
normalized_species_charge = -1 # For electrons

#Fixing weight from a given density and nppc
#ne = 1e18 #Target density in number of electrons per cm-3.
#normalized_density = ne*1e6/ncrit
#npcc = 8 #equivalent of 8 particle per cell
#weight = normalized_density/npcc

#Fixing weight from a given charge
#Qpart = -5e-17 #C
#cell_volume = dx*dtrans*dtrans*(c/omega0)**3
#weight = Qpart/(cell_volume*ncrit*e*normalized_species_charge)



Q_part  = Q_bunch/npart
cell_volume = dx*dtrans*dtrans*(c/omega0)**3
weight = Q_part/(cell_volume*ncrit*e*normalized_species_charge)


#npart=10000
array_position = np.zeros((4,npart))   # positions x,y,z, weight
array_momentum = np.zeros((3,npart))   # momenta x,y,z
#y = numpy.ones(npart)
#x = numpy.random.rand(npart)

random_number_x_plane = np.random.normal(0., 1., npart)  # generate random number from gaussian distribution for x position
random_number_y_plane = np.random.normal(0., 1., npart)  # generate random number from gaussian distribution for y position
random_number_z_plane = np.random.normal(0., 1., npart)  # generate random number from gaussian distribution for z position
random_number_px_plane = np.random.normal(0., 1., npart) # generate random number from gaussian distribution for px position
random_number_py_plane = np.random.normal(0., 1., npart) # generate random number from gaussian distribution for py position
random_number_pz_plane = np.random.normal(0., 1., npart) # generate random number from gaussian distribution for pz position


array_position[0,:] = np.multiply(random_number_x_plane,sigma_x)+np.multiply(center_bunch,np.ones(npart))
array_position[1,:] = np.multiply(random_number_y_plane,sigma_r)+np.multiply(Main.grid_length[1]/2,np.ones(npart))
array_position[2,:] = np.multiply(random_number_z_plane,sigma_r)+np.multiply(Main.grid_length[2]/2,np.ones(npart))
array_momentum[0,:] = np.multiply(gamma_bunch,np.ones(npart))+np.multiply(random_number_px_plane,bunch_energy_spread*gamma_bunch)
array_momentum[1,:] = np.multiply(random_number_py_plane,bunch_normalized_emittance/sigma_r)
array_momentum[2,:] = np.multiply(random_number_pz_plane,bunch_normalized_emittance/sigma_r)



#array_position[0,:] = x*1.6
#array_position[1,:] = y
#array_position[2,:] = 0.
array_position[3,:] = np.multiply(np.ones(npart),weight)
#array_momentum[0,:] = 500.
#array_momentum[1,:] = 0.
#array_momentum[2,:] = 0.

Species( 
    name = "electron_bunch",
    position_initialization = array_position,
    momentum_initialization = array_momentum,
    #particles_per_cell = 1,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    relativistic_field_initialization = True,
    #charge_density = 0.000494,       # comment when importing numpy particles
    #mean_velocity = [0.0, 0.0, 0.0], # comment when importing numpy particles
    #temperature = [0.2],             # comment when importing numpy particles
    pusher = "boris",    
    time_frozen = center_bunch+center_laser+86.,    # bunch starts moving when laser peak is at 86.c/omega0 of distance, half of the plasma wave 
    boundary_conditions = [
    	["remove", "remove"],
    	["remove", "remove"],
    	["remove", "remove"],
    ],
)


Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = polygonal(xpoints=[focus_laser+20,focus_laser+30.,15000.,15100.],xvalues=[0.,n0,n0,0.]),
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.,0.,0.],
    pusher = "boris",
    time_frozen = 0.,
    boundary_conditions = [
       ["remove", "remove"],
       ["remove", "remove"],
       ["remove", "remove"],
    ],
)






# We build a gaussian laser from scratch instead of using LaserGaussian3D
# The goal is to test the space_time_profile attribute
time_envelope = tgaussian(center=center_laser, fwhm=laser_fwhm)
focus = [focus_laser, Main.grid_length[1]/2., Main.grid_length[2]/2.]
Zr = omega * waist**2/2.
w  = math.sqrt(1./(1.+(focus[0]/Zr)**2))
invWaist2 = (w/waist)**2
coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
def By(y,z,t):
    return 0.
def Bz(y,z,t):
    r2 = (y-focus[1])**2 + (z-focus[2])**2
    omegat = omega*t - coeff*r2
    return a0 * w * math.exp( -invWaist2*r2  ) * time_envelope( omegat/omega ) * math.sin( omegat )
Laser(
    box_side = "xmin",
    space_time_profile = [By, Bz]
)

#LaserGaussian3D(
#    box_side         = "xmin",
#    a0              = 2.,
#    focus           = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
#    waist           = 10.,
#    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
#)

Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ["Ex","Ey","Ez","Bx","By","Bz","Rho","Jx","Jy","Jz"] #all fields

DiagFields(
    every = 100,
    fields = list_fields
)

DiagProbe(
	every = 100,
	origin = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
	corners = [
	    [Main.grid_length[0], Main.grid_length[1]/2., Main.grid_length[2]/2.]
	],
	number = [nx],
	fields = list_fields
)

DiagProbe(
	every = 100,
	origin = [0., 0. , Main.grid_length[2]/2.],
	corners = [
	    [Main.grid_length[0], 0., Main.grid_length[2]/2.],
	    [0., Main.grid_length[1], Main.grid_length[2]/2.],
	],
	number = [nx, ntrans],
	fields = list_fields
)

DiagScalar(every = 10, vars=['Uelm','Ukin_electron','EyMax','EyMaxCell'])

#DiagParticleBinning(
#	deposited_quantity = "weight_charge",
#	every = 100,
#	species = ["electron_bunch"],
#	axes = [
#		["moving_x", 0, Main.grid_length[0], nx],
#		["px", -1, 2., 100]
#	]
#)


