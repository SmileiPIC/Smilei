############################# Input namelist for laser guiding in a plasma channel
import math 
import numpy as np
import scipy.constants

##### Physical constants
lambda0                            = 0.8e-6                    # laser wavelength, m
c                                  = scipy.constants.c         # lightspeed, m/s
omega0                             = 2*math.pi*c/lambda0       # laser angular frequency, rad/s
eps0                               = scipy.constants.epsilon_0 # Vacuum permittivity, F/m
e                                  = scipy.constants.e         # Elementary charge, C
me                                 = scipy.constants.m_e       # Electron mass, kg
ncrit                              = eps0*omega0**2*me/e**2    # Plasma critical number density, m-3
c_over_omega0                      = lambda0/2./math.pi        # converts from c/omega0 units to m
reference_frequency                = omega0                    # reference frequency, s-1
E0                                 = me*omega0*c/e             # reference electric field, V/m

##### Variables used for unit conversions
c_normalized                       = 1.                        # speed of light in vacuum in normalized units
um                                 = 1.e-6/c_over_omega0       # 1 micron in normalized units
mm                                 = 1.e-3/c_over_omega0       # 1 mm in normalized units
fs                                 = 1.e-15*omega0             # 1 femtosecond in normalized units

#########################  Simulation parameters

# This mesh is particularly coarse and can cause artefacts near the axis.
# A more accurate simulation would use a finer mesh

##### mesh resolution
dr                                 = 1*um                      # transverse mesh resolution
dx                                 = 0.03*um                   # longitudinal mesh resolution

##### simulation window size
nr                                 = 80                        # number of mesh points in the transverse direction
Lr                                 = nr * dr                   # transverse size of the simulation window

nx                                 = 1792                      # number of mesh points in the longitudinal direction
Lx                                 = nx * dx                   # longitudinal size of the simulation window  

##### Total simulation time
dt                                 = 0.99*dx/c_normalized      # integration timestep
T_sim                              = 3*Lx                        # total simulation time                       

##### patches parameters (parallelization)
npatch_r                           = 4
npatch_x                           = 128                       # patches in the x direction (parallelization)



######################### Main simulation definition block

Main(
    geometry                       = "AMcylindrical",

    interpolation_order            = 1,

    timestep                       = dt,
    simulation_time                = T_sim,

    cell_length                    = [dx, dr],
    grid_length                    = [ Lx,  Lr],

    number_of_AM                   = 2,

    number_of_patches              = [npatch_x,npatch_r],
 
    EM_boundary_conditions         = [["silver-muller"],["buneman"],],
    
 
    solve_poisson                  = False,
    print_every                    = 100,
    use_BTIS3_interpolation        = True,

    reference_angular_frequency_SI = omega0,
    # to see the difference with the Yee solver, increase a lot the simulation time
    # and change the maxwell_solver to "Yee"
    maxwell_solver                 = "Terzani", 
    

    random_seed                    = smilei_mpi_rank
)

######################### Define the laser pulse

### laser parameters
laser_fwhm                         = 28*math.sqrt(2)*fs            # laser fwhm duration in field
laser_w_0                          = 14*um                         # laser waist
center_laser                       = 2*laser_fwhm*c_normalized     # time when the peak of the laser is injected from the border
a0                                 = 1.                            # peak normalized transverse field of the laser
focus                              = [2*um]                        # laser focus, [x] position
omega                              = omega0/omega0                 # normalized laser frequency (=1 since its frequency omega0 is also the normalizing frequency)


LaserGaussianAM(
    box_side                       = "xmin",
    a0                             = a0,
    omega                          = omega,
    focus                          = focus,
    waist                          = laser_w_0,
    ellipticity                    = 0., # linear polarization
    polarization_phi               = 0., # in the y direction
    time_envelope                  = tgaussian(fwhm=laser_fwhm, center=center_laser),
)


######################### Define a moving window

MovingWindow(
    time_start                     = Lx, 
    velocity_x                     = c_normalized,
)

########################## Define the plasma channel for guiding

###### plasma parameters
Radius_plasma                      = 60.*um                        # Radius of plasma

Lramp                              = 2.*um                         # Plasma density upramp length
Lplateau                           = 300*mm                        # Length of density plateau
Ldownramp                          = 0.*um                         # Length of density downramp

x_begin_upramp                     = 0.
x_begin_plateau                    = x_begin_upramp+Lramp          # x coordinate of the end of the density upramp / start of density plateau
x_end_plateau                      = x_begin_plateau+Lplateau      # x coordinate of the end of the density plateau start of the density downramp
x_end_downramp                     = x_end_plateau+Ldownramp       # x coordinate of the end of the density downramp

##### plasma density profile, uniform along x, parabolic along r
longitudinal_profile               = polygonal(xpoints=[x_begin_upramp,x_begin_plateau,x_end_plateau,x_end_downramp],xvalues=[0.,1.,1.,0.])
lambda_p                           = 20.e-6                        # plasma wavelength, m
omega_p                            = 2*math.pi*c/lambda_p          # plasma angular frequency, rad/s
n0_center_SI                       = eps0*omega_p**2*me/e**2       # Plasma number density at the center, m-3
n0                                 = n0_center_SI/ncrit            # plateau number density at the center, code units
kp                                 = np.sqrt(e**2 * n0_center_SI/me/eps0 ) / c  * c_over_omega0
R_channel                          = laser_w_0                     # matched channel
kp_sq_times_R_fourth_power         = ((kp**2))*((R_channel)**4)

def plasma_density(x,r):
	profile_r                      = 0.
	if (r<Radius_plasma):
		profile_r                  = 1.+ 4.* r**2 / kp_sq_times_R_fourth_power
	return n0*profile_r*longitudinal_profile(x,r)

###### define the plasma electrons
Species(
 name                              = "plasmaelectrons",
 position_initialization           = "regular",
 momentum_initialization           = "cold",
 particles_per_cell                = 8,
 regular_number                    = [1,1,8],
 mass                              = 1.0,
 charge                            = -1.0,
 number_density                    = plasma_density,
 mean_velocity                     = [0.0, 0.0, 0.0],
 temperature                       = [0.,0.,0.],
 pusher                            = "borisBTIS3",
 time_frozen                       = 0.0,
 boundary_conditions               = [["remove", "remove"],["remove", "remove"],],
)


######################### Diagnostics

output_every = int(50*um/dt)

list_fields_diagnostic             = ['Ex','Ey','Rho','Bz','BzBTIS3','Rho_plasmaelectrons']

##### 1D Probe diagnostic on the x axis
DiagProbe(
        every                      = output_every,
        origin                     = [ 0.                  , 1.*dr, 1.*dr],
        corners                    = [ [Main.grid_length[0], 1.*dr, 1.*dr]],
        number                     = [ nx ],
        fields                     = list_fields_diagnostic
)

##### 2D Probe diagnostics on the xy plane
DiagProbe(
    every                          = output_every,
    origin                         = [ 0.                  , -nr*dr,0.],
    corners                        = [ [nx*dx,-nr*dr,0.]   , [0,nr*dr,0.] ],
    number                         = [ nx                  , int(2*nr)],
    fields                         = list_fields_diagnostic
)
