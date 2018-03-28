##### SMILEI INPUT FILE FOR SURFACE WAVE EXCITATION ####################################################################


###### MY OWN PYTHON VARIABLES #########################################################################################

import math

# Nb of particles per cells
nppc   = 8
ScaleY = 1

resx = 256
Nb_patches = [1,16]

Nb_of_diag_output_per_optcycle = 0.5

# normalizations
l0     = 2.*math.pi                    # normalized laser wavelength
t0     = 2.*math.pi                    # normalized laser optical cycle
one_eV = 1./511.e3                    # 1eV in normalized units (me c^2)

# -----------------------------------
# laser & grating resonant wavelength
# -----------------------------------

aL = 0.1                                       # normalized laser vector potential

theta_inc = 30. * math.pi/180.                 # laser angle of incidence [in rad.]
lg_ov_l0  = 2.                                 # grating wavelength in units of laser wavelength

Tramp     = 2.  *t0                            # ramp-up of laser temporal profile
Tsim      = 14. *t0                            # simulation time

# -----------------
# Plasma properties
# -----------------

n0 = 100.
T0 = 500. * one_eV

x0 = 4.*l0
Lp = 4.*l0
lg = lg_ov_l0 * l0
dg = 0.1*l0
kg = 2.*math.pi/lg

def n_(x,y):
    if ( x > x0 + dg*math.sin(kg*y) ):
        return n0
    else:
        return 0.

# -----------------------------------------------
# Prepares the simulation box, time & resolutions
# -----------------------------------------------

dx   = l0/resx
dy   = dx
dt   = 0.95 * dx/math.sqrt(2.)

Ly   = 4 * l0 * ScaleY
Lx   = int( (x0+Lp)/dx ) * dx


###### SMILEI VARIABLES ################################################################################################

# Main()
Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,

    simulation_time   = Tsim,
    timestep          = dt,
    grid_length       = [Lx,Ly],
    cell_length       = [dx,dy],
    number_of_patches = Nb_patches,

    EM_boundary_conditions       = [ ["silver-muller"],["periodic"] ],
    EM_boundary_conditions_k = [[math.cos(theta_inc),math.sin(theta_inc)],[-1.,0.],[0.,1.],[0.,-1.]],

    print_every = int(Tsim/dt/100.),
    random_seed = smilei_mpi_rank
)

# Laser (oblique plane-wave)
LaserGaussian2D(
    a0              = aL,
    omega           = 1.,
    focus           = [0,2*Ly],
    waist           = 1.e10,
    incidence_angle = theta_inc,
    time_envelope   = tsin2plateau(slope1=Tramp,plateau=10.*Tsim)
)

# Ion species (immobile)
Species(
    name      = "ion",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = nppc,
    mass = 1.,
    charge = 1.,
    number_density = n_,
    boundary_conditions = [["reflective","thermalize"],["periodic"]],
    thermal_boundary_temperature = [T0],
    thermal_boundary_velocity = [0.,0.,0.],
    time_frozen = 2.*Tsim
)

# Electron species
Species(
    name      = "eon",
    position_initialization = "ion",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = n_,
    temperature = [T0],
    boundary_conditions = [["reflective","thermalize"],["periodic"]],
    thermal_boundary_temperature = [T0],
    thermal_boundary_velocity = [0.,0.,0.],
    time_frozen = x0 - l0
)

### DIAGNOSTICS
globalEvery = int(t0/dt/Nb_of_diag_output_per_optcycle)

DiagScalar(every=globalEvery)

DiagFields(every=globalEvery)


