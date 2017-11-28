# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

# -------------------
# MY PYTHON VARIABLES
# -------------------

import math as m 

v0 = 0.9   # electron flow velocity
k  = 2.    # seeded wave number 

#T0  = 1.e-4                 # electron temperature
#Lde = m.sqrt(T0)            # Debye length [non-relativistic for ne=1]
g0  = 1./m.sqrt(1.-v0**2)   # Lorentz factor corresponding to v0

Lx  = 8*2*m.pi/k 	        # simulation length (8 seeded wavelength)

dx   = 2.*m.pi/k / 64.; resx = 1./dx   # spatial resolution in the x-direction = 500 punti per lunghezza d onda
dt   = 0.95*dx                         # time-step
t_sim     = 2000*dt                    # simulation time
nt        = t_sim/dt                   # nb of timesteps

nppc    = 200

# define the function used to seed the Weibel unstable mode
By0     = 0.001
def By(x):
	return By0*m.sin(k*x)

# --------------------------------------
# SMILEI's VARIABLES (defined in blocks)
# --------------------------------------

# MAIN BLOCK
Main(
    geometry = "1Dcartesian",

    interpolation_order = 2,

    timestep = dt,
    simulation_time = t_sim,

    cell_length = [dx],
    grid_length  = [Lx],

    number_of_patches = [ 8 ],

    EM_boundary_conditions = [ ['periodic'] ] ,

    random_seed = smilei_mpi_rank
)

# DEFINE ALL SPECIES (ADDING BLOCKS)

Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = 1836.,
    charge = 1.0,
    number_density = 1.,
    boundary_conditions = [
        ['periodic'],
    ],
    time_frozen = 2.*t_sim
)

Species(
    name = 'eon1',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.0,
    number_density = 0.5,
    mean_velocity = [0.,0.,v0],
    boundary_conditions = [
        ['periodic'],
    ]
)

Species(
    name = 'eon2',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.0,
    number_density = 0.5,
    mean_velocity = [0.,0.,-v0],
    boundary_conditions = [
        ['periodic'],
    ]
)

# APPLY EXTERNAL FIELD (to seeded unstable mode)

ExternalField(
	field='By',
    profile = By
)

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

DiagScalar(every=1)

DiagFields(
    every = 4,
    fields = ['Ex','Ey','Ez','By_m','Rho','Rho_eon1','Rho_eon2','Jz','Jz_eon1','Jz_eon2']
)


