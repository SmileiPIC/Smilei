# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# CAUTION:  never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField
#

# MY PYTHON VARIABLES
# here are defined some useful python variables
# 
import math

Te_keV = 1.              # electron temperature in keV
Te  = Te_keV/511.        # Te normalised in mec^2 (code units)
vth = math.sqrt(Te)      # normalised thermal velocity
Ld    = vth              # Debye length in normalised units
dx  = Ld/10.             # spatial resolution
Lsim = 100.*Ld           # simulation length
tsim = 50.               # duration of the simulation

mi = 100.0               # ion mass (use reduced one to accelerate computation)
cs = math.sqrt(Te/mi)    # ion acoustic velocity (normalised to c)
Uion = 5.*mi*cs*cs       # max. energy used to compute the ion spectrum


# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
#
geometry = "1d3v"

# order of interpolation
#
interpolation_order = 2

# SIMULATION TIME 
# either use the resolution (res_time) or time-step (timestep)
#
# res_time  : temporal resolution 
# time step : time step
# sim_time  : duration of the simulation 
#
timestep = 0.95*dx
sim_time = tsim

# SIMULATION BOX : for all space directions (in 2D & 3D use vector of doubles)
# either use the resolution (res_space) or cell-length (cell_length)
#
# res_space   : spatial resolution 
# sim_length  : length of the simulation 
# cell_length : cell-length 
#
cell_length = [dx]
sim_length  = [Lsim]

number_of_patches = [ 1 ] # or 8

# ELECTROMAGNETIC BOUNDARY CONDITIONS
# bc_em_type_long/trans : boundary conditions used for EM fields 
#                         in the longitudinal or transverse directions
#                         periodic      = periodic BC (using MPI topology)
#                         silver-muller = injecting/absorbing
#
bc_em_type_x = ['silver-muller','silver-muller'] 

# RANDOM seed 
# this is used to randomize the random number generator
#
random_seed = 0

# DEFINE ALL SPECIES
# species_type       = string, given name to the species (e.g. ion, electron, positron, test ...)
# initPosition_type  = string, "regular" or "random"
# initMomentum_type  = string "cold", "maxwell-juettner" or "rectangular"
# c_part_max         = float, factor on the memory reserved for the total number of particles
# mass               = float, particle mass in units of the electron mass
# dynamics_type      = string, type of species dynamics = "norm" or "rrLL"
# time_frozen        = float, time during which particles are frozen in units of the normalization time
# radiating          = boolean, if true, incoherent radiation calculated using the Larmor formula 
# n_part_per_cell    = integer or function, number of particles/cell
# charge             = float or function, particle charge in units of the electron charge
# charge_density     = float or function, species charge density in units of the "critical" density
#     or nb_density for number density
# mean_velocity      = list of floats or functions, mean velocity in units of the speed of light
# temperature        = list of floats or functions, temperature in units of m_e c^2
# Predefined functions: constant, trapezoidal, gaussian, polygonal, cosine
#
Species(
    species_type = 'ion',
    initPosition_type = 'regular',
    initMomentum_type = 'mj',
    n_part_per_cell = 10,
    mass = mi, 
    charge = 1.0,
    nb_density = trapezoidal(1., xplateau=20.*Ld),
    temperature = [1.e-6],
	thermT = [1.e-6],
	thermVelocity = [0.,0.,0.],
    bc_part_type_west = 'thermalize',
    bc_part_type_east = 'refl'
)
Species(
    species_type = 'eon',
    initPosition_type = 'regular',
    initMomentum_type = 'maxwell-juettner',
    n_part_per_cell = 10,
    mass = 1.0,
    charge = -1.0,
    nb_density = trapezoidal(1., xplateau=20.*Ld),
    temperature = [Te],
    thermT = [Te],
    thermVelocity = [0.,0.,0.],
    bc_part_type_west = 'thermalize',
    bc_part_type_east = 'refl'
)


# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

every=100

# SCALAR DIAGNOSTICS
# every       = integer, nb of timesteps between each output
# tmin & tmax = floats, min & max times between which scalars are computed (optional)
# precision   = integer, nb of digits for the outputs (default=10)
DiagScalar(every = every)#, vars=['Utot','Ubal_norm','Uelm','Ukin','Ukin_ion','Ukin_eon'])    


# FIELD DUMPS
DiagFields(
    every = every,
    fields = ['Ex','Rho_ion','Rho_eon']
)

# PHASE-SPACE DIAGNOSTICS (new version from DiagParticles)
DiagParticles(
    output = "density",
    every = every,
    species = ["ion"],
    axes = [
        ["x", 0., Lsim, 50],
        ["px", -2.*cs, 3.*cs, 100]
    ]
)

DiagParticles(
    output = "density",
    every = every,
    time_average = 1,
    species = ["ion"],
    axes = [
        ["ekin", 0.0001, Uion, 100, "logscale", "edge_inclusive"]
    ]
)
