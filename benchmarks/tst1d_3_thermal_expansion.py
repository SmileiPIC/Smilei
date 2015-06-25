# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
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

Te_keV = 1.				# electron temperature in keV
Te  = Te_keV/511.		# Te normalised in mec^2 (code units)
vth = math.sqrt(Te)		# normalised thermal velocity
Ld	= vth				# Debye length in normalised units
dx  = Ld/10.			# spatial resolution
Lsim = 100.*Ld			# simulation length
tsim = 50.;				# duration of the simulation

mi = 100.0				# ion mass (use reduced one to accelerate computation)
cs = math.sqrt(Te/mi)	# ion acoustic velocity (normalised to c)
Uion = 5.*mi*cs*cs		# max. energy used to compute the ion spectrum

# step-like density profile
def f(x):
    if x < 20.*Ld:
        return 1.0
    else :
        return 0.0

# wavelength_SI: used by Fred Diags. (MG: should be removed at some point)
#
wavelength_SI = 1.e-6

# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
#
dim = '1d3v'

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
# species_type: ion, electron, positron, test ...
# initialization_type: regular, cold or (isotrop) Maxwell-Juettner distribution
# n_part_per_cell: number of particle-per-cell
# c_part_max: factor on the memory reserved for the total number of particles
# mass: particle mass in units of the electron mass
# charge: particle charge in units of e (-e is the electron charge)
# density: species density in units of the normalization density
# mean_velocity: mean velocity of the species (3D vector) in units of the light velocity
# temperature: temperature of the species in units of m_e c^2
# dynamics_type: species type of dynamics = norm or rrLL
# time_frozen: time during which the particles are frozen in units of the normalization time
# radiating: boolean, if true incoherent radiation are calculated using the Larmor formula 
#
Species(
	species_type = 'ion',
	dens_profile = f,
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	n_part_per_cell = 10,
	mass = mi, 
	charge = 1.0,
	nb_density = 1.,
	temperature = 1.e-6,
	#time_frozen = 2.*tsim,
	bc_part_type_west = 'thermalize',
	bc_part_type_east = 'refl'
)
Species(
	species_type = 'eon',
	dens_profile = f,
	initPosition_type = 'regular',
	initMomentum_type = 'maxwell-juettner',
	n_part_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	nb_density = 1.,
	temperature = [Te],
	bc_part_type_west = 'thermalize',
	bc_part_type_east = 'refl'
)


# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

every=10

# SCALAR DIAGNOSTICS
# every       = integer, nb of timesteps between each output
# tmin & tmax = floats, min & max times between which scalars are computed (optional)
# precision   = integer, nb of digits for the outputs (default=10)
DiagScalar(every = every)#, vars=['Utot','Ubal_norm','Uelm','Ukin','Ukin_ion','Ukin_eon'])	


# FIELD DUMPS
# fieldDump_every = integer, nb of timesteps between each output
# fieldsToDump    = ('string'), name of the fields to dump
fieldDump_every = every
fieldsToDump = ('Ex','Rho_ion','Rho_eon');

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
