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
l0 = 2.0*math.pi	# wavelength in normalized units
t0 = l0				# optical cycle in normalized units
rest = 200.0		# nb of timestep in 1 optical cycle
resx = 100.0		# nb cells in 1 wavelength
Lsim = l0		# simulation length
Tsim = 5.0*t0		# duration of the simulation


# wavelength_SI: used by Fred Diags. (MG: should be removed at some point)
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
 
# SIMULATION BOX : for all space directions (use vector)
# cell_length: length of the cell
# sim_length: length of the simulation in units of the normalization wavelength 
#
cell_length = [l0/resx]
sim_length  = [Lsim]

# SIMULATION TIME
# timestep: duration of the timestep
# sim_time: duration of the simulation in units of the normalization period 
#
timestep = t0/rest
sim_time = Tsim
 
# ELECTROMAGNETIC BOUNDARY CONDITIONS
# bc_em_type_x/y/z : boundary conditions used for EM fields 
#                    periodic = periodic BC (using MPI topology)
#                    silver-muller = injecting/absorbing BC
#                    reflective = consider the ghost-cells as a perfect conductor
#
bc_em_type_x = ['silver-muller']
 
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
fieldDump_every=1

Species(
	species_type = 'helium',
	ionization_model = 'tunnel',
	atomic_number = 2,
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 100,
	mass = 1836.0,
	charge = 0.0,
	nb_density = trapezoidal(1.0,xvacuum=0.49*l0,xplateau=0.02*l0),
	bc_part_type_west = 'none',
	bc_part_type_east = 'none'
)
Species(
	species_type = 'electron',
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	bc_part_type_west = 'none',
	bc_part_type_east = 'none'
)

# LASER PROPERTIES
# for each laser define:
# a0: maximum amplitude of the laser electric field (in units of the normalization field)
# angle: angle (in degree) at which the laser enters the simulation box
# delta: polarization parameter, (0:y) (1:z) (0.707106781:circ)
# time_profile: string defining the time profile
# double_params: vector of real parameters used by the different time-profiles
# 
Laser(
    a0 = 0.1,
    boxSide = 'west',
    delta = 1.0,
    time_profile = 'constant'
)
 
# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------
#
# SCALAR DIAGNOSTICS
# every       = integer, nb of timesteps between each output
# tmin & tmax = floats, min & max times between which scalars are computed (optional)
# precision   = integer, nb of digits for the outputs (default=10)
DiagScalar(every = 10)
#print_every=10000000
 
# --------- 
# DUMP INFO
# ---------
#
# dump_step = 2500
# dump_minutes = 1
# random_seed = 13121977


print_every=1

# DiagParticles(
# 	output = "density",
# 	every = 50,
# 	species = ["electron"],
# 	axes = [
# 		["x",  0.45*Lsim, 0.55*Lsim, 200],
# 		["px", -0.1, 0.1, 200]
# 	]
# )
# 
