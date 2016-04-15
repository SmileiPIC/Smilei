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
L  = 1.12			# wavelength=simulation box length
dn = 0.001			# amplitude of the perturbation


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
cell_length = [0.01]
sim_length  = [L]

number_of_patches = [ 16 ] # n_space_x = 103

# SIMULATION TIME
# timestep: duration of the timestep
# sim_time: duration of the simulation in units of the normalization period 
#
timestep = 0.0095
sim_time = 50.
 
# ELECTROMAGNETIC BOUNDARY CONDITIONS
# bc_em_type_x/y/z : boundary conditions used for EM fields 
#                    periodic = periodic BC (using MPI topology)
#                    silver-muller = injecting/absorbing BC
#                    reflective = consider the ghost-cells as a perfect conductor
#
bc_em_type_x = ['periodic']
 
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
	species_type = "ion",
	initPosition_type = "regular",
	initMomentum_type = "cold",
	n_part_per_cell = 10,
	mass = 1836.0,
	charge = 1.0,
	nb_density = 1.,
	time_frozen = 10000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)
Species(
	species_type = "eon1",
	initPosition_type = "regular",
	initMomentum_type = "cold",
	n_part_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	nb_density = cosine(0.5,xamplitude=dn,xlength=L),
	mean_velocity = [-0.1,0.0,0.0],
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)
Species(
	species_type = "eon2",
	initPosition_type = "regular",
	initMomentum_type = "cold",
	n_part_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	nb_density = cosine(0.5,xamplitude=dn,xlength=L),
	mean_velocity = [0.1,0.0,0.0],
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)


# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

#every = 1000000000000
every = 100

# SCALAR DIAGNOSTICS
# every       = integer, nb of timesteps between each output
# tmin & tmax = floats, min & max times between which scalars are computed (optional)
# precision   = integer, nb of digits for the outputs (default=10)
DiagScalar(every = every, vars=['Utot','Ubal_norm','Uelm','Ukin'])	


# FIELD DUMPS
# fieldDump_every = integer, nb of timesteps between each output
# fieldsToDump    = ('string'), name of the fields to dump
fieldDump_every = every
fieldsToDump = ('Ex','Ey','Ez','By_m','Bz_m');


# PHASE-SPACE DIAGNOSTICS (new version from DiagParticles)
DiagParticles(
	output = "density",
	every = every,
	species = ["eon1","eon2"],
	axes = [
		["x", 0., L, 50],
		["px", -0.4, 0.4, 100]
	]
)
