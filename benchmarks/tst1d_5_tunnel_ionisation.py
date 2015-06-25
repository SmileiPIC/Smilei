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
Lsim = 2.0*l0		# simulation length
Tsim = 10.0*t0		# duration of the simulation

def f(x):
	if (0.99*l0 < x < 1.01*l0):
		return 1.0
	else:
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
# species_type: ion, electron, positron, test ...
# initialization_type: regular, cold or (isotrop) Maxwell−Juettner distribution
# n_part_per_cell: number of particle−per−cell
# c_part_max: factor on the memory reserved for the total number of particles
# mass: particle mass in units of the electron mass
# charge: particle charge in units of e (−e is the electron charge)
# density: species density in units of the normalization density
# mean_velocity: mean velocity of the species (3D vector) in units of the light velocity
# temperature: temperature of the species in units of m_e c^2
# dynamics_type: species type of dynamics = norm or rrLL
# time_frozen: time during which the particles are frozen in units of the normalization time
# radiating: boolean, if true incoherent radiation are calculated using the Larmor formula 
#
Species(
	species_type = 'helium',
	ionization_model = 'tunnel',
	atomic_number = 2,
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	dens_profile = f,
	n_part_per_cell = 1000,
	mass = 1836.0,
	charge = 0.0,
	nb_density = 1.0,
	bc_part_type_west = 'none',
	bc_part_type_east = 'none'
)
Species(
	species_type = 'electron',
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	dens_profile = f,
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
