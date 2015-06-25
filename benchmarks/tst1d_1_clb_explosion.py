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
# resolution
resx = 100.0
rest = 110.0
# plasma length
L = 2.0*math.pi

# step-like density profile
def f(x):
    if L < x < 2.0*L:
        return 1.0
    else :
        return 0.0

# wavelength_SI: used by Fred Diags. should be removed
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
cell_length = [L/resx]
sim_length  = [3.0*L]

# SIMULATION TIME
# timestep: duration of the timestep
# sim_time: duration of the simulation in units of the normalization period 
#
timestep = L/rest
sim_time = 10.0 * math.pi

# PARALLELISATION
clrw = 10

# ELECTROMAGNETIC BOUNDARY CONDITIONS
# bc_em_type_x/y/z : boundary conditions used for EM fields 
#                    periodic = periodic BC (using MPI topology)
#                    silver-muller = injecting/absorbing BC
#                    reflective = consider the ghost-cells as a perfect conductor
#
bc_em_type_x = ['silver-muller']

# RANDOM seed 
# this is used to randomize the random number generator
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
	species_type = "charges",
	initPosition_type = "random",
	initMomentum_type = "cold",
	dens_profile = f,
	n_part_per_cell = 100,
	mass = 1836.0,
	charge = 1.0,
	nb_density = 1.0,
	bc_part_type_west = "stop",
	bc_part_type_east = "stop"
)


# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

fieldDump_every = int(rest/10.)
fieldsToDump=('Ex','Rho_charges')

DiagScalar(
	every = 1
)


