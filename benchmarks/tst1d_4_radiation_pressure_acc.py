# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField
#
import math

l0 = 2.0*math.pi	# laser wavelength
t0 = l0				# optical cicle
Lsim = 10.*l0		# length of the simulation
Tsim = 40.*t0		# duration of the simulation
resx = 500.			# nb of cells in on laser wavelength
rest = resx/0.95	# time of timestep in one optical cycle (0.95 * CFL)

# plasma slab
def f(x):
    if l0 < x < 2.0*l0:
        return 1.0
    else :
        return 0.0

# wavelength_SI: used by Fred Diags. should be removed
wavelength_SI = 1.e-6

# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
#
dim = '1d3v'

# order of interpolation
interpolation_order = 2 

# SIMULATION BOX : for all space directions (use vector)
# cell_length: length of the cell
# sim_length: length of the simulation in units of the normalization wavelength 
#
cell_length = [l0/resx]
sim_length  = [Lsim]

number_of_patches = [ 8 ]

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
Species(
	species_type = 'ion',
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 10,
	mass = 1836.0,
	charge = 1.0,
	nb_density = trapezoidal(10.,xvacuum=l0,xplateau=l0),
	temperature = [0.],
	bc_part_type_west = 'refl',
	bc_part_type_east = 'refl'
)
Species(
	species_type = 'eon',
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	nb_density = trapezoidal(10.,xvacuum=l0,xplateau=l0),
	temperature = [0.],
	bc_part_type_west = 'refl',
	bc_part_type_east = 'refl'
)

# LASER PROPERTIES
LaserPlanar1D(
	boxSide = 'west',
	a0 = 10.,
    omega = 1.,
    ellipticity = 1.,
    time_envelope = tconstant(),
)

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

every = int(rest/2.)

fieldDump_every	= every
fieldsToDump	= ('Ex','Ey','Ez','Rho_ion','Rho_eon')

DiagScalar(every=every)

# PHASE-SPACE DIAGNOSTICS (new version from DiagParticles)
DiagParticles(
	output = "density",
	every = every,
	species = ["ion"],
	axes = [
		["x",  0.,   Lsim, 200],
		["px", -10., 1000., 200]
	]
)

DiagParticles(
	output = "density",
	every = every,
	species = ["ion"],
	axes = [
		["ekin", 0., 200., 200, "edge_inclusive"]
	]
)


