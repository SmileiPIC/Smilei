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
# initialization_type: regular, cold or (isotrop) Maxwell Juettner distribution
# n_part_per_cell: number of particle per cell
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
	initMomentum_type = 'cold',
	n_part_per_cell = 10,
	mass = 1836.0,
	charge = 1.0,
	nb_density = 10.,
	temperature = 0.,
	bc_part_type_west = 'refl',
	bc_part_type_east = 'refl'
)
Species(
	species_type = 'eon',
	dens_profile = f,
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	nb_density = 10.,
	temperature = 0.,
	bc_part_type_west = 'refl',
	bc_part_type_east = 'refl'
)

# LASER PROPERTIES
#
# for each laser define:
# a0: maximum amplitude of the laser electric field (in units of the normalization field)
# angle: angle (in degree) at which the laser enters the simulation box
# delta: polarization parameter, (0:y) (1:z) (0.707106781:circ)
# time_profile: string defining the time profile
# double_params: vector of real parameters used by the different time-profiles
#
Laser(
	boxSide = 'west',
	a0=10.0,
	delta=0.707106781,                 
	time_profile = 'constant'
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
		["px", -10., 500., 200]
	]
)

DiagParticles(
	output = "density",
	every = every,
	species = ["ion"],
	axes = [
		["ekin", 0., 100., 200, "edge_inclusive"]
	]
)


