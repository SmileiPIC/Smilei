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
t0 = l0			# optical cicle
Lsim = [40.*l0,5.*l0]	# length of the simulation
Tsim = 30.*t0		# duration of the simulation
resx = 50.		# nb of cells in on laser wavelength
rest = 75.		# time of timestep in one optical cycle 

def fp(x,y):
	if (l0 < x < 11.*l0):
		return 1.
	else:
		return 0.
def fm(x,y):
	if (11.*l0 < x < 39.5*l0):
		return 1.
	else:
		return 0.


# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
#
dim = '2d3v'

# order of interpolation
#
interpolation_order = 2 

# SIMULATION BOX : for all space directions (use vector)
# cell_length: length of the cell
# sim_length: length of the simulation in units of the normalization wavelength 
#
cell_length = [l0/resx,l0/resx]
sim_length  = Lsim

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
bc_em_type_y = ['periodic']

# RANDOM seed 
# this is used to randomize the random number generator
random_seed = 0



# species_type: ion, electron, positron, test ...
# initialization_type: regular, cold or (isotrop) Maxwell−Juettner distribution# n_part_per_cell: number of particle−per−cell
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
	species_type = 'pon1',
	dens_profile = fp,  
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	ionization_model = 'none',
	n_part_per_cell = 2,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	nb_density = 1.0,
	mean_velocity = [0.,0.,0.],
	temperature = 0.001,
	time_frozen = 100000000.0,
	bc_part_type_west  = 'stop',
	bc_part_type_east  = 'refl',
	bc_part_type_south = 'none',
	bc_part_type_north = 'none'
)
Species(
	species_type = 'eon1',
	dens_profile = fp,  
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	ionization_model = 'none',
	n_part_per_cell = 2,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	nb_density = 1.0,
	mean_velocity = [0.,0.,0.],
	temperature = 0.001,
	time_frozen = 0.0,
	bc_part_type_west  = 'stop',
	bc_part_type_east  = 'refl',
	bc_part_type_south = 'none',
	bc_part_type_north = 'none'
)
Species(
	species_type = 'pon2',
	dens_profile = fm,  
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	ionization_model = 'none',
	n_part_per_cell = 2,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	nb_density = 1.0,
	mean_velocity = [-0.5,0.,0.],
	temperature = 0.001,
	time_frozen = 0.0,
	bc_part_type_west  = 'stop',
	bc_part_type_east  = 'refl',
	bc_part_type_south = 'none',
	bc_part_type_north = 'none'
)
Species(
	species_type = 'eon2',
	dens_profile = fm,  
	initPosition_type = 'regular',
	initMomentum_type = 'mj',
	ionization_model = 'none',
	n_part_per_cell = 2,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	nb_density = 1.0,
	mean_velocity = [-0.5,0.,0.],
	temperature = 0.001,
	time_frozen = 0.0,
	bc_part_type_west  = 'stop',
	bc_part_type_east  = 'refl',
	bc_part_type_south = 'none',
	bc_part_type_north = 'none'
)


# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

every = 50
globalEvery = int(rest/2.)

DiagScalar(every=globalEvery)
fieldDump_every = globalEvery
fieldsToDump = ('Ex','Ey','Ez','Bx','By','Bz','Rho_pon1','Rho_eon1','Rho_pon2','Rho_eon2')