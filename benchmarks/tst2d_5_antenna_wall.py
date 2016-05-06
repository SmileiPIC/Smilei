# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField
#
import math

l0 = 2.0*math.pi		# laser wavelength
t0 = l0					# optical cicle
Lsim = [20.*l0,20.*l0]	# length of the simulation
Tsim = 20.*t0			# duration of the simulation
resx = 5.				# nb of cells in on laser wavelength
rest = 8.				# time of timestep in one optical cycle 

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

number_of_patches = [ 4, 4  ]

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
bc_em_type_y = ['silver-muller']

# RANDOM seed 
# this is used to randomize the random number generator
random_seed = 0

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

globalEvery = 20
    
Species(
	initPosition_type = 'random',
	initMomentum_type = 'cold',
	n_part_per_cell = 2,
	mass = 1.0,
	charge = -1.0,
	nb_density = trapezoidal(1.0,xvacuum=1.*l0,xplateau=4.*l0,yvacuum=5.*l0,yplateau=10.*l0),
	bc_part_type_west  = 'refl',
	bc_part_type_east  = 'refl',
	bc_part_type_south = 'none',
	bc_part_type_north = 'none',
	mean_velocity=[0.9,0.01,0]
)

Species(
	initPosition_type = 'random',
	initMomentum_type = 'cold',
	n_part_per_cell = 2,
	mass = 1.0,
	charge = 1.0,
	nb_density = trapezoidal(1.0,xvacuum=1.*l0,xplateau=4.*l0,yvacuum=5.*l0,yplateau=10.*l0),
	bc_part_type_west  = 'refl',
	bc_part_type_east  = 'refl',
	bc_part_type_south = 'none',
	bc_part_type_north = 'none',
	mean_velocity=[0.9,0.01,0]
)

Antenna(
    field='Jz',
    time_profile= lambda t: math.sin(2*t/t0),
    space_profile=gaussian(0.2, xfwhm=l0, yfwhm=l0, xcenter=sim_length[0]*0.6, ycenter=sim_length[1]*0.5)
)

PartWall (
    x= 15.0*l0,
    kind="refl"
)

# DIAG ON SCALARS
# every = number of time-steps between each output
#
DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ez','Rho_species0','Rho_species1']
)

