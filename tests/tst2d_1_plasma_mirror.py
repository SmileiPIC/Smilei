mysim=smilei()
myspec1=species(mysim)
myspec2=species(mysim)
laser1=laser(mysim)

# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

# sim_units: normalisation units for the input data
#            it is used only in the input data & log file
#            codes outputs are always in "normalised" units
#            wavelength = input data are in wavelength-related units
#            normalized = input data are put in code (relativistic) units
#

mysim.sim_units = 'wavelength'


# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
#
mysim.dim = '2d3v'

# order of interpolation
mysim.interpolation_order = 2 

# SIMULATION TIME
# res_time: temporal resolution (integer = number of time−steps within one normalization period)
# sim_time: duration of the simulation in units of the normalization period 
#
mysim.res_time = 30
mysim.sim_time = 40

# SIMULATION BOX : for all space directions (use vector)
# res_space: spatial resolution (vector of integer = number of cells in one normalization wavelength )
# sim_length: length of the simulation in units of the normalization wavelength 
#
mysim.res_space  = (20,  20)  
mysim.sim_length = (20,  50)

mysim.bc_em_type_long  = 'silver-muller'
mysim.bc_em_type_trans = 'periodic'

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


myspec1.dens_profile = 'constant'
myspec1.vacuum_length   = (10.0,  0.0) 
myspec1.dens_length_x   = 4.0  
myspec1.dens_length_y   = 50.0
myspec1.species_type = 'ion'
myspec1.initPosition_type = 'random'
myspec1.initMomentum_type = 'cold'
myspec1.ionization_model = 'none'
myspec1.n_part_per_cell = 5
myspec1.c_part_max = 1.0
myspec1.mass = 1836.0
myspec1.charge = 1
myspec1.density = 2.0
myspec1.mean_velocity = 0.0
myspec1.temperature = 0.0
myspec1.dynamics_type = 'norm'
myspec1.time_frozen = 30.0
myspec1.radiating = False
myspec1.bc_part_type_west  = 'refl'
myspec1.bc_part_type_east  = 'refl'
myspec1.bc_part_type_south = 'none'
myspec1.bc_part_type_north = 'none'


myspec2.dens_profile = 'constant'
myspec2.vacuum_length   = (10.0,  0.0) 
myspec2.dens_length_x   = 4.0  
myspec2.dens_length_y   = 50.0
myspec2.species_type = 'electron'
myspec2.initPosition_type = 'random'
myspec2.initMomentum_type = 'maxj'
myspec2.n_part_per_cell=5
myspec2.c_part_max=1.0
myspec2.mass = 1.0
myspec2.charge = -1
myspec2.density = 2.0
myspec2.mean_velocity = 0.0
myspec2.temperature = 0.01
myspec2.dynamics_type = 'norm'
myspec2.time_frozen = 0.0
myspec2.radiating = False
myspec2.bc_part_type_west  = 'refl'
myspec2.bc_part_type_east  = 'refl'
myspec2.bc_part_type_south = 'none'
myspec2.bc_part_type_north = 'none'



# ----------------
# LASER PROPERTIES
# ----------------
#
# for each laser define:
# a0: maximum amplitude of the laser electric field (in units of the normalization field)
# angle: angle (in degree) at which the laser enters the simulation box
# delta: polarization parameter, (0:y) (1:z) (0.707106781:circ)
# time_profile: string defining the time profile
# double_params: vector of real parameters used by the different time-profiles
#

laser1.boxSide = 'west'
laser1.a0=0.1
laser1.focus=(10.0,  25.0)
laser1.angle=20.0 
laser1.delta=0.0              
laser1.time_profile = 'sin2'
laser1.double_params = 5
laser1.transv_profile = 'focused'
laser1.double_params_transv = 5.0 

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

# print_every (on screen text output) 
# print_every = 60


# DIAG ON SCALARS
# every = number of time-steps between each output
#
scalar1=diag_scalar(mysim)
scalar1.every = 10

mysim.fieldDump_every = 15

