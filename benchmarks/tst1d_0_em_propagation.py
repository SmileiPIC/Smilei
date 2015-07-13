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
rest = 102.0		# nb of timestep in 1 optical cycle
resx = 100.0		# nb cells in 1 wavelength

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
sim_length  = [5.0*l0]

# SIMULATION TIME
# timestep: duration of the timestep
# sim_time: duration of the simulation in units of the normalization period 
#
timestep = t0/rest
sim_time = 10.0*t0
 
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
 
# LASER PROPERTIES
# for each laser define:
# a0: maximum amplitude of the laser electric field (in units of the normalization field)
# angle: angle (in degree) at which the laser enters the simulation box
# delta: polarization parameter, (0:y) (1:z) (0.707106781:circ)
# time_profile: string defining the time profile
# double_params: vector of real parameters used by the different time-profiles
# 
Laser(
    a0 = 1.0,
    boxSide = 'west',
    delta = 1.0,
    time_profile = 'constant',
    double_params = [3.0*t0, 0.]
)
 
# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------
#

# log output frequency
print_every = int(rest/2.0)

# SCALAR DIAGNOSTICS
# every       = integer, nb of timesteps between each output
# tmin & tmax = floats, min & max times between which scalars are computed (optional)
# precision   = integer, nb of digits for the outputs (default=10)
DiagScalar(every = 1)

# FIELD DUMPS
# fieldDump_every = integer, nb of timesteps between each output
# fieldsToDump    = ('string'), name of the fields to dump
fieldDump_every = int(rest/2.0)
fieldsToDump = ('Ex','Ey','Ez','By_m','Bz_m');
 
# PROBE DIAGNOSTICS 
DiagProbe(
    every = 1, 
    pos = [0.0]   
)
 
DiagProbe(
    every = 1,
    pos = [0.0],
    pos_first = [10.0],
    number = [10]
)
 
# --------- 
# DUMP INFO
# ---------
#
# dump_step = 2500
# dump_minutes = 1
# random_seed = 13121977
