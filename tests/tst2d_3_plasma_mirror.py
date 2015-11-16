# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField

print "----------------------- RANK -----------------------",smilei_mpi_rank

import math 

part_per_cell=5
t_sim=40
thickness=1

res=8

dx, dy = 20, 50

twopi=2*math.pi

def profile_dens(x,y):
    if x < twopi*(dx-thickness)/2.0 or x  > twopi*(dx+thickness)/2.0:
        return 0.0
    else :
        return 1.0



# sim_units: normalisation units for the input data
#            it is used only in the input data & log file
#            codes outputs are always in "normalised" units
#            wavelength = input data are in wavelength-related units
#            normalized = input data are put in code (relativistic) units
#

sim_units = 'wavelength'
raise Exception("TODO : modify the input file without sim_units")

# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
#
dim = '2d3v'

# order of interpolation
interpolation_order = 2 

# SIMULATION TIME
# res_time: temporal resolution integer  number of time-steps within one normalization period
# sim_time: duration of the simulation in units of the normalization period 
#
res_time = 1.5*res
sim_time = t_sim

# SIMULATION BOX : for all space directions (use vector)
# res_space: spatial resolution (vector of integer = number of cells in one normalization wavelength )
# sim_length: length of the simulation in units of the normalization wavelength 
#
res_space  = [res, res]
sim_length = [dx,  dy]

bc_em_type_x  = 'silver-muller'
bc_em_type_y = 'periodic'

# RANDOM seed 
# this is used to randomize the random number generator
random_seed = 0

fieldDump_every = 24

avgfieldDump_every =24
ntime_step_avg=24

fieldsToDump = ("Bz", "Rho_electron", "Rho_ion", "Bz_avg")


Species(
dens_profile = profile_dens,
vacuum_length   = [(dx-thickness)/2.0,  0.0] ,
dens_length_x   = thickness,
dens_length_y   = dy,
species_type = 'ion',
initPosition_type = 'random',
initMomentum_type = 'cold',
ionization_model = 'none',
n_part_per_cell = part_per_cell,
c_part_max = 1.0,
mass = 1836.0,
charge = 1,
density = 2.0,
mean_velocity = [0.0,0.0,0.0],
temperature = 0.0,
dynamics_type = 'norm',
time_frozen = t_sim,
radiating = False,
bc_part_type_west  = 'refl',
bc_part_type_east  = 'refl',
bc_part_type_south = 'none',
bc_part_type_north = 'none'
)


Species(
dens_profile = profile_dens,
vacuum_length   = [(dx-thickness)/2.0,  0.0] ,
dens_length_x   = thickness ,
dens_length_y   = dy,
species_type = 'electron' ,
initPosition_type = 'random' ,
initMomentum_type = 'maxj' ,
n_part_per_cell = part_per_cell ,
c_part_max=1.0 ,
mass = 1.0 ,
charge = -1 ,
density = 2.0 ,
mean_velocity =  [0.0,0.0,0.0],
temperature = 0.0001 ,
dynamics_type = 'norm' ,
time_frozen = 0.0 ,
radiating = False ,
bc_part_type_west  = 'refl' ,
bc_part_type_east  = 'refl' ,
bc_part_type_south = 'none' ,
bc_part_type_north = 'none'
)


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
Laser(
boxSide = 'west' ,
a0=0.1 ,
focus=[10.0,  25.0] ,
angle=20.0 ,
delta=0.0 ,
time_profile = 'sin2' ,
double_params = 5 ,
transv_profile = 'focused' ,
double_params_transv = 5.0 
)

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

# print_every (on screen text output) 
# print_every = 60


# DIAG ON SCALARS
# every = number of time-steps between each output
#

DiagScalar(every = 10)


