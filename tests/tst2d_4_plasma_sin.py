import math 

part_per_cell=1
t_sim=40
thickness=0.5

res=4

dx, dy = 20, 50

twopi=2*math.pi

def my_func_density(codex,codey):
    x,y=codex/twopi,codey/twopi

    valx=math.exp(-((x-dx/2))**2)
    valy=math.cos(y/2)+1.0
        
    return valx*valy


mysim=Smilei()

mysim.output_script='my_test.py'


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
# res_time: temporal resolution integer  number of time-steps within one normalization period
# sim_time: duration of the simulation in units of the normalization period 
#
mysim.res_time = 1.5*res
mysim.sim_time = t_sim

# SIMULATION BOX : for all space directions (use vector)
# res_space: spatial resolution (vector of integer = number of cells in one normalization wavelength )
# sim_length: length of the simulation in units of the normalization wavelength 
#
mysim.res_space  = (res, res)  
mysim.sim_length = (dx,  dy)

mysim.bc_em_type_long  = 'silver-muller'
mysim.bc_em_type_trans = 'periodic'

# RANDOM seed 
# this is used to randomize the random number generator
mysim.random_seed = 0

mysim.fieldDump_every = 20

mysim.print_every = 10



myspec1=Species()
myspec1.dens_profile = my_func_density
myspec1.vacuum_length   = ((dx-thickness)/2.0,  0.0) 
myspec1.dens_length_x   = thickness
myspec1.dens_length_y   = dy
myspec1.species_type = 'ion'
myspec1.initPosition_type = 'random'
myspec1.initMomentum_type = 'cold'
myspec1.ionization_model = 'none'
myspec1.n_part_per_cell = part_per_cell
myspec1.c_part_max = 1.0
myspec1.mass = 1836.0
myspec1.charge = 1
myspec1.density = 2.0
myspec1.mean_velocity = 0.0
myspec1.temperature = 0.0
myspec1.dynamics_type = 'norm'
myspec1.time_frozen = t_sim
myspec1.radiating = False
myspec1.bc_part_type_west  = 'refl'
myspec1.bc_part_type_east  = 'refl'
myspec1.bc_part_type_south = 'none'
myspec1.bc_part_type_north = 'none'



Species(
dens_profile = my_func_density,
vacuum_length   = ((dx-thickness)/2.0,  0.0) ,
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
mean_velocity = 0.0 ,
temperature = 0.0001 ,
dynamics_type = 'norm' ,
time_frozen = 0.0 ,
radiating = False ,
bc_part_type_west  = 'refl' ,
bc_part_type_east  = 'refl' ,
bc_part_type_south = 'none' ,
bc_part_type_north = 'none'
)



def my_func_laser_time(codet) :
    t=codet*twopi
    return math.exp(-(t-t_sim/(twopi*2.0))**2)
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
focus=(10.0,  25.0) ,
angle=20.0 ,
delta=0.0 ,

time_profile = my_func_laser_time ,
double_params = 5 ,
transv_profile = 'focused' ,
double_params_transv = 5.0 
)


# DIAG ON SCALARS
# every = number of time-steps between each output
#

DiagScalar(every = 1)


