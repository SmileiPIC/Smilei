import math 


# import sys
# print sys.path

part_per_cell=100
t_sim=50
res=10

position=9
thickness=1


density=10

dx, dy = 15, 40

twopi=2*math.pi

wavelength=0.2

def my_real_func(x,y):
    return 1 if x>position and x< position+thickness and y>5 and y<35 else 0

check_stop_file=True

# 
# def my_function_easy(x,y,leng=1):
#     """ my function 2 """
#     return my_real_func(x,y,leng)
# 
# 
# def my_function_nicer(leng=1):
#     """ my function 1 """
#     return lambda x,y: my_real_func(x,y,leng)      
# 
# 
# class my_function_easy_and_nice():
#     def __init__(self, leng=1):
#         """ my function 3 """
#         self.leng=leng
#         
#     def __call__(self, x, y):
#         return my_real_func(x,y,self.leng)


output_script='my_test.py'


# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

# sim_units: normalisation units for the input data
#            it is used only in the input data & log file
#            codes outputs are always in "normalised" units
#            wavelength = input data are in wavelength-related units
#            normalized = input data are put in code (relativistic) units
#

sim_units = 'wavelength'


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
res_space  = (res, res)  
sim_length = (dx,  dy)

bc_em_type_x  = 'silver-muller'
bc_em_type_y = 'periodic'

# RANDOM seed 
# this is used to randomize the random number generator
random_seed = 0

fieldDump_every = res

avgfieldDump_every =res
ntime_step_avg=res

fieldsToDump = ("Bz", "Rho_electron", "Rho_ion", "Bz_avg")

print_every = 10


Species(
    dens_profile = lambda x, y: my_real_func(x,y),
    species_type = 'ion',
    initPosition_type = 'random',
    initMomentum_type = 'cold',
    ionization_model = 'none',
    n_part_per_cell = part_per_cell,
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1,
    density = density,
    mean_velocity = (0.0,0.0,0.0),
    temperature = 0.0,
    dynamics_type = 'norm',
    time_frozen = 0,
    radiating = False,
    bc_part_type_west  = 'refl',
    bc_part_type_east  = 'refl',
    bc_part_type_south = 'none',
    bc_part_type_north = 'none',
    mvel_x_profile='constant',
    mvel_y_profile='constant',
    mvel_z_profile='constant'
)

Species(
    dens_profile = Species[0].dens_profile,
    species_type = 'electron' ,
    initPosition_type = 'random' ,
    initMomentum_type = 'maxj' ,
    n_part_per_cell = part_per_cell ,
    c_part_max=1.0 ,
    mass = 1.0 ,
    charge = -1 ,
    density = density ,
    mean_velocity = (0.0,0.0,0.0),
    mvel_x_profile='constant',
    mvel_y_profile='constant',
    mvel_z_profile='constant',
    temperature = 0.0001 ,
    dynamics_type = 'norm' ,
    time_frozen = 0.0 ,
    radiating = False ,
    bc_part_type_west  = 'refl' ,
    bc_part_type_east  = 'refl' ,
    bc_part_type_south = 'none' ,
    bc_part_type_north = 'none'
)


def my_func_laser_profile(t,y):
    val = math.exp(-t**2)*math.exp(-(y)**2)
    return val
    
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

# import numpy as np
# large_memory_variable=np.arange(0,200000000,0.9)

def cleanup():
    print "cleanup"

Laser(
    boxSide = 'west' ,
    a0=10 ,
    focus=(position,  dy/2.) ,
    angle=30 ,
    delta=0.0 ,
    time_profile = 'sin2' ,
    double_params = 2 ,
    transv_profile = 'focused' ,
    double_params_transv = 1.0 
)


# DIAG ON SCALARS
# every = number of time-steps between each output
#

DiagScalar(every = 1)


