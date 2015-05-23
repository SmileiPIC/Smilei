import math 

part_per_cell=5
t_sim=30
res=10

position=8
thickness=1


density=0.8

dx, dy = 10, 50

twopi=2*math.pi


def my_function1(pos=0,leng=1,thick=1):
    """ my function 1 """
    return lambda x,y: (math.cos(2*y/twopi)+2.0)*math.exp(-((x/twopi-pos)/leng)**2)/3 if x/twopi<pos else 1 if x/twopi < pos+thick else 0       

def my_function2(codex,codey, pos=0,leng=1,thick=1):
    """ my function 2 """
    x,y=codex/twopi,codey/twopi
    return (math.cos(2*y)+2.0)*math.exp(-((x-pos)/leng)**2)/3 if x<pos else 1 if x < pos+thick else 0


class my_function3():
    def __init__(self, pos=0,leng=1,thick=1):
        """ my function 3 """
        self.pos=pos
        self.leng=leng
        self.thick=thick
        
    def __call__(self, codex, codey):
        x,y=codex/twopi,codey/twopi
        return (math.cos(2*y)+2.0)*math.exp(-((x-self.pos)/self.leng)**2)/3 if x<self.pos else 1 if x < self.pos+self.thick else 0



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

mysim.fieldDump_every = 10

mysim.print_every = 10


myspec1=Species()
myspec1.dens_profile = my_function3(pos=8,leng=1,thick=1)
myspec1.vacuum_length   = ((dx-thickness)/2.0,  dy/4.0) 
myspec1.dens_length_x   = thickness
myspec1.dens_length_y   = dy/2
myspec1.species_type = 'ion'
myspec1.initPosition_type = 'random'
myspec1.initMomentum_type = 'cold'
myspec1.ionization_model = 'none'
myspec1.n_part_per_cell = part_per_cell
myspec1.c_part_max = 1.0
myspec1.mass = 1836.0
myspec1.charge = 1
myspec1.density = density
myspec1.mean_velocity = 0.0
myspec1.temperature = 0.0
myspec1.dynamics_type = 'norm'
myspec1.time_frozen = 0.0
myspec1.radiating = False
myspec1.bc_part_type_west  = 'refl'
myspec1.bc_part_type_east  = 'refl'
myspec1.bc_part_type_south = 'none'
myspec1.bc_part_type_north = 'none'
myspec1.mvel_x_profile='constant'
myspec1.mvel_y_profile='constant'
myspec1.mvel_z_profile='constant'

Species(
dens_profile = lambda x,y : my_function2(x,y,pos=8,leng=1,thick=1),
vacuum_length   = ((dx-thickness)/2.0,  dy/4,0) ,
dens_length_x   = thickness ,
dens_length_y   = dy/2,
species_type = 'electron' ,
initPosition_type = 'random' ,
initMomentum_type = 'maxj' ,
n_part_per_cell = part_per_cell ,
c_part_max=1.0 ,
mass = 1.0 ,
charge = -1 ,
density = density ,
mean_velocity = 0.0 ,
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
Laser(
boxSide = 'west' ,
a0=0.2 ,
focus=(position,  dy/2) ,
angle=20 ,
delta=0.0 ,
time_profile = 'sin2' ,
double_params = 2 ,
transv_profile = my_func_laser_profile ,
double_params_transv = 10.0 
)


# DIAG ON SCALARS
# every = number of time-steps between each output
#

DiagScalar(every = 1)


