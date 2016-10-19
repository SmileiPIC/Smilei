"""@package pyinit
    Definition of Smilei components
"""

class SmileiComponentType(type):
    """Metaclass to all Smilei components"""
    
    # Constructor of classes
    def __init__(self, name, bases, attrs):
        self._list = []
        self.verify = True
        # Run standard metaclass init
        super(SmileiComponentType, self).__init__(name, bases, attrs)
    
    # Functions to define the iterator
    def __iter__(self):
        self.current = 0
        return self
    def next(self):
        if self.current >= len(self._list):
            raise StopIteration
        self.current += 1
        return self._list[self.current - 1]
    
    # Function to return one given instance, for example DiagParticles[0]
    # Special case: species can also be indexed by their name: Species["ion1"]
    def __getitem__(self, key):
        if self.__name__ == "Species" and type(key) is str:
            for obj in self._list:
                if obj.species_type == key:
                    return obj
        else:
            return self._list[key]
    
    # Function to return the number of instances, for example len(Species)
    def __len__(self):
        return len(self._list)
    
    # Function to display the content of the component
    def __repr__(self):
        if len(self._list)==0:
            return "<Empty list of "+self.__name__+">"
        else:
            l = []
            for obj in self._list: l.append(str(obj))
            return "["+", ".join(l)+"]"


class SmileiComponent(object):
    """Smilei component generic class"""
    __metaclass__ = SmileiComponentType
    
    # This constructor is used always for all child classes
    def __init__(self, **kwargs):
        if kwargs is not None: # add all kwargs as internal class variables
            for key, value in kwargs.iteritems():
                if key=="_list":
                    print "Python warning: in "+type(self).__name__+": cannot have argument named '_list'. Discarding."
                else:
                    setattr(self, key, value)
        type(self)._list.append(self) # add the current object to the static list "list"


class Species(SmileiComponent):
    """Species parameters"""
    species_type = None
    initPosition_type = None
    initMomentum_type = ""
    n_part_per_cell = None
    c_part_max = 1.0
    charge_density = None
    nb_density = None
    mean_velocity = [0.]
    temperature = [1e-10]
    dynamics_type = "norm"
    time_frozen = 0.0
    radiating = False
    bc_part_type_xmin = None
    bc_part_type_xmax = None
    bc_part_type_ymax = None
    bc_part_type_ymin = None
    ionization_model = "none"
    atomic_number = None
    isTest = False
    track_every = 0

class Laser(SmileiComponent):
    """Laser parameters"""
    boxSide = None
    a0 = None
    omega0 = 1.
    delta = 1.
    tchirp = 0.
    focus = []
    angle = 0.
    delay = 0.
    time_profile = None
    int_params = []
    double_params = []
    transv_profile = None
    int_params_transv = []
    double_params_transv = []

class Collisions(SmileiComponent):
    """Collisions parameters"""
    species1 = None
    species2 = None
    coulomb_log = 0.
    debug_every = 0


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    every = 0
    time_range = [None, None]
    number = []
    pos = []
    pos_first = []
    pos_second = []
    pos_third = []
    fields = []

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    output = None
    every = None
    time_average = 1
    species = None
    axes = []

class DiagPhase(SmileiComponent):
    """Diagnostic phase"""
    every=None
    first=[]
    second=[]
    time_range = []
    deflate = 0
    pass

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    every = None
    time_range = []
    precision = 10
    vars = []

# external fields
class ExtField(SmileiComponent):
    """External Field"""
    field = []
    profile = None

# external current (antenna)
class Antenna(SmileiComponent):
    """Antenna"""
    field = []
    time_profile  = None
    space_profile = None

# Particle wall
class PartWall(SmileiComponent):
    """Particle Wall"""
    kind = "none"
    x = None
    y = None
    z = None

# default simulation values
output_dir = None
smilei_mpi_rank = 0
smilei_mpi_size = 1
smilei_rand_max = 2**31-1
dump_step = 0
dump_minutes = 0.0
exit_after_dump = True
restart = False
dump_file_sequence = 2
dump_deflate = 0
restart_dir = None
sim_units = ""
wavelength_SI = 0.
dim = ""
interpolation_order = 2
timestep = None
cell_length = []
sim_time = None
sim_length = []
maxwell_sol = 'Yee'
bc_em_type_x = []
bc_em_type_y = []
time_fields_frozen = 0.0
nspace_win_x = 0
t_move_win = 0.0
vx_win = 1.
clrw = 1
every = 0
number_of_procs = [None]
print_every = None
fieldDump_every = None
fieldsToDump = []
avgfieldDump_every = None
ntime_step_avg = 0
time_fields_frozen = 0.
random_seed = None
# Some predefined profiles (see doc)

def constant(value, xvacuum=0., yvacuum=0.):
    global dim, sim_length
    if dim == "1d3v": return lambda x: value #if x>=xvacuum else 0.
    if dim == "2d3v": return lambda x,y: value #if (x>=xvacuum and y>=yvacuum) else 0.

def trapezoidal(max,
                xvacuum=0., xplateau=None, xslope1=0., xslope2=0.,
                yvacuum=0., yplateau=None, yslope1=0., yslope2=0. ):
    global dim, sim_length
    if not dim or not sim_length:
        raise Exception("trapezoidal profile has been defined before `dim` or `sim_length`")
    if len(sim_length)>0 and xplateau is None: xplateau = sim_length[0]-xvacuum
    if len(sim_length)>1 and yplateau is None: yplateau = sim_length[1]-yvacuum
    def fx(x):
        # vacuum region
        if x < xvacuum: return 0.
        # linearly increasing density
        elif x < xvacuum+xslope1: return max*(x-xvacuum) / xslope1
        # density plateau
        elif x < xvacuum+xslope1+xplateau: return max
        # linearly decreasing density
        elif x < xvacuum+xslope1+xplateau+xslope2:
            return max*(1. - ( x - (xvacuum+xslope1+xslope2) ) / xslope2)
        # beyond the plasma
        else: return 0.0
    if dim == "1d3v": return fx
    def fy(y):
        # vacuum region
        if y < yvacuum: return 0.
        # linearly increasing density
        elif y < yvacuum+yslope1: return (y-yvacuum) / yslope1
        # density plateau
        elif y < yvacuum+yslope1+yplateau: return 1.
        # linearly decreasing density
        elif y < yvacuum+yslope1+yplateau+yslope2:
            return 1. - ( y - (yvacuum+yslope1+yslope2) ) / yslope2
        # beyond
        else: return 0.0
    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def gaussian(max,
             xvacuum=0., xlength=None, xfwhm=None, xcenter=None, xorder=2,
             yvacuum=0., ylength=None, yfwhm=None, ycenter=None, yorder=2 ):
    import math
    global dim, sim_length
    if not dim or not sim_length:
        raise Exception("gaussian profile has been defined before `dim` or `sim_length`")
    if len(sim_length)>0:
        if xlength is None: xlength = sim_length[0]-xvacuum
        if xfwhm   is None: xfwhm   = xlength/3.
        if xcenter is None: xcenter = xvacuum + xlength/2.
    if len(sim_length)>1: 
        if ylength is None:ylength = sim_length[1]-yvacuum
        if yfwhm   is None:yfwhm   = ylength/3.
        if ycenter is None:ycenter = yvacuum + ylength/2.
    sigmax = (0.5*xfwhm)**xorder/math.log(2.0)
    def fx(x):
        # vacuum region
        if x < xvacuum: return 0.
        # gaussian
        elif x < xvacuum+xlength: return max*math.exp( -(x-xcenter)**xorder / sigmax )
        # beyond
        else: return 0.0
    if dim == "1d3v": return fx
    sigmay = (0.5*yfwhm)**yorder/math.log(2.0)
    def fy(y):
        if yorder == 0: return 1.
        # vacuum region
        if y < yvacuum: return 0.
        # gaussian
        elif y < yvacuum+ylength: return math.exp( -(y-ycenter)**yorder / sigmay )
        # beyond
        else: return 0.0
    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def polygonal(xpoints=[], xvalues=[]):
    global dim, sim_length
    if not dim or not sim_length:
        raise Exception("polygonal profile has been defined before `dim` or `sim_length`")
    if len(xpoints)!=len(xvalues):
        raise Exception("polygonal profile requires as many points as values")
    if len(sim_length)>0 and len(xpoints)==0:
        xpoints = [0., sim_length[0]]
        xvalues = [1., 1.]
    N = len(xpoints)
    xpoints = [float(x) for x in xpoints]
    xvalues = [float(x) for x in xvalues]
    xslopes = [0. for i in range(1,N)]
    for i in range(1,N):
        if xpoints[i] == xpoints[i-1]: continue
        xslopes[i-1] = (xvalues[i]-xvalues[i-1])/(xpoints[i]-xpoints[i-1])
    def f(x,y=0.):
        if x < xpoints[0]: return 0.0;
        for i in range(1,N):
            if x<xpoints[i]: return xvalues[i-1] + xslopes[i-1] * ( x-xpoints[i-1] )
        return 0.
    return f

def cosine(base,
           xamplitude=1., xvacuum=0., xlength=None, xphi=0., xnumber=1,
           yamplitude=1., yvacuum=0., ylength=None, yphi=0., ynumber=1):
    import math
    global dim, sim_length
    if not dim or not sim_length:
        raise Exception("cosine profile has been defined before `dim` or `sim_length`")
    
    if len(sim_length)>0 and xlength is None: xlength = sim_length[0]-xvacuum
    if len(sim_length)>1 and ylength is None: ylength = sim_length[1]-yvacuum
    
    def fx(x):
        #vacuum region
        if x < xvacuum: return 0.
        # profile region
        elif x < xvacuum+xlength:
            return base + xamplitude * math.cos(xphi + 2.*math.pi * xnumber * (x-xvacuum)/xlength)
        # beyond
        else: return 0.
    if dim == "1d3v": return fx
    def fy(y):
        #vacuum region
        if y < yvacuum: return 0.
        # profile region
        elif y < yvacuum+ylength:
            return base + yamplitude * math.cos(yyphi + 2.*math.pi * ynumber * (y-yvacuum)/ylength)
        # beyond
        else: return 0.

    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def tconstant(start=0.):
    return lambda t: 1. if t>=start else 0.

def ttrapezoidal(start=0., plateau=None, slope1=0., slope2=0.):
    global sim_time
    if not sim_time:
        raise Exception("ttrapezoidal profile has been defined before `sim_time`")
    if plateau is None: plateau = sim_time - start
    def ft(t):
        if t < start: return 0.
        elif t < start+slope1: return (t-start) / slope1
        elif t < start+slope1+plateau: return 1.
        elif t < start+slope1+plateau+slope2:
            return 1. - ( t - (start+slope1+slope2) ) / slope2
        else: return 0.0
    return ft

def tgaussian(start=0., duration=None, fwhm=None, center=None, order=2):
    import math
    global sim_time
    if not sim_time:
        raise Exception("tgaussian profile has been defined before `sim_time`")
    if duration is None: duration = sim_time-start
    if fwhm     is None: fwhm     = duration/3.
    if center   is None: center   = start + duration/2.
    sigma = (0.5*fwhm)**order/math.log(2.0)
    def ft(t):
        if t < start: return 0.
        elif t < start+duration: return math.exp( -(t-center)**order / sigma )
        else: return 0.0

def tpolygonal(points=[], values=[]):
    global sim_time
    if not sim_time:
        raise Exception("tpolygonal profile has been defined before `sim_time`")
    if len(points)==0:
        points = [0., sim_time]
        values = [1., 1.]
    N = len(points)
    points = [float(x) for x in points]
    values = [float(x) for x in values]
    slopes = [0. for i in range(1,N)]
    for i in range(1,N):
        if points[i] == points[i-1]: continue
        slopes[i-1] = (values[i]-values[i-1])/(points[i]-points[i-1])
    def f(t):
        if t < points[0]: return 0.0;
        for i in range(1,N):
            if t<points[i]: return values[i-1] + slopes[i-1] * ( t-points[i-1] )
        return 0.
    return f

def tcosine(base=0., amplitude=1., start=0., duration=None, phi=0., freq=1.):
    import math
    global sim_time
    if not sim_time:
        raise Exception("tcosine profile has been defined before `sim_time`")
    if duration is None: duration = sim_time-start
    def f(t):
        if t < start: return 0.
        elif t < start+duration:
            return base + amplitude * math.cos(phi + freq*(t-start))
        else: return 0.
    return f
############### BEGIN USER NAMELISTS/COMMANDS ###############
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
L  = 1.03			# wavelength=simulation box length
dn = 0.001			# amplitude of the perturbation



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
cell_length = [0.01]
sim_length  = [L]

# SIMULATION TIME
# timestep: duration of the timestep
# sim_time: duration of the simulation in units of the normalization period 
#
timestep = 0.0095
sim_time = 50.
 
# ELECTROMAGNETIC BOUNDARY CONDITIONS
# bc_em_type_x/y/z : boundary conditions used for EM fields 
#                    periodic = periodic BC (using MPI topology)
#                    silver-muller = injecting/absorbing BC
#                    reflective = consider the ghost-cells as a perfect conductor
#
bc_em_type_x = ['periodic']
 
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
	species_type = "ion",
	initPosition_type = "regular",
	initMomentum_type = "cold",
	n_part_per_cell = 10,
	mass = 1836.0,
	charge = 1.0,
	nb_density = 1.,
	time_frozen = 10000.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
)
Species(
	species_type = "eon1",
	initPosition_type = "regular",
	initMomentum_type = "cold",
	n_part_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	nb_density = cosine(0.5,xamplitude=dn,xlength=L),
	mean_velocity = [-0.1,0.0,0.0],
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
)
Species(
	species_type = "eon2",
	initPosition_type = "regular",
	initMomentum_type = "cold",
	n_part_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	nb_density = cosine(0.5,xamplitude=dn,xlength=L),
	mean_velocity = [0.1,0.0,0.0],
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none"
)


# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

#every = 1000000000000
every = 100

# SCALAR DIAGNOSTICS
# every       = integer, nb of timesteps between each output
# tmin & tmax = floats, min & max times between which scalars are computed (optional)
# precision   = integer, nb of digits for the outputs (default=10)
DiagScalar(every = every, vars=['Utot','Ubal_norm','Uelm','Ukin'])	


# FIELD DUMPS
# fieldDump_every = integer, nb of timesteps between each output
# fieldsToDump    = ('string'), name of the fields to dump
fieldDump_every = every
fieldsToDump = ('Ex','Ey','Ez','By_m','Bz_m');


# PHASE-SPACE DIAGNOSTICS (older version working with TPUPMC/SmileiQt.py)
# kind of projection: 1D) xPx xPy xPz xLor PxPy PxPz PyPz
# kind of projection: 2D) xPx xPy xPz xLor yPx yPy yPz yLor PxPy PxPz PyPz
#
DiagPhase(
	every	= every,
 	species = ['eon1','eon2'],
 	kind    = ['xpx'],
 	first = [0., 0., 50],
 	second = [-0.4, 0.4, 100],
 	deflate=5
)

# PHASE-SPACE DIAGNOSTICS (new version from DiagParticles)
DiagParticles(
	output = "density",
	every = every,
	species = ["eon1","eon2"],
	axes = [
		["x", 0., L, 50],
		["px", -0.4, 0.4, 100]
	]
)

################ END USER NAMELISTS/COMMANDS  ################
"""@package pycontrol
    here we check if the namelist is clean and without errors
"""

import gc 
gc.collect()

import os

def _smilei_check():
    """Do checks over the script"""
    # Verify classes were not overriden
    for CheckClassName,CheckClass in {"SmileiComponent":SmileiComponent,"Species":Species,
            "Laser":Laser,"Collisions":Collisions,"DiagProbe":DiagProbe,"DiagParticles":DiagParticles,
            "DiagScalar":DiagScalar,"DiagPhase":DiagPhase,"ExtField":ExtField}.iteritems():
        try:
            if not CheckClass.verify: raise
        except:
            raise Exception("ERROR in the namelist: it seems that the name `"+CheckClassName+"` has been overriden")

    if smilei_mpi_rank == 0 and output_dir:
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except OSError as exception:
                raise Exception("ERROR in the namelist: output_dir "+output_dir+" does not exists and cannot be created")
        elif not os.path.isdir(output_dir):
                raise Exception("ERROR in the namelist: output_dir "+output_dir+" exists and is not a dir")

    if restart and restart_dir:
        if not os.path.isdir(restart_dir):
            raise Exception("ERROR in the namelist: restart_dir "+restart_dir+" is not a dir")



# this function will be called after initialising the simulation, just before entering the time loop
# if it returns false, the code will call a Py_Finalize();
def _keep_python_running():
    for las in Laser:
        for prof in (las.time_profile, las.transv_profile):
            if callable(prof): return True
    for ant in Antenna:
        if callable(ant.time_profile): return True

    if not nspace_win_x == 0:
        return True

    return False

# Prevent creating new components (by mistake)
def _noNewComponents(cls, *args, **kwargs):
    print "Please do not create a new "+cls.__name__
    return None
SmileiComponent.__new__ = staticmethod(_noNewComponents)



