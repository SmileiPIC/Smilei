"""@package pyinit
    Definition of Smilei components
"""

# Since the pytohn interpreter grabs key keyboards,
# we have to filter the ctrl-c kill command:
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

class SmileiComponentType(type):
    """Metaclass to all Smilei components"""
    
    # Constructor of classes
    def __init__(self, name, bases, attrs):
        self.list = []
        self.verify = True
        self.current = 0
        # Run standard metaclass init
        super(SmileiComponentType, self).__init__(name, bases, attrs)
    
    # Functions to define the iterator
    def __iter__(self):
        return self
    def next(self):
        if self.current >= len(self.list):
            raise StopIteration
        self.current += 1
        return self.list[self.current - 1]
    
    # Function to return one given instance, for example DiagParticles[0]
    # Special case: species can also be indexed by their name: Species["ion1"]
    def __getitem__(self, key):
        if self.__name__ == "Species" and type(key) is str:
            for obj in self.list:
                if obj.species_type == key:
                    return obj
        else:
            return self.list[key]
    
    # Function to return the number of instances, for example len(Species)
    def __len__(self):
        return len(self.list)
    
    # Function to display the content of the component
    def __repr__(self):
        if len(self.list)==0:
            return "<Empty list of "+self.__name__+">"
        else:
            l = []
            for obj in self.list: l.append(str(obj))
            return "["+", ".join(l)+"]"


class SmileiComponent(object):
    """Smilei component generic class"""
    __metaclass__ = SmileiComponentType
    
    # This constructor is used always for all child classes
    def __init__(self, **kwargs):
        if kwargs is not None: # add all kwargs as internal class variables
            for key, value in kwargs.iteritems():
                if key=="list":
                    print "Python warning: in "+type(self).__name__+": cannot have argument named 'list'. Discarding."
                else:
                    setattr(self, key, value)
        type(self).list.append(self) # add the current object to the static list "list"


class Species(SmileiComponent):
    """Species parameters"""
    species_type = None
    initPosition_type = None
    initMomentum_type = ""
    n_part_per_cell = None
    c_part_max = 1.0
    charge_density = None
    nb_density = None
    density = None
    mean_velocity = None
    temperature = None
    dynamics_type = "norm"
    time_frozen = 0.0
    radiating = False
    bc_part_type_west = None
    bc_part_type_east = None
    bc_part_type_north = None
    bc_part_type_south = None
    ionization_model = "none"
    atomic_number = None
    vacuum_length = []
    for prefix in ["dens","mvel_x","mvel_y","mvel_z","temp_x","temp_y","temp_z"]:
        exec prefix+"_profile = None"
        for suffix in ["length_x","length_y","length_z","dbl_params","int_params"]:
            exec prefix+"_"+suffix+" = []"

class Laser(SmileiComponent):
    """Laser parameters"""
    pass

class Collisions(SmileiComponent):
    """Collisions parameters"""
    debug_every = 0


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    pass

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    time_average = 1

class DiagPhase(SmileiComponent):
    """Diagnostic phase"""
    pass

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    every = None
    time_range = [None]
    precision = 10
    vars = [None]

# external fields
class ExtField(SmileiComponent):
    """External Field"""
    pass

# default simulation values
output_script = "smilei.py"
dump_step = 0
dump_minutes = 0.0
exit_after_dump = True
restart = False
check_stop_file = False
dump_file_sequence = 2
sim_units = ""
wavelength_SI = 0.
dim = ""
interpolation_order = None
res_time = None
res_space = [None]
timestep = None
cell_length = [None]
sim_time = None
sim_length = [None]
bc_em_type_long = None
bc_em_type_trans = None
nspace_win_x = 0
t_move_win = 0.0
vx_win = 1.
clrw = 1
every = 0
number_of_procs = [None]
print_every = None
fieldDump_every = 0
fieldsToDump = [None]
avgfieldDump_every = None
ntime_step_avg = 0
particleDump_every = None # for backwards-compatibility



# Some predefined profiles

def constant(value):
    return lambda x: value

def trapezoidal(max=1.,
                xvacuum=0., xplateau=None, xslope1=0., xslope2=0.,
                yvacuum=0., yplateau=None, yslope1=0., yslope2=0. ):
    global dim, sim_length
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
    if dim == "1d3v": return fx
    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def gaussian(max=1.,
             xvacuum=0., xlength=None, xfwhm=None, xcenter=None, xorder=2,
             yvacuum=0., ylength=None, yfwhm=None, ycenter=None, yorder=2 ):
    import math
    global dim, sim_length
    if len(sim_length)>0:
        if xlength is None: xlength = sim_length[0]-xvacuum
        if xfwhm   is None: xfwhm   = xlength/3.
        if xcenter is None: xcenter = xvacuum + xlength/2.
    if len(sim_length)>1: 
        if ylength is None:ylength = sim_length[1]-yvacuum
        if yfwhm   is None:yfwhm   = ylength/3.
        if ycenter is None:ycenter = yvacuum + ylength/2.
    def fx(x):
        sigmaN = (0.5*xfwhm)**xorder/math.log(2.0)
        # vacuum region
        if x < xvacuum: return 0.
        # gaussian
        elif x < xvacuum+xlength: return max*math.exp( -(x-xcenter)**xorder / sigmaN )
        # beyond
        else: return 0.0
    def fy(y):
        if yorder == 0: return 1.
        sigmaN = (0.5*yfwhm)**yorder/math.log(2.0)
        # vacuum region
        if y < yvacuum: return 0.
        # gaussian
        elif y < yvacuum+ylength: return math.exp( -(y-ycenter)**yorder / sigmaN )
        # beyond
        else: return 0.0
    if dim == "1d3v": return fx
    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def polygonal(xvacuum=0., xpoints=[], xvalues=[]):
    global dim, sim_length
    if len(xpoints)!=len(xvalues):
        raise Exception("polygonal profile requires as many points as values")
    if len(sim_length)>0 and len(xpoints)==0:
        xpoints = [0., sim_length[0]]
        xvalues = [1., 1.]
    N = len(xpoints)
    def f(x,y=0.):
        # vacuum region
        if x < xvacuum: return 0.0;
        # polygon region (defined over N segments)
        elif x < xpoints[-1]:
            for i in range(1,len(xpoints)):
                if xpoints[i-1]==xpoints[i]: return xvalues[i-1]
                if x>=xpoints[i-1] and x<xpoints[i]:
                    m = (xvalues[i]-xvalues[i-1])/(xpoints[i]-xpoints[i-1])
                    return xvalues[i-1] + m * ( x-xpoints[i-1] )
        # beyond
        else: return 0.
    return f

def cosine(base=1., amplitude=1.,
           xvacuum=0., xlength=None, xnumber=1):
    import math
    global sim_length
    if len(sim_length)>0 and xlength is None: xlength = sim_length[0]-xvacuum
    def f(x,y=0.):
        #vacuum region
        if x < xvacuum: return 0.
        # profile region
        elif x < xvacuum+xlength:
            return base + amplitude * math.cos(2.*math.pi * xnumber * (x-xvacuum)/xlength)
        # beyond
        else: return 0.
    return f

