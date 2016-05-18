"""@package pyinit
    Definition of Smilei components
"""



class SmileiComponentType(type):
    """Metaclass to all Smilei components"""
    
    # Constructor of classes
    def __init__(self, name, bases, attrs):
        self._list = []
        self._verify = True
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
    
    # Function to initialize new components
    def _fillObjectAndAddToList(self, cls, obj, **kwargs):
        # add all kwargs as internal class variables
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                if key=="_list":
                    print "Python warning: in "+cls.__name__+": cannot have argument named '_list'. Discarding."
                elif not hasattr(cls, key):
                    raise Exception("ERROR in the namelist: cannot define `"+key+"` in block "+cls.__name__+"()");
                else:
                    setattr(obj, key, value)
        # add the new component to the "_list"
        cls._list.append(obj)
    
    # Constructor for all SmileiComponents
    def __init__(self, **kwargs):
        self._fillObjectAndAddToList(type(self), self, **kwargs)


class SmileiSingletonType(SmileiComponentType):
    """Metaclass to all Smilei singletons"""
    
    def __repr__(self):
        return "<Smilei "+str(self.__name__)+">"

class SmileiSingleton(SmileiComponent):
    """Smilei singleton generic class"""
    __metaclass__ = SmileiSingletonType
    
    # Prevent having two instances
    def __new__(cls, **kwargs):
        if len(cls._list) >= 1:
            raise Exception("ERROR in the namelist: cannot define block "+cls.__name__+"() twice")
        return super(SmileiSingleton, cls).__new__(cls)
    
    # Constructor for all SmileiSingletons
    def __init__(self, **kwargs):
        self._fillObjectAndAddToList(type(self), type(self), **kwargs)


class Main(SmileiSingleton):
    """Main parameters"""
    
    # Default geometry info
    geometry = None
    cell_length = []
    sim_length = []
    timestep = None
    sim_time = None
    interpolation_order = 2
    number_of_patches = None
    clrw = 1
    
    # Default fields
    maxwell_sol = 'Yee'
    bc_em_type_x = []
    bc_em_type_y = []
    time_fields_frozen = 0.
    
    # Default Misc
    referenceAngularFrequency_SI = 0.
    print_every = None
    output_dir = None
    random_seed = None


class LoadBalancing(SmileiSingleton):
    """Load balancing parameters"""
    
    every = None
    coef_cell = 1.0
    coef_frozen = 0.1


class MovingWindow(SmileiSingleton):
    """Moving window parameters"""
    
    delay = 0.
    velocity_x = 1.
    nspace_x = 0


class DumpRestart(SmileiSingleton):
    """Dump and restart parameters"""
    
    restart_dir = None
    dump_step = 0
    dump_minutes = 0.
    dump_file_sequence = 2
    dump_deflate = 0
    exit_after_dump = True


class Species(SmileiComponent):
    """Species parameters"""
    species_type = None
    initPosition_type = None
    initMomentum_type = ""
    n_part_per_cell = None
    c_part_max = 1.0
    mass = None
    charge = None
    charge_density = None
    nb_density = None
    mean_velocity = [0.]
    temperature = [1e-10]
    dynamics_type = "norm"
    time_frozen = 0.0
    radiating = False
    bc_part_type_west = None
    bc_part_type_east = None
    bc_part_type_north = None
    bc_part_type_south = None
    ionization_model = "none"
    ionization_electrons = None
    atomic_number = None
    isTest = False
    track_every = 0

class Laser(SmileiComponent):
    """Laser parameters"""
    boxSide = "west"
    omega = 1.
    chirp_profile = 1.
    time_envelope = 1.
    space_envelope = [1., 0.]
    phase = [0., 0.]
    space_time_profile = None

class Collisions(SmileiComponent):
    """Collisions parameters"""
    species1 = None
    species2 = None
    coulomb_log = 0.
    debug_every = 0
    ionizing = False


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    every = None
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

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    every = None
    precision = 10
    vars = []

class DiagFields(SmileiComponent):
    """Diagnostic Fields"""
    every = None
    fields = []
    time_average = 1

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

# Smilei-defined
smilei_mpi_rank = 0
smilei_mpi_size = 1
smilei_rand_max = 2**31-1
