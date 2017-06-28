"""
    GENERAL DEFINITIONS FOR SMILEI
"""

import math, os, gc, operator

def _add_metaclass(metaclass):
    """Class decorator for creating a class with a metaclass."""
    # Taken from module "six" for compatibility with python 2 and 3
    def wrapper(cls):
        orig_vars = cls.__dict__.copy()
        slots = orig_vars.get('__slots__')
        if slots is not None:
            if isinstance(slots, str): slots = [slots]
            for slots_var in slots: orig_vars.pop(slots_var)
        orig_vars.pop('__dict__', None)
        orig_vars.pop('__weakref__', None)
        return metaclass(cls.__name__, cls.__bases__, orig_vars)
    return wrapper


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
    __next__ = next #python3

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

    # Function to display the list of instances
    def __repr__(self):
        if len(self._list)==0:
            return "<Empty list of "+self.__name__+">"
        else:
            l = []
            for obj in self._list: l.append(str(obj))
            return "["+", ".join(l)+"]"

@_add_metaclass(SmileiComponentType)
class SmileiComponent(object):
    """Smilei component generic class"""

    # Function to initialize new components
    def _fillObjectAndAddToList(self, cls, obj, **kwargs):
        # add all kwargs as internal class variables
        if kwargs is not None:
            for key, value in kwargs.items():
                if key=="_list":
                    print("Python warning: in "+cls.__name__+": cannot have argument named '_list'. Discarding.")
                elif not hasattr(cls, key):
                    raise Exception("ERROR in the namelist: cannot define `"+key+"` in block "+cls.__name__+"()");
                else:
                    setattr(obj, key, value)
        # add the new component to the "_list"
        cls._list.append(obj)

    # Constructor for all SmileiComponents
    def __init__(self, **kwargs):
        self._fillObjectAndAddToList(type(self), self, **kwargs)

    def __repr__(self):
        return "<Smilei "+type(self).__name__+">"


class SmileiSingletonType(SmileiComponentType):
    """Metaclass to all Smilei singletons"""

    def __repr__(self):
        return "<Smilei "+str(self.__name__)+">"

@_add_metaclass(SmileiSingletonType)
class SmileiSingleton(SmileiComponent):
    """Smilei singleton generic class"""

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
    number_of_cells = []
    timestep = None
    sim_time = None
    number_of_timesteps = None
    interpolation_order = 2
    number_of_patches = None
    clrw = 1
    every_clean_particles_overhead = 100
    timestep = None
    timestep_over_CFL = None

    # Poisson tuning
    solve_poisson = True
    poisson_iter_max = 50000
    poisson_error_max = 1.e-14

    # Default fields
    maxwell_sol = 'Yee'
    bc_em_type_x = []
    bc_em_type_y = []
    bc_em_type_z = []
    time_fields_frozen = 0.
    currentFilter_int = 0
    Friedman_filter = False
    Friedman_theta = 0.

    # Default Misc
    referenceAngularFrequency_SI = 0.
    print_every = None
    random_seed = None

    def __init__(self, **kwargs):
        super(Main, self).__init__(**kwargs)
        #initialize timestep if not defined based on timestep_over_CFL
        if Main.timestep is None:
            if Main.timestep_over_CFL is None:
                raise Exception("timestep and timestep_over_CFL not defined")
            else:
                if Main.cell_length is None:
                    raise Exception("Need cell_length to calculate timestep")

                # Yee solver
                if Main.maxwell_sol == 'Yee':
                    dim = int(Main.geometry[0])
                    if dim<1 or dim>3:
                        raise Exception("timestep_over_CFL not implemented in geometry "+Main.geometry)
                    Main.timestep = Main.timestep_over_CFL / math.sqrt(sum([1./l**2 for l in Main.cell_length]))

                # Grassi
                elif Main.maxwell_sol == 'Grassi':
                    if Main.geometry == '2d3v':
                        Main.timestep = Main.timestep_over_CFL * 0.7071067811*Main.cell_length[0];
                    else:
                        raise Exception("timestep_over_CFL not implemented in geometry "+Main.geometry)

                # GrassiSpL
                elif Main.maxwell_sol == 'GrassiSpL':
                    if Main.geometry == '2d3v':
                        Main.timestep = Main.timestep_over_CFL * 0.6471948469*Main.cell_length[0];
                    else:
                        raise Exception("timestep_over_CFL not implemented in geometry "+Main.geometry)

                # None recognized solver
                else:
                    raise Exception("timestep: maxwell_sol not implemented "+Main.maxwell_sol)

        #initialize sim_length if not defined based on number_of_cells and cell_length
        if len(Main.sim_length) is 0:
            if len(Main.number_of_cells) is 0:
                raise Exception("sim_length and number_of_cells not defined")
            elif len(Main.number_of_cells) != len(Main.cell_length):
                raise Exception("sim_length and number_of_cells not defined")
            else :
                Main.sim_length = [a*b for a,b in zip(Main.number_of_cells, Main.cell_length)]

        #initialize sim_time if not defined based on number_of_timesteps and timestep
        if Main.sim_time is None:
            if Main.number_of_timesteps is None:
                raise Exception("sim_time and number_of_timesteps not defined")
            else:
                Main.sim_time = Main.number_of_timesteps * Main.timestep

class LoadBalancing(SmileiSingleton):
    """Load balancing parameters"""

    every = 150
    initial_balance = True
    coef_cell = 1.0
    coef_frozen = 0.1


class MovingWindow(SmileiSingleton):
    """Moving window parameters"""

    time_start = 0.
    velocity_x = 1.


class DumpRestart(SmileiSingleton):
    """Dump and restart parameters"""

    restart_dir = None
    restart_number = None
    dump_step = 0
    dump_minutes = 0.
    dump_file_sequence = 2
    dump_deflate = 0
    exit_after_dump = True
    file_grouping = None


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
    thermT = None
    thermVelocity = None
    dynamics_type = "norm"
    radiation_type = "none"
    time_frozen = 0.0
    radiating = False
    bc_part_type_xmin = None
    bc_part_type_xmax = None
    bc_part_type_ymax = None
    bc_part_type_ymin = None
    bc_part_type_zmin = None
    bc_part_type_zmax = None
    ionization_model = "none"
    ionization_electrons = None
    atomic_number = None
    isTest = False
    track_every = 0
    track_flush_every = 1
    track_filter = None

class Laser(SmileiComponent):
    """Laser parameters"""
    boxSide = "xmin"
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
    flush_every = 1

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    output = None
    time_average = 1
    species = None
    axes = []
    every = None
    flush_every = 1

class DiagScreen(SmileiComponent):
    """Diagnostic particles"""
    shape = None
    point = None
    vector = None
    direction = "both"
    output = None
    species = None
    axes = []
    time_average = 1
    every = None
    flush_every = 1

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
    flush_every = 1

# external fields
class ExtField(SmileiComponent):
    """External Field"""
    field = None
    profile = None

# external current (antenna)
class Antenna(SmileiComponent):
    """Antenna"""
    field = None
    time_profile  = None
    space_profile = None

# Particle wall
class PartWall(SmileiComponent):
    """Particle Wall"""
    kind = "none"
    x = None
    y = None
    z = None

# Radiation loss (continuous and MC algorithms)
class RadiationLoss(SmileiComponent):
    """
    Synchrotron-like radiation loss 
    (classical continuous, quantum correction, stochastics and MC algorithms)
    """
    # First table parameters
    chipa_integfochi_min = 1e-5
    chipa_integfochi_max = 1e2
    integfochi_dim = 128
    # Second table parameters
    chipa_xip_min = 1e-5
    chipa_xip_max = 1e2
    xip_power = 4
    xip_threshold = 1e-3
    chipa_xip_dim = 128
    chiph_xip_dim = 128
    # Output format, can be "ascii", "binary", "hdf5"
    output_format = "hdf5"
    # Threshold on chipa between the continuous and
    # the discontinuous approaches
    chipa_disc_min_threshold = 1e-2
    # Path the tables/databases
    table_path = "./"

# Smilei-defined
smilei_mpi_rank = 0
smilei_mpi_size = 1
smilei_rand_max = 2**31-1
