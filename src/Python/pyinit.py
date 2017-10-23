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
    
    # Function to return one given instance, for example DiagParticleBinning[0]
    # Special case: species can also be indexed by their name: Species["ion1"]
    def __getitem__(self, key):
        if self.__name__ == "Species" and type(key) is str:
            for obj in self._list:
                if obj.name == key:
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
            deprecated = {
                "maxwell_sol":"maxwell_solver",
                "referenceAngularFrequency_SI":"reference_angular_frequency_SI",
                "initPosition_type":"position_initialization",
                "initMomentum_type":"momentum_initialization",
                "thermT":"thermal_boundary_temperature",
                "thermVelocity":"thermal_boundary_velocity",
                "isTest":"is_test",
                "boxSide":"box_side",
                "polarizationPhi":"polarization_phi",
                "dump_file_sequence":"keep_n_dumps",
                "bc_em_type_x":"EM_boundary_conditions",
                "bc_em_type_y":"EM_boundary_conditions",
                "bc_em_type_z":"EM_boundary_conditions",
                "bc_part_type_xmin":"boundary_conditions",
                "bc_part_type_xmax":"boundary_conditions",
                "bc_part_type_ymin":"boundary_conditions",
                "bc_part_type_ymax":"boundary_conditions",
                "bc_part_type_zmin":"boundary_conditions",
                "bc_part_type_zmax":"boundary_conditions",
                "pos"       :"origin",
                "pos_first" :"corners or vectors",
                "pos_second":"corners or vectors",
                "pos_third" :"corners or vectors",
                "track_every"      :"the block DiagTrackParticles()",
                "track_flush_every":"the block DiagTrackParticles()",
                "track_filter"     :"the block DiagTrackParticles()",
                "species_type"     :"name",
                "dynamics_type"    :"pusher",
                "coef_cell"        :"cell_load",
                "coef_frozen"      :"frozen_particle_load",
                "currentFilter_int":"the block CurrentFilter()",
                "Friedman_filter"  :"the block FieldFilter()",
                "Friedman_theta"   :"the block FieldFilter()",
                "n_part_per_cell"  :"particles_per_cell",
                "poisson_iter_max" :"poisson_max_iteration",
                "poisson_error_max":"poisson_max_error",
                "nb_density"       :"number_density",
                "sim_length"       :"grid_length",
                "sim_time"         :"simulation_time",
                "output"           :"deposited_quantity",
            }
            for key, value in kwargs.items():
                if key in deprecated:
                    raise Exception("Deprecated `"+key+"` parameter should be replaced by `"+deprecated[key]+"`")
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

class ParticleData(object):
    """Container for particle data at run-time (for exposing particles in numpy)"""
    _verify = True

class Main(SmileiSingleton):
    """Main parameters"""

    # Default geometry info
    geometry = None
    cell_length = []
    grid_length = []
    number_of_cells = []
    timestep = None
    simulation_time = None
    number_of_timesteps = None
    interpolation_order = 2
    number_of_patches = None
    clrw = -1
    every_clean_particles_overhead = 100
    timestep = None
    timestep_over_CFL = None

    # Poisson tuning
    solve_poisson = True
    poisson_max_iteration = 50000
    poisson_max_error = 1.e-14

    # Default fields
    maxwell_solver = 'Yee'
    EM_boundary_conditions = [["periodic"]]
    time_fields_frozen = 0.
    
    # Default Misc
    reference_angular_frequency_SI = 0.
    print_every = None
    random_seed = None

    def __init__(self, **kwargs):
        # Load all arguments to Main()
        super(Main, self).__init__(**kwargs)
        
        # Deprecation error for the "geometry" argument
        if Main.geometry in ["1d3v", "2d3v", "3d3v"]:
            raise Exception("Deprecated geometry = \""+Main.geometry+"\". Use \""+Main.geometry[0]+"Dcartesian\" instead")
        
        # Initialize timestep if not defined based on timestep_over_CFL
        if Main.timestep is None:
            if Main.timestep_over_CFL is None:
                raise Exception("timestep and timestep_over_CFL not defined")
            else:
                if Main.cell_length is None:
                    raise Exception("Need cell_length to calculate timestep")

                # Yee solver
                if Main.maxwell_solver == 'Yee':
                    dim = int(Main.geometry[0])
                    if dim<1 or dim>3:
                        raise Exception("timestep_over_CFL not implemented in geometry "+Main.geometry)
                    Main.timestep = Main.timestep_over_CFL / math.sqrt(sum([1./l**2 for l in Main.cell_length]))

                # Grassi
                elif Main.maxwell_solver == 'Grassi':
                    if Main.geometry == '2Dcartesian':
                        Main.timestep = Main.timestep_over_CFL * 0.7071067811*Main.cell_length[0];
                    else:
                        raise Exception("timestep_over_CFL not implemented in geometry "+Main.geometry)

                # GrassiSpL
                elif Main.maxwell_solver == 'GrassiSpL':
                    if Main.geometry == '2Dcartesian':
                        Main.timestep = Main.timestep_over_CFL * 0.6471948469*Main.cell_length[0];
                    else:
                        raise Exception("timestep_over_CFL not implemented in geometry "+Main.geometry)

                # None recognized solver
                else:
                    raise Exception("timestep: maxwell_solver not implemented "+Main.maxwell_solver)
                
        #initialize grid_length if not defined based on number_of_cells and cell_length
        if len(Main.grid_length) is 0:
            if len(Main.number_of_cells) is 0:
                raise Exception("grid_length and number_of_cells not defined")
            elif len(Main.number_of_cells) != len(Main.cell_length):
                raise Exception("grid_length and number_of_cells not defined")
            else :
                Main.grid_length = [a*b for a,b in zip(Main.number_of_cells, Main.cell_length)]

        #initialize simulation_time if not defined based on number_of_timesteps and timestep
        if Main.simulation_time is None:
            if Main.number_of_timesteps is None:
                raise Exception("simulation_time and number_of_timesteps not defined")
            else:
                Main.simulation_time = Main.number_of_timesteps * Main.timestep

class LoadBalancing(SmileiSingleton):
    """Load balancing parameters"""

    every = 150
    initial_balance = True
    cell_load = 1.0
    frozen_particle_load = 0.1


class MovingWindow(SmileiSingleton):
    """Moving window parameters"""

    time_start = 0.
    velocity_x = 1.


class Checkpoints(SmileiSingleton):
    """Checkpoints parameters"""
    
    restart_dir = None
    restart_number = None
    dump_step = 0
    dump_minutes = 0.
    keep_n_dumps = 2
    dump_deflate = 0
    exit_after_dump = True
    file_grouping = None
    restart_files = []

class CurrentFilter(SmileiSingleton):
    """Current filtering parameters"""
    model = "binomial"
    passes = 0

class FieldFilter(SmileiSingleton):
    """Fields filtering parameters"""
    model = "Friedman"
    theta = 0.

class Species(SmileiComponent):
    """Species parameters"""
    name = None
    position_initialization = None
    momentum_initialization = ""
    particles_per_cell = None
    c_part_max = 1.0
    mass = None
    charge = None
    charge_density = None
    number_density = None
    mean_velocity = [0.]
    temperature = [1e-10]
    thermal_boundary_temperature = []
    thermal_boundary_velocity = []
    pusher = "boris"
    radiation_model = "none"
    radiation_photon_species = None
    radiation_photon_sampling = 1
    radiation_photon_gamma_threshold = 2
    multiphoton_Breit_Wheeler = [None,None]
    multiphoton_Breit_Wheeler_sampling = [1,1]
    time_frozen = 0.0
    boundary_conditions = [["periodic"]]
    ionization_model = "none"
    ionization_electrons = None
    atomic_number = None
    is_test = False

class Laser(SmileiComponent):
    """Laser parameters"""
    box_side = "xmin"
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
    """Probe diagnostic"""
    every = None
    number = []
    origin = []
    corners = []
    vectors = []
    fields = []
    flush_every = 1

class DiagParticleBinning(SmileiComponent):
    """Particle Binning diagnostic"""
    deposited_quantity = None
    time_average = 1
    species = None
    axes = []
    every = None
    flush_every = 1

class DiagScreen(SmileiComponent):
    """Screen diagnostic"""
    shape = None
    point = None
    vector = None
    direction = "both"
    deposited_quantity = None
    species = None
    axes = []
    time_average = 1
    every = None
    flush_every = 1

class DiagScalar(SmileiComponent):
    """Scalar diagnostic"""
    every = None
    precision = 10
    vars = []

class DiagFields(SmileiComponent):
    """Field diagnostic"""
    every = None
    fields = []
    time_average = 1
    flush_every = 1

class DiagTrackParticles(SmileiComponent):
    """Track diagnostic"""
    species = None
    every = 0
    flush_every = 1
    filter = None

# external fields
class ExternalField(SmileiComponent):
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


# Radiation reaction configuration (continuous and MC algorithms)
class RadiationReaction(SmileiComponent):
    """
    Fine-tuning of synchrotron-like radiation loss
    (classical continuous, quantum correction, stochastics and MC algorithms)
    """
    # Table h parameters
    h_chipa_min = 1e-3
    h_chipa_max = 1e1
    h_dim = 128
    # Table integfochi parameters
    integfochi_chipa_min = 1e-3
    integfochi_chipa_max = 1e1
    integfochi_dim = 128
    # Table xip_chiphmin and xip parameters
    xip_chipa_min = 1e-3
    xip_chipa_max = 1e1
    xip_power = 4
    xip_threshold = 1e-3
    xip_chipa_dim = 128
    xip_chiph_dim = 128
    # Output format, can be "ascii", "binary", "hdf5"
    output_format = "hdf5"
    # Threshold on chipa between the continuous and
    # the discontinuous approaches
    chipa_disc_min_threshold = 1e-2
    # Threshold on chipa: if chipa < 1E-3 no radiation reaction
    chipa_radiation_threshold = 1e-3
    # Path the tables/databases
    table_path = "./"

# MutliphotonBreitWheeler pair creation
class MultiphotonBreitWheeler(SmileiComponent):
    """
    Photon decay into electron-positron pairs
    """
    # Output format, can be "ascii", "binary", "hdf5"
    output_format = "hdf5"
    # Path the tables/databases
    table_path = "./"
    # Table T parameters
    T_chiph_min = 1e-2
    T_chiph_max = 1e1
    T_dim = 128
    # Table xip parameters
    xip_chiph_min = 1e-2
    xip_chiph_max = 1e1
    xip_power = 4
    xip_threshold = 1e-3
    xip_chipa_dim = 128
    xip_chiph_dim = 128

# Smilei-defined
smilei_mpi_rank = 0
smilei_mpi_size = 1
smilei_rand_max = 2**31-1

# DEPRECATION ERRORS
class DiagParticles(object):
    def __init__(self, *args, **kwargs):
        raise Exception("Deprecated `DiagParticles()` must be replaced by `DiagParticleBinning()`")
class DumpRestart(object):
    def __init__(self, *args, **kwargs):
        raise Exception("Deprecated `DumpRestart()` must be replaced by `Checkpoints()`")
class ExtField(object):
    def __init__(self, *args, **kwargs):
        raise Exception("Deprecated `ExtField()` must be replaced by `ExternalField()`")

