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
                "chipa_radiation_threshold":"minimum_chi_continuous",
                "chipa_disc_min_threshold":"minimum_chi_discontinuous",
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
        cls = type(self)
        self._fillObjectAndAddToList(cls, cls, **kwargs)
        # Change all methods to static
        for k,v in cls.__dict__.items():
            if k[0]!='_' and hasattr(v,"__get__"):
               setattr(cls, k, staticmethod(v))

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
    patch_arrangement = "hilbertian"
    clrw = -1
    every_clean_particles_overhead = 100
    timestep = None
    number_of_AM = 2
    timestep_over_CFL = None

    # PXR tuning
    global_factor = []
    norder = []
    is_spectral = False
    is_pxr = False

    # Poisson tuning
    solve_poisson = True
    poisson_max_iteration = 50000
    poisson_max_error = 1.e-14

    # Relativistic Poisson tuning
    solve_relativistic_poisson = False
    relativistic_poisson_max_iteration = 50000
    relativistic_poisson_max_error = 1.e-22

    # Default fields
    maxwell_solver = 'Yee'
    EM_boundary_conditions = [["periodic"]]
    EM_boundary_conditions_k = []
    save_magnectic_fields_for_SM = True
    time_fields_frozen = 0.
    Laser_Envelope_model = False

    # Default Misc
    reference_angular_frequency_SI = 0.
    print_every = None
    random_seed = None
    print_expected_disk_usage = True

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

        # Initialize simulation_time if not defined by the user
        if Main.simulation_time is None:
            if Main.number_of_timesteps is None:
                raise Exception("simulation_time and number_of_timesteps are not defined")
            Main.simulation_time = Main.timestep * Main.number_of_timesteps

        # Initialize grid_length if not defined based on number_of_cells and cell_length
        if (    len(Main.grid_length + Main.number_of_cells) == 0
             or len(Main.grid_length + Main.cell_length) == 0
             or len(Main.number_of_cells + Main.cell_length) == 0
             or len(Main.number_of_cells) * len(Main.grid_length) * len(Main.cell_length) != 0
           ):
                raise Exception("Main: you must define two (and only two) between grid_length, number_of_cells and cell_length")

        if len(Main.grid_length) == 0:
            Main.grid_length = [a*b for a,b in zip(Main.number_of_cells, Main.cell_length)]

        if len(Main.cell_length) == 0:
            Main.cell_length = [a/b for a,b in zip(Main.grid_length, Main.number_of_cells)]

        if len(Main.number_of_cells) == 0:
            Main.number_of_cells = [int(round(float(a)/float(b))) for a,b in zip(Main.grid_length, Main.cell_length)]
            old_grid_length = Main.grid_length
            Main.grid_length = [a*b for a,b in zip(Main.number_of_cells, Main.cell_length)]
            if smilei_mpi_rank == 0:
                different = [abs((a-b)/(a+b))>1e-10 for a,b in zip(Main.grid_length, old_grid_length)]
                if any(different):
                    print("\t[Python WARNING] Main.grid_length="+str(Main.grid_length)+" (was "+str(old_grid_length)+")")

class LoadBalancing(SmileiSingleton):
    """Load balancing parameters"""

    every                = 150
    initial_balance      = True
    cell_load            = 1.0
    frozen_particle_load = 0.1

# Radiation reaction configuration (continuous and MC algorithms)
class Vectorization(SmileiSingleton):
    """
    Vectorization parameters
    """
    mode                = "off"
    reconfigure_every   = 20
    initial_mode        = "off"


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
    mean_velocity = []  # Default value is     0, set in createParticles function in species.cpp
    temperature = []    # Default value is 1e-10, set in createParticles function in species.cpp
    thermal_boundary_temperature = []
    thermal_boundary_velocity = [0.,0.,0.]
    pusher = "boris"

    # Radiation species parameters
    radiation_model = "none"
    radiation_photon_species = None
    radiation_photon_sampling = 1
    radiation_photon_gamma_threshold = 2

    # Multiphoton Breit-Wheeler parameters
    multiphoton_Breit_Wheeler = [None,None]
    multiphoton_Breit_Wheeler_sampling = [1,1]

    # Particle merging species Parameters
    merging_method = "none"
    merge_every = 0

    time_frozen = 0.0
    radiating = False
    relativistic_field_initialization = False
    time_relativistic_initialization = 0.0
    boundary_conditions = [["periodic"]]
    ionization_model = "none"
    ionization_electrons = None
    ionization_rate = None
    atomic_number = None
    maximum_charge_state = None
    is_test = False
    relativistic_field_initialization = False
    ponderomotive_dynamics = False

class Laser(SmileiComponent):
    """Laser parameters"""
    box_side = "xmin"
    omega = 1.
    chirp_profile = 1.
    time_envelope = 1.
    space_envelope = [1., 0.]
    phase = [0., 0.]
    delay_phase = [0., 0.]
    space_time_profile = None
    file = None
    _offset = None

class LaserEnvelope(SmileiSingleton):
    """Laser Envelope parameters"""
    omega = 1.
    #time_envelope = 1.
    #space_envelope = [1., 0.]
    envelope_solver = "explicit"
    envelope_profile = None
    Envelope_boundary_conditions = [["reflective"]]


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
    subgrid = None
    flush_every = 1

class DiagTrackParticles(SmileiComponent):
    """Track diagnostic"""
    species = None
    every = 0
    flush_every = 1
    filter = None
    attributes = ["x", "y", "z", "px", "py", "pz"]

class DiagPerformances(SmileiSingleton):
    """Performances diagnostic"""
    every = 0
    flush_every = 1
    patch_information = True

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
    Fine-tuning of synchrotron-like radiation reaction
    (classical continuous, quantum correction, stochastics and MC algorithms)
    """
    # Table h parameters
    h_chipa_min = 1e-3
    h_chipa_max = 1e1
    h_dim = 128
    h_computation_method = "table"
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
    # Minimum particle_chi value for the discontinuous radiation
    # Under this value, the discontinuous approach is not applied
    minimum_chi_discontinuous = 1e-2
    # Threshold on particle_chi: if particle_chi < 1E-3 no radiation reaction
    minimum_chi_continuous = 1e-3
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

# Variable to set to False for the actual run (useful for the test mode)
_test_mode = True
