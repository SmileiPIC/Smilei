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
    # Special case: diagnostics can also be indexed by their label: DiagParticleBinning["x-px"]
    def __getitem__(self, key):
        try:
            for obj in self._list:
                if obj.name == key:
                    return obj
        except:
            pass
        return self._list[key]
    
    def has(self, key):
        if type(key) is int and key < len(self._list):
            return True
        for obj in self._list:
            if obj.name == key:
                return True
        return False
    
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
                "output_format":"See documentation for radiation reaction",
                "clrw":"See https://smileipic.github.io/Smilei/namelist.html#main-variables",
                "h_chipa_min":"See documentation for radiation reaction",
                "h_chipa_max":"See documentation for radiation reaction",
                "h_dim":"See documentation for radiation reaction",
                "h_computation_method":"See documentation for radiation reaction",
                "integfochi_chipa_min":"See documentation for radiation reaction",
                "integfochi_chipa_max":"See documentation for radiation reaction",
                "integfochi_dim":"See documentation for radiation reaction",
                "xip_chipa_min":"See documentation for radiation reaction",
                "xip_chipa_max":"See documentation for radiation reaction",
                "xip_power":"See documentation for radiation reaction or Breit-Wheeler",
                "xip_threshold":"See documentation for radiation reaction or Breit-Wheeler",
                "xip_chipa_dim":"See documentation for radiation reaction or Breit-Wheeler",
                "xip_chiph_dim":"See documentation for radiation reaction or Breit-Wheeler",
                "compute_table":"See documentation for Breit-Wheeler",
                "T_chiph_min":"See documentation for Breit-Wheeler",
                "T_chiph_max":"See documentation for Breit-Wheeler",
                "T_dim":"See documentation for Breit-Wheeler",
                "xip_chiph_min":"See documentation for Breit-Wheeler",
                "xip_chiph_max":"See documentation for Breit-Wheeler",
            }
            for key, value in kwargs.items():
                if key in deprecated:
                    raise Exception("Deprecated `"+key+"` parameter. "+deprecated[key])
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
    interpolator = "momentum-conserving"
    custom_oversize = 2
    number_of_patches = None
    patch_arrangement = "hilbertian"
    cluster_width = -1
    every_clean_particles_overhead = 100
    timestep = None
    number_of_AM = 2
    number_of_AM_relativistic_field_initialization = 1
    number_of_AM_classical_Poisson_solver = 1
    timestep_over_CFL = None
    cell_sorting = None
    gpu_computing = False                      # Activate the computation on GPU
    
    # PXR tuning
    spectral_solver_order = []
    initial_rotational_cleaning = False

    # Poisson tuning
    solve_poisson = True
    poisson_max_iteration = 50000
    poisson_max_error = 1.e-14

    # Relativistic Poisson tuning
    solve_relativistic_poisson = False
    relativistic_poisson_max_iteration = 50000
    relativistic_poisson_max_error = 1.e-22
    
    # BTIS3 interpolator
    use_BTIS3_interpolation = False

    # Default fields
    maxwell_solver = 'Yee'
    EM_boundary_conditions = [["periodic"]]
    EM_boundary_conditions_k = []
    save_magnectic_fields_for_SM = True
    number_of_pml_cells = [[10]]

    def default_sigma(x):
        return 20. * x**2  
    def default_integrate_sigma(x):
        return 20./3. * x**3  
    def default_kappa(x):
        return 1 + 79. * x**4  
    def default_integrate_kappa(x):
        return x + 79./5. * x**5  

    pml_sigma = [default_sigma, default_sigma, default_sigma, default_integrate_sigma]
    pml_kappa = [default_kappa, default_kappa, default_kappa, default_integrate_kappa]
    time_fields_frozen = 0.
    Laser_Envelope_model = False

    # Default Misc
    reference_angular_frequency_SI = 0.
    print_every = None
    random_seed = None
    print_expected_disk_usage = True

    terminal_mode = True

    def __init__(self, **kwargs):
        # Load all arguments to Main()
        super(Main, self).__init__(**kwargs)

        # Initialize timestep if not defined based on timestep_over_CFL
        if Main.timestep is None:
            if Main.timestep_over_CFL is None:
                raise Exception("ERROR in the namelist in Main: timestep and timestep_over_CFL not defined")
            else:
                if Main.cell_length is None:
                    raise Exception("ERROR in the namelist in Main: cell_length must be defined in order to calculate timestep from timestep_over_CFL")

                # Yee solver or Terzani solver
                if Main.maxwell_solver in ['Yee','Terzani']:
                    if (Main.geometry=="AMcylindrical"):
                        alpha = [0.210486, 0.591305, 3.5234, 8.51041, 15.5059]
                        if (Main.number_of_AM < 6):
                            alpha = 1. + alpha[Main.number_of_AM-1]
                        else:
                            alpha = (Main.number_of_AM-1)**2
                        Main.timestep = Main.timestep_over_CFL / math.sqrt(1./Main.cell_length[0]**2 + alpha/Main.cell_length[1]**2 )
                        if (Main.maxwell_solver == "Terzani"):
                            print("WARNING CFL: The timestep has been set as for a Yee solver but the maximum stable timestep for the Terzani solver may be slightly smaller.")
                    else:
                        dim = int(Main.geometry[0])
                        Main.timestep = Main.timestep_over_CFL / math.sqrt(sum([1./l**2 for l in Main.cell_length]))

                # Grassi
                elif Main.maxwell_solver == 'Grassi':
                    if Main.geometry == '2Dcartesian':
                        Main.timestep = Main.timestep_over_CFL * 0.7071067811*Main.cell_length[0];
                    else:
                        raise Exception("ERROR in the namelist in Main: timestep_over_CFL for the Grassi solver not implemented in geometry "+Main.geometry)

                # GrassiSpL
                elif Main.maxwell_solver == 'GrassiSpL':
                    if Main.geometry == '2Dcartesian':
                        Main.timestep = Main.timestep_over_CFL * 0.6471948469*Main.cell_length[0];
                    else:
                        raise Exception("ERROR in the namelist in Main: timestep_over_CFL for the GrassiSpL solver not implemented in geometry "+Main.geometry)

                # M4
                elif Main.maxwell_solver == 'M4':
                    if Main.geometry == '1Dcartesian':
                        Main.timestep = Main.timestep_over_CFL * Main.cell_length[0];
                    elif Main.geometry == '2Dcartesian':
                        Main.timestep = Main.timestep_over_CFL * min(Main.cell_length[0:2])
                    elif Main.geometry == '3Dcartesian':
                        Main.timestep = Main.timestep_over_CFL * min(Main.cell_length)
                    else:
                        raise Exception("ERROR in the namelist in Main: timestep_over_CFL for the M4 solver not implemented in geometry "+Main.geometry)

                # None recognized solver
                else:
                    raise Exception("ERROR in the namelist in Main: timestep_over_CFL not implemented for the solver "+Main.maxwell_solver)

        # Constraint on timestep for WT interpolation
        if Main.interpolator.lower() == "wt":
            if Main.geometry == '1Dcartesian':
                if Main.timestep > 0.5 * Main.cell_length[0]:
                    raise Exception("timestep for WT cannot be larger than 0.5*dx")
            elif Main.geometry == '2Dcartesian':
                if Main.timestep > 0.5 * min(Main.cell_length[0:2]):
                    raise Exception("timestep for WT cannot be larger than 0.5*min(dx,dy)")
            elif Main.geometry == '3Dcartesian':
                if Main.timestep > 0.5 * min(Main.cell_length):
                    raise Exception("timestep for WT cannot be larger than 0.5*min(dx,dy,dz)")
            else:
                raise Exception("WT interpolation not implemented in geometry "+Main.geometry)

        # Initialize simulation_time if not defined by the user
        if Main.simulation_time is None:
            if Main.number_of_timesteps is None:
                raise Exception("ERROR in the namelist in Main: simulation_time and number_of_timesteps are not defined")
            Main.simulation_time = Main.timestep * Main.number_of_timesteps

        # Initialize grid_length if not defined based on number_of_cells and cell_length
        if (    len(Main.grid_length + Main.number_of_cells) == 0
             or len(Main.grid_length + Main.cell_length) == 0
             or len(Main.number_of_cells + Main.cell_length) == 0
             or len(Main.number_of_cells) * len(Main.grid_length) * len(Main.cell_length) != 0
           ):
                raise Exception("ERROR in the namelist in Main: you must define exactly two parameters among grid_length, number_of_cells and cell_length")

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

class MultipleDecomposition(SmileiSingleton):
    """Multiple Decomposition parameters"""
    region_ghost_cells   = 2

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
    number_of_additional_shifts = 0
    additional_shifts_time = 0.


class Checkpoints(SmileiSingleton):
    """Checkpoints parameters"""

    restart_dir = None
    restart_number = None
    dump_step = 0
    dump_minutes = 0.
    keep_n_dumps = 2
    dump_deflate = 0
    exit_after_dump = True
    file_grouping = 0
    restart_files = []

class CurrentFilter(SmileiSingleton):
    """Current filtering parameters"""
    model = "binomial"
    passes = [0]
    kernelFIR = [0.25,0.5,0.25]

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
    regular_number = []
    c_part_max = 1.0
    mass = None
    charge = None
    charge_density = None
    number_density = None
    mean_velocity = []  # Default value is     0, set in ParticleCreator function in species.cpp
    mean_velocity_AM = []
    temperature = []    # Default value is 1e-10, set in ParticleCreator function in species.cpp
    thermal_boundary_temperature = []
    thermal_boundary_velocity = [0.,0.,0.]
    pusher = "boris"

    # Radiation species parameters
    radiation_model = "none"
    radiation_photon_species = None
    radiation_photon_sampling = 1
    radiation_photon_gamma_threshold = 2.
    radiation_max_emissions = 10

    # Multiphoton Breit-Wheeler parameters
    multiphoton_Breit_Wheeler = [None,None]
    multiphoton_Breit_Wheeler_sampling = [1,1]

    # Particle merging species Parameters
    merging_method = "none"
    merge_every = 0
    merge_min_packet_size = 4
    merge_max_packet_size = 4
    merge_min_particles_per_cell = 4
    merge_min_momentum_cell_length = [1e-10,1e-10,1e-10]
    merge_momentum_cell_size = [16,16,16]
    merge_accumulation_correction = True
    merge_discretization_scale = "linear"
    merge_min_momentum = 1e-5

    time_frozen = 0.0
    radiating = False
    relativistic_field_initialization = False
    boundary_conditions = [["periodic"]]
    ionization_model = "none"
    ionization_electrons = None
    ionization_rate = None
    atomic_number = None
    maximum_charge_state = 0
    ionization_tl_parameter = 6
    is_test = False
    relativistic_field_initialization = False
    keep_interpolated_fields = []

class ParticleInjector(SmileiComponent):
    """Parameters for particle injection at boundaries"""
    name = None
    species = None
    box_side = "xmin"
    position_initialization = "species"
    momentum_initialization = "species"
    mean_velocity = []  # Default value is     0, set in ParticleCreator function
    temperature = []    # Default value is 1e-10, set in ParticleCreator function
    charge_density = None
    number_density = None
    particles_per_cell = None
    regular_number = []
    time_envelope = 1
    
    
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
    space_time_profile_AM = None
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
    '''
    Defaut parameters tested in AM geometry for
    dx = 0.5 ; dy = 2*dx ; dt = 0.8*dx with 'envelope_reduced_dispersion' solver
    PML cells = 16
    Give good results for a propagation during 40 000 dt
    '''
    Env_pml_sigma_parameters = [[0.00,2],[20.0,2],[20.0,2]]
    Env_pml_kappa_parameters = [[1.00,1.00,2],[1.00,1.00,2],[1.00,1.00,2]]
    Env_pml_alpha_parameters = [[0.00,0.00,1],[0.05,0.05,1],[0.05,0.05,1]]
    polarization_phi = 0.
    ellipticity = 0.


class Collisions(SmileiComponent):
    """Collisions parameters"""
    species1 = None
    species2 = None
    coulomb_log = 0.
    coulomb_log_factor = 1.
    every = 1
    debug_every = 0
    time_frozen = 0
    ionizing = False
    nuclear_reaction = None
    nuclear_reaction_multiplier = 0.


#diagnostics
class DiagProbe(SmileiComponent):
    """Probe diagnostic"""
    name = ""
    every = None
    number = []
    origin = []
    corners = []
    vectors = []
    fields = []
    flush_every = 1
    time_integral = False
    datatype = "double"

class DiagParticleBinning(SmileiComponent):
    """Particle Binning diagnostic"""
    name = ""
    deposited_quantity = None
    time_average = 1
    species = None
    axes = []
    every = None
    flush_every = 1

class DiagRadiationSpectrum(SmileiComponent):
    """Radiation Spectrum diagnostic"""
    name = ""
    time_average = 1
    species = None
    photon_energy_axis = None
    axes = []
    every = None
    flush_every = 1

class DiagScreen(SmileiComponent):
    """Screen diagnostic"""
    name = ""
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
    name = ""
    every = None
    fields = []
    time_average = 1
    subgrid = None
    flush_every = 1
    datatype = "double"

class DiagTrackParticles(SmileiComponent):
    """Track diagnostic"""
    name = ""
    species = None
    every = 0
    flush_every = 1
    filter = None
    attributes = ["x", "y", "z", "px", "py", "pz", "w"]

class DiagNewParticles(SmileiComponent):
    """Track diagnostic"""
    name = ""
    species = None
    every = 0
    flush_every = 1
    filter = None
    attributes = ["x", "y", "z", "px", "py", "pz", "w"]

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

# external time fields
class PrescribedField(SmileiComponent):
    """External Time Field"""
    field = None
    profile = None

# external current (antenna)
class Antenna(SmileiComponent):
    """Antenna"""
    field = None
    time_profile  = None
    space_profile = None
    space_time_profile = None

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
    # Minimum particle_chi value for the discontinuous radiation
    # Under this value, the discontinuous approach is not applied
    minimum_chi_discontinuous = 1e-2
    # Threshold on particle_chi: if particle_chi < 1E-3 no radiation reaction
    minimum_chi_continuous = 1e-3

    # Path to read or write the tables/databases
    table_path = ""

    # Parameters for computing the tables
    Niel_computation_method = "table"

# MultiphotonBreitWheeler pair creation
class MultiphotonBreitWheeler(SmileiComponent):
    """
    Photon decay into electron-positron pairs
    """
    # Path the tables/databases
    table_path = ""

# Smilei-defined
smilei_mpi_rank = 0
smilei_mpi_size = 1
smilei_omp_threads = 1
smilei_total_cores = 1

# Variable to set to False for the actual run (useful for the test mode)
_test_mode = True
