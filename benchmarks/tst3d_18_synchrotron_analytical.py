# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Synchrotron: radiation loss of electrons rotating in a constant
#              magnetic field using the MC model
#
# In this tests case, the MC model is compared to analytics after 1 iteration.
#
# Validation:
# - Monte-Carlo radiation loss
# - Species scalar diagnostics
# - External fields
# - Particle binning with the quantum parameter
# ----------------------------------------------------------------------------------------

import math
import datetime

# ----------------------------------------------------------------------------------------
# Main parameters

alpha = 1./137.
HBar = 1.0545718e-34                              # (*m^2kg/s*)
c = 299792458
electron_mass = 9.10938356e-31                    # (*kg*)
electron_charge = 1.60217662e-19                  # (*coulombs*)
lambdar = 1e-6                                    # Normalization wavelength
wr = 2*math.pi*c/lambdar                          # Normalization frequency

Schwinger_E_field= 1.3E18                         # Schwinger electric field
Enorm = electron_mass*wr*c/electron_charge        # Normalization electric field at lambda = 1e-6

l0 = 2.0*math.pi                                  # laser wavelength

chi = 1.                                          # Initial quantum parameter
B = 270.                                          # Magnetic field strength
gamma = chi* Schwinger_E_field/(Enorm*B)          # Initial gamma factor
Rsync = math.sqrt(gamma**2 - 1.)/B                # Synchrotron radius without radiation
v = math.sqrt(1.-1./gamma**2)                     # Initial particle velocity

field_vector = np.array([1.6, 0.7, 0.8])
field_vector_norm = np.sqrt(np.sum(np.power(field_vector,2)))
field_vector /= field_vector_norm

other_vector = np.array([0.8, 1.7, 1.4])

orthogonal_vector = np.array([0., 0., 0.])
orthogonal_vector[0] = field_vector[1]*other_vector[2] - field_vector[2]*other_vector[1]
orthogonal_vector[1] = field_vector[2]*other_vector[0] - field_vector[0]*other_vector[2]
orthogonal_vector[2] = field_vector[0]*other_vector[1] - field_vector[1]*other_vector[0]
orthogonal_vector_norm = np.sqrt(np.sum(np.power(orthogonal_vector,2)))
orthogonal_vector /= orthogonal_vector_norm

mean_velocity = v * orthogonal_vector

Lx = 4.*Rsync
Ly = 4.*Rsync
Lz = 4.*Rsync

n0 = 1e-5                                         # particle density

res = 32.                                        # nb of cells in one synchrotron radius

dx = Rsync/res                                    # space step
dy = Rsync/res                                    # space step
dz = Rsync/res                                    # space step
dt  = 1./math.sqrt(1./(dx*dx) + 1./(dy*dy) + 1./(dz*dz))       # timestep given by the CFL

dt_factor = 0.9                                   # factor on dt
dt *= dt_factor                                   # timestep used for the simulation
Tsim = 1*dt                          # duration of the simulation

pusher = "vay"                                    # type of pusher
radiation_list = ["Monte-Carlo",]    # List of radiation models for species
species_name_list = ["MC"]               # List of names for species

datetime = datetime.datetime.now()
random_seed = datetime.microsecond

# ----------------------------------------------------------------------------------------
# Functions

# Density profile for inital location of the particles
def n0_(x,y,z):
    return n0

# ----------------------------------------------------------------------------------------
# Namelists

Main(
    geometry = "3Dcartesian",

    interpolation_order = 4,

    cell_length = [dx,dy, dz],
    grid_length  = [Lx,Ly, Lz],

    number_of_patches = [8,8,8],

    time_fields_frozen = Tsim,

    timestep = dt,
    simulation_time = Tsim,

    EM_boundary_conditions = [['periodic'],['periodic'],['periodic']],

    random_seed = random_seed,
    
    solve_poisson = False,

    reference_angular_frequency_SI = wr

)

# ----------------------------------------------------------------------------------------
# Initialization of the constant external field

ExternalField(
    field = "Bx",
    profile = constant(B*field_vector[0])
)

ExternalField(
    field = "By",
    profile = constant(B*field_vector[1])
)

ExternalField(
    field = "Bz",
    profile = constant(B*field_vector[2])
)

# ----------------------------------------------------------------------------------------

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0,
    mean_velocity = mean_velocity,
    temperature = [0.],
    pusher = pusher,
    boundary_conditions = [["periodic", "periodic"],["periodic", "periodic"],["periodic", "periodic"]],
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_gamma_threshold = 0,
)

#   The mc synchrotron photon emitted will be stored here.
Species(
  name = "photon",
  position_initialization = "random",
  momentum_initialization = "cold",
  particles_per_cell = 0,
  mass = 0.0,
  charge = 0.0,
  number_density = 0.0,
  boundary_conditions = [["periodic", "periodic"],["periodic", "periodic"],["periodic", "periodic"]],
)
# ----------------------------------------------------------------------------------------
# Global parameters for the radiation reaction models
RadiationReaction(
    minimum_chi_discontinuous = 1e-4,
    minimum_chi_continuous = 1e4,
)

# ----------------------------------------------------------------------------------------
# Scalar diagnostics
DiagScalar(
    every = 1,
)

# ----------------------------------------------------------------------------------------
# Particle Binning

# # Energy distributions of synchrotron photons
DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    species =["photon"],
    axes = [["gamma", 0., gamma, 256]]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    species =["photon"],
    axes = [["gamma", 1e-3, gamma, 256, "logscale"]]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    species =["electron"],
    axes = [["chi", 1e-3, 1., 256, "logscale"]]
)
