# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Synchrotron: radiation loss of electrons rotating in a constant
#              magnetic field
#
# In this tests case, an electron bunch is initialized per radiation
# loss models at the same positions. The magnetic field and the initial energy
# is computed so that the initial quantum parameter is equal to 0.1
#
# Validation:
# - Corrected Landau-Lifshitz radiation loss
# - Niel stochastic radiation loss
# - Monte-Carlo model
# - Species scalar diagnostics
# - External fields
# - Particle binning diagnostics for the quantum parameter
# ----------------------------------------------------------------------------------------

import math
import datetime

# ----------------------------------------------------------------------------------------
# User defined parameters

c = 299792458
electron_mass = 9.10938356e-31
electron_charge = 1.60217662e-19
lambdar = 1e-6                                    # Normalization wavelength
wr = 2*math.pi*c/lambdar                          # Normalization frequency

Schwinger_E_field= 1.3E18                         # Schwinger electric field
Enorm = electron_mass*wr*c/electron_charge        # Normalization electric field at lambda = 1e-6

l0 = 2.0*math.pi                                  # laser wavelength

chi = 0.1                                         # Initial quantum parameter
B = 90.                                           # Magnetic field strength
gamma = chi* Schwinger_E_field/(Enorm*B)          # Initial gamma factor
Rsync = math.sqrt(gamma**2 - 1.)/B                # Synchrotron radius without radiation
v = math.sqrt(1.-1./gamma**2)                     # Initial particle velocity

Lx = 4.*Rsync                                     # Dimension of the box in x
Ly = 4.*Rsync                                     # Dimension of the box in y

n0 = 1e-5                                         # particle density

res = 128.                                        # nb of cells in one synchrotron radius

dx = Rsync/res                                    # space step
dy = Rsync/res                                    # space step
dt  = 1./math.sqrt(1./(dx*dx) + 1./(dy*dy))       # timestep given by the CFL

dt_factor = 0.9                                   # factor on dt
dt *= dt_factor                                   # timestep used for the simulation

Tsim = 5000*dt/dt_factor                          # duration of the simulation

pusher = "vay"                                    # type of pusher
radiation_list = ["corrected-Landau-Lifshitz","Niel","Monte-Carlo"]    # List of radiation models for species
species_name_list = ["CLL","Niel","MC"]           # List of names for species

datetime = datetime.datetime.now()
random_seed = datetime.microsecond

# ----------------------------------------------------------------------------------------
# User-defined functions

# Density profile for inital location of the particles
def n0_(x,y):
        if ((x-0.75*Lx)**2 + (y-0.5*Ly)**2 <= 0.5*dx):
                return n0
        else:
                return 0.

# ----------------------------------------------------------------------------------------
# Namelists

Main(
    geometry = "2Dcartesian",

    interpolation_order = 2,

    cell_length = [dx,dy],
    grid_length  = [Lx,Ly],

    number_of_patches = [4,4],
    time_fields_frozen = Tsim,

    timestep = dt,
    simulation_time = Tsim,

    EM_boundary_conditions = [['periodic'],['periodic']],

    random_seed = smilei_mpi_rank,

    reference_angular_frequency_SI = wr

)

Vectorization(
    mode = "normal",
)

# ----------------------------------------------------------------------------------------
# Initialization of the constant external field
ExternalField(
    field = "Bz",
    profile = constant(B)
)

# ----------------------------------------------------------------------------------------
# Loop to create all the species
# One species per radiation implementations
for i,radiation in enumerate(radiation_list):

    Species(
        name = "electron_" + species_name_list[i],
        position_initialization = "regular",
        momentum_initialization = "cold",
        particles_per_cell = 16,
        c_part_max = 1.0,
        mass = 1.0,
        charge = -1.0,
        charge_density = n0_,
        mean_velocity = [0.0 ,v, 0.0],
        temperature = [0.],
        pusher = pusher,
        radiation_model = radiation,
        boundary_conditions = [
            ["periodic", "periodic"],
            ["periodic", "periodic"],
        ],
    )

# ----------------------------------------------------------------------------------------
# Global parameters for the radiation reaction models
RadiationReaction(
    chipa_disc_min_threshold = 1e-3,
    table_path = "./"
)

# ----------------------------------------------------------------------------------------
# Scalar diagnostics
DiagScalar(
    every = 100
)

# ----------------------------------------------------------------------------------------
# Particle Binning

# Loop to create all the species particle binning diagnostics
# One species per radiation implementations
for i,radiation in enumerate(radiation_list):

    # Weight spatial-distribution
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = 500,
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["x", 0., Lx, 200],
            ["y", 0., Ly, 200],
        ]
    )


for i,radiation in enumerate(radiation_list):
    # Weight x chi spatial-distribution
    DiagParticleBinning(
        deposited_quantity = "weight_chi",
        every = 500,
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["x", 0., Lx, 200],
            ["y", 0., Ly, 200],
        ]
    )


for i,radiation in enumerate(radiation_list):
    # Chi-distribution
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = 500,
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["chi", 1e-3, 1., 1000,"logscale"],
        ]
    )
