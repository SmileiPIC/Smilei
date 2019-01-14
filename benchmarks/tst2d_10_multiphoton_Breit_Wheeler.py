# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Descrition:
# Interaction between a bunch of photons and a constant magnetic field
#
# Purpose:
# During the interaction,the photons will progressively decay into pairs
# of electron-positron via the multiphoton Breit-Wheeler process
# The radiation reaction of particles is not activated.
#
# Validation:
# - Propagation of macro-photons
# - multiphoton Breit-Wheeler pair creation process
#
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Importation

import math

# ----------------------------------------------------------------------------------------
# User defined parameters

c = 299792458                              # Speed of light
lambdar = 1e-6                             # Wavelength for normalization
wr = 2*math.pi*c/lambdar                   # Normalization angular frequency

Schwinger_E_field= 1.3E18                  # Schwinger electric field
Enorm = 3.2E12                             # Normalization electric field at lambda = 1

l0 = 2.0*math.pi                           # laser wavelength

chi = 20.                                  # Requested quantum parameter
B = 1000                                    # Requested magnetic field amplitude
gamma = chi* Schwinger_E_field/(Enorm*B)   # Initial normalized energy of the photons
Rsync = math.sqrt(gamma**2 - 1.)/B         # Synchrotron radius for a particle of energy gamma in a magnetic field B

Lx = 2.*Rsync                               # x length
Ly = 2.*Rsync                               # y length

n0 = 1e-5                                   # particle density

res = 128.                                  # nb of cells in one synchrotron radius

dt_factor = 0.9                             # factor on dt
dx = Rsync/res                              # space step
dy = Rsync/res                              # space step
dt  = 1./math.sqrt(1./(dx*dx) + 1./(dy*dy)) # timestep (CFL)
dt *= dt_factor

Tsim = 1000*dt                              # duration of the simulation

pusher = "boris"                            # dynamic type

part_cond = "periodic"                      # Particle and photon boundary conditions
field_cond = ['periodic','periodic']        # Field boundary conditions

# ----------------------------------------------------------------------------------------
# User-defined functions

# Density profile for inital location of the photons
def n0_photon(x,y):
        if ((x-0.5*Lx)**2 + (y-0.5*Ly)**2 <= (5*dx)**2):
                return n0
        else:
                return 0.

# Density profile for inital location of the electrons
def n0_electron(x,y):
                return 0.

# Density profile for inital location of the positrons
def n0_positron(x,y):
                return 0.

# ----------------------------------------------------------------------------------------
# Main parameters

Main(
    geometry = "2Dcartesian",

    interpolation_order = 4 ,

    cell_length = [dx,dy],
    grid_length  = [Lx,Ly],

    number_of_patches = [4,4],

    time_fields_frozen = Tsim,

    timestep = dt,
    simulation_time = Tsim,

    EM_boundary_conditions = [field_cond, field_cond],

    random_seed = smilei_mpi_rank,

    reference_angular_frequency_SI = wr
)

# ----------------------------------------------------------------------------------------
# External field

ExternalField(
    field = "Bz",
    profile = constant(B)
)

# ----------------------------------------------------------------------------------------
# List of species

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "rectangular",
    particles_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = pusher,
    boundary_conditions = [[part_cond], [part_cond]],
)

Species(
    name = "positron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = n0_positron,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = pusher,
    boundary_conditions = [[part_cond], [part_cond]],
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 128,
    c_part_max = 1.0,
    mass = 0,
    charge = 0.,
    number_density = n0_photon,
    mean_velocity = [0.0 ,-gamma, 0.0],
    temperature = [0.],
    pusher = "norm",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = [[part_cond], [part_cond]],
)

# ----------------------------------------------------------------------------------------
# QED parameters

RadiationReaction(
    minimum_chi_discontinuous = 1e-3,
    table_path = "./"
)

MultiphotonBreitWheeler(
    table_path = "./"
)

# ----------------------------------------------------------------------------------------
# Diagnostics

DiagScalar(
    every = 100,
    vars=['Uelm','Ukin','Utot','Uexp','Ubal',
          'Urad',
          'UmBWpairs',
          'Ukin_electron',
          'Ukin_photon',
          'Ukin_positron',
          'Ntot_electron',
          'Ntot_photon',
          'Ntot_positron']
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1000,
    time_average = 1,
    species = ["electron"],
    axes = [ ["gamma",    0.,  gamma,  50] ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1000,
    time_average = 1,
    species = ["positron"],
    axes = [ ["gamma",    0.,  gamma,  50] ]
)
