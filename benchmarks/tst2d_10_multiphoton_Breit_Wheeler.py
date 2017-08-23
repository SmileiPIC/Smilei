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
# Main parameters

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

Tsim = 2000*dt                              # duration of the simulation

pusher = "norm"                             # dynamic type

part_cond = "none"                          # Particle and photon boundary conditions
field_cond = ['periodic','periodic']        # Field boundary conditions

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
    geometry = "2d3v",

    interpolation_order = 4 ,

    cell_length = [dx,dy],
    sim_length  = [Lx,Ly],

    number_of_patches = [4,4],

    time_fields_frozen = Tsim,

    timestep = dt,
    sim_time = Tsim,

    bc_em_type_x = field_cond,
    bc_em_type_y = field_cond,

    random_seed = 0,

    referenceAngularFrequency_SI = wr
)

# ----------------------------------------------------------------------------------------
# External field

ExtField(
    field = "Bz",
    profile = constant(B)
)

# ----------------------------------------------------------------------------------------
# List of species

Species(
    species_type = "electron",
    initPosition_type = "random",
    initMomentum_type = "rectangular",
    n_part_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    dynamics_type = pusher,
    bc_part_type_xmin  = part_cond,
    bc_part_type_xmax  = part_cond,
    bc_part_type_ymin = part_cond,
    bc_part_type_ymax = part_cond,
    isTest = False
)

Species(
    species_type = "positron",
    initPosition_type = "random",
    initMomentum_type = "cold",
    n_part_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = n0_positron,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    dynamics_type = pusher,
    bc_part_type_xmin  = part_cond,
    bc_part_type_xmax  = part_cond,
    bc_part_type_ymin = part_cond,
    bc_part_type_ymax = part_cond,
    isTest = False
)

Species(
    species_type = "photon",
    initPosition_type = "random",
    initMomentum_type = "cold",
    n_part_per_cell = 128,
    c_part_max = 1.0,
    mass = 0,
    charge = 0.,
    nb_density = n0_photon,
    mean_velocity = [0.0 ,-gamma, 0.0],
    temperature = [0.],
    dynamics_type = "norm",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    bc_part_type_xmin  = part_cond,
    bc_part_type_xmax  = part_cond,
    bc_part_type_ymin = part_cond,
    bc_part_type_ymax = part_cond,
    isTest = False
)

# ----------------------------------------------------------------------------------------
# QED parameters

RadiationReaction(
    chipa_disc_min_threshold = 1e-2,
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

DiagParticles(
    output = "density",
    every = 2000,
    time_average = 1,
    species = ["electron"],
    axes = [ ["gamma",    0.,  gamma,  50] ]
)

DiagParticles(
    output = "density",
    every = 2000,
    time_average = 1,
    species = ["positron"],
    axes = [ ["gamma",    0.,  gamma,  50] ]
)
