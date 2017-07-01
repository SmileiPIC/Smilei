# _____________________________________________________________________________
#
# Synchrotron: radiation loss of electrons rotating in a constant
#              magnetic field
#
# In this tests case, an electron bunch is initialized per radiation
# loss models at the same positions. The magnetic field and the initial energy
# is computed so that the initial quantum parameter is equal to 1.
#
# Validation:
# - Discontinuous radiation loss
# - Continuous radiation loss
# _____________________________________________________________________________

import math

# _____________________________________________________________________________
# Main parameters

c = 299792458
electron_mass = 9.10938356e-31
electron_charge = 1.60217662e-19
lambdar = 1e-6                # Normalization wavelength
wr = 2*math.pi*c/lambdar      # Normalization frequency

Schwinger_E_field= 1.3E18     # Schwinger electric field
Enorm = electron_mass*wr*c/electron_charge        # Normalization electric field at lambda = 1e-6

l0 = 2.0*math.pi              # laser wavelength

chi = 1.                      # Initial quantum parameter
B = 270.                      # Magnetic field strength
gamma = chi* Schwinger_E_field/(Enorm*B) # Initial gamma factor
Rsync = math.sqrt(gamma**2 - 1.)/B            # Synchrotron radius without radiation
v = math.sqrt(1.-1./gamma**2)                 # Initial particle velocity

Lx = 4.*Rsync
Ly = 4.*Rsync

n0 = 1e-5                              # particle density

res = 128.                             # nb of cells in one synchrotron radius

dx = Rsync/res                            # space step
dy = Rsync/res                            # space step
dt  = 1./math.sqrt(1./(dx*dx) + 1./(dy*dy)) # timestep (0.95 x CFL)

dt_factor = 0.9                           # factor on dt
dx = Rsync/res                            # space step
dy = Rsync/res                            # space step
dt  = 1./math.sqrt(1./(dx*dx) + 1./(dy*dy)) # timestep (CFL)
dt *= dt_factor

Tsim = 5000*dt/dt_factor                 # duration of the simulation

pusher = "vay"                         # dynamic type
radiation_list = ["Monte-Carlo","corrected-Landau-Lifshitz"]
species_name_list = ["disc","cont"]

# ______________________________________________________________________________
# Functions

# Density profile for inital location of the particles
def n0_(x,y):
        if ((x-0.75*Lx)**2 + (y-0.5*Ly)**2 <= 0.25*dx):
                return n0
        else:
                return 0.

# ______________________________________________________________________________
# Namelists

Main(
    geometry = "2d3v",

    interpolation_order = 4,

    cell_length = [dx,dy],
    sim_length  = [Lx,Ly],

    number_of_patches = [4,4],

    time_fields_frozen = Tsim,

    timestep = dt,
    sim_time = Tsim,

    bc_em_type_x = ['periodic'],
    bc_em_type_y = ['periodic'],

    random_seed = 0,

    referenceAngularFrequency_SI = wr

)

ExtField(
    field = "Bz",
    profile = constant(B)
)

# Loop to create all the species
# One species per radiation implementations
for i,radiation in enumerate(radiation_list):

    Species(
        species_type = "electron_" + species_name_list[i],
        initPosition_type = "centered",
        initMomentum_type = "cold",
        n_part_per_cell = 10,
        c_part_max = 1.0,
        mass = 1.0,
        charge = -1.0,
        charge_density = n0_,
        mean_velocity = [0.0 ,v, 0.0],
        temperature = [0.],
        dynamics_type = pusher,
        radiation_model = radiation,
        bc_part_type_xmin  = "none",
        bc_part_type_xmax  = "none",
        bc_part_type_ymin = "none",
        bc_part_type_ymax = "none",
        bc_part_type_zmin = "none",
        bc_part_type_zmax = "none",
#        track_every = 1,
#        track_flush_every = 1000,
#        isTest = False
    )

RadiationLoss(
    chipa_integfochi_min = 1e-4,
    chipa_integfochi_max = 1e1,
    integfochi_dim = 10,

    chipa_xip_min = 1e-4,
    chipa_xip_max = 1e1,
    xip_power = 4,
    xip_threshold = 1e-3,
    chipa_xip_dim = 10,
    chiph_xip_dim = 10,

    chipa_disc_min_threshold = 1e-2
)

DiagScalar(
    every = 100,
    vars=['Ukin_electron_disc',
          'Ukin_electron_cont',
          'Urad_electron_disc',
          'Urad_electron_cont']
)
