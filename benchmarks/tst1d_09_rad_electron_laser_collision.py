# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Benchmark for the radiation in the collision of
# a GeV electron bunch with a counter-propagatin circularly polarized wave
#
# In this tests case, an electron bunch is initialized per radiation
# loss models at the same positions with an energy of 1 GeV near the right
# boundary of the box. They propagate to the left of the box where a circularly
# polarized laser plane wave is injected. This wave has an hyper-guassian
# profile of wavelength \lambda.
#
# Validation:
# - Continuous radiation reaction model: Landau-Lifshitz
#   with and without quantum corrections
# - Niel stochastic radiation reaction model
# - Monte-Carlo radiation reaction model without photon creation
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Importation

import math
import datetime

# ----------------------------------------------------------------------------------------
# User-defined parameters

c = 299792458
electron_mass = 9.10938356e-31
electron_charge = 1.60217662e-19
lambdar = 1e-6                          # Normalization wavelength
wr = 2*math.pi*c/lambdar                # Normalization frequency

l0 = 2.0*math.pi                        # laser wavelength
t0 = l0                                 # optical cicle
Lx = 30*l0                              # Domain length

n0 = 1e-5                               # bunch density

Tsim = 50.*t0                           # duration of the simulation
resx = 128.                             # nb of cells in one laser wavelength

dx = l0/resx                            # space step
dt  = 0.95 * dx                         # timestep (0.95 x CFL)

start = 0                               # Laser start
fwhm = 10*t0                            # Gaussian time fwhm
duration = 50*t0                        # Laser duration
center = duration*0.5                   # Laser profile center
order = 4                               # Laser order

gamma = 1000./0.511                     # Electron bunch initial energy
v = math.sqrt(1 - 1./gamma**2)          # electron bunch initial velocity

pusher = "vay"                         # type of pusher
radiation_list = ["Landau-Lifshitz","corrected-Landau-Lifshitz","Niel","Monte-Carlo"] # List of radiation models
species_name_list = ["LL","CLL","Niel","MC"]            # Name of the species

datetime = datetime.datetime.now()
random_seed = datetime.microsecond

# ----------------------------------------------------------------------------------------
# User-defined functions

# Density profile for inital location of the particles
def n0_(x):
        if (Lx - 10*dx < x < Lx - dx):
                return n0
        else:
                return 0.

# ----------------------------------------------------------------------------------------
# Namelists

Main(
    geometry = "1Dcartesian",

    interpolation_order = 4 ,

    cell_length = [dx],
    grid_length  = [Lx],

    number_of_patches = [16],

    timestep = dt,
    simulation_time = Tsim,

    EM_boundary_conditions = [['silver-muller']],

    reference_angular_frequency_SI = wr,

    random_seed = random_seed

)

# ----------------------------------------------------------------------------------------
# Laser definition
LaserPlanar1D(
    box_side         = "xmin",
    a0              = 100.,
    omega           = 1.,
    polarization_phi = 0.,
    ellipticity     = 1,
    time_envelope  = tgaussian(start=start,duration=duration,
                               fwhm=fwhm,
                               center=center,
                               order=order)
)

# ----------------------------------------------------------------------------------------
# Species

# Loop to create all the species
# One species per radiation implementations
for i,radiation in enumerate(radiation_list):

    Species(
        name = "electron_" + species_name_list[i],
        position_initialization = "regular",
        momentum_initialization = "cold",
        particles_per_cell = 64,
        c_part_max = 1.0,
        mass = 1.0,
        charge = -1.0,
        charge_density = n0_,
        mean_velocity = [-v, 0.0, 0.0],
        temperature = [0.],
        pusher = pusher,
        radiation_model = radiation,
        time_frozen = 29*t0,
        boundary_conditions = [
            ["periodic", "periodic"],
        ],
    )

# ----------------------------------------------------------------------------------------
# Radiation parameters

RadiationReaction(
    minimum_chi_discontinuous = 1e-3,
    #table_path = "./"
)

# ----------------------------------------------------------------------------------------
# Scalar Diagnostics

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
            ["x", 0., Lx, 1000],
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
            ["x", 0., Lx, 1000],
        ]
    )


for i,radiation in enumerate(radiation_list):
    # Chi-distribution
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = [5000,6500,100],
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["chi", 1e-3, 1., 256,"logscale"],
        ]
    )

for i,radiation in enumerate(radiation_list):
    # Gamma-distribution
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = [5000,6500,100],
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["gamma", 1., 1.1*gamma, 256,"logscale"],
        ]
    )
