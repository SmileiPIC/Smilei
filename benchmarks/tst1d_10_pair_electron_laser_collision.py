# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Descrition:
# Head-on collision between an electron bunch and a counter-propagating laser wave.
# The electron beam has an initial energy of 4 GeV.
# The laser has an intensity close to 10^23 W/cm^2
#
# Purpose:
# During the interaction,the electrons will radiate high-energy photons
# that in turn will decay onto pairs via the multiphoton Breit-Wheeler process.
#
# Validation:
# - Radiation reaction Monte-Carlo process
# - Emission and propagation of macro-photons
# - multiphoton Breit-Wheeler pair creation process
#
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Importation

import math

# ----------------------------------------------------------------------------------------
# Main parameters

c = 299792458                           # Speed of light
lambdar = 1e-6                          # Wavelength for normalization
wr = 2*math.pi*c/lambdar                # Normalization angular frequency

l0 = 2.0*math.pi                        # laser wavelength
t0 = l0                                 # optical cicle

Lx = 30*l0                              # Domain length

n0 = 1e-5                               # beam density

Tsim = 50.*t0                           # duration of the simulation
resx = 100.                             # nb of cells in one laser wavelength

dx = l0/resx                            # space step
dt  = 0.95 * dx                 		# timestep (0.95 x CFL)

start = 0                               # Laser start
fwhm = 10*t0                            # Gaussian time fwhm
duration = 50*t0                        # Laser duration
center = duration*0.5                   # Laser profile center
order = 4                               # Laser order

gamma0 = 4000./0.511                    # Initial electron beam normalized energy (Lorentz factor)
v0 = math.sqrt(1 - 1./gamma0**2)         # Initial velocity

radiation = "Monte-Carlo"               # Radiation algorithm
pusher = "vay"                          # dynamic type

part_cond = "supp"                      # Particle boundary conditions
field_cond = ['silver-muller']          # Field boundary conditions

# Density profile for inital electron location
def n0_electron(x):
        if (0.98*Lx < x < 0.99*Lx):
                return n0
        else:
                return 0.

# Density profile for inital positron location
def n0_positron(x):
                return 0.

# Density profile for inital photon location
def n0_photon(x):
                return 0.

# ----------------------------------------------------------------------------------------
# Main parameters

Main(
    geometry = "1d3v",

    interpolation_order = 4 ,

    cell_length = [dx],
    sim_length  = [Lx],

    number_of_patches = [4],

    timestep = dt,
    sim_time = Tsim,

    bc_em_type_x = ['silver-muller'],

    referenceAngularFrequency_SI = wr,

    random_seed = 0

)

# ----------------------------------------------------------------------------------------
# Laser profile

LaserPlanar1D(
    boxSide         = "xmin",
    a0              = 270.,
    omega           = 1.,
    polarizationPhi = 0.,
    ellipticity     = 1,
    time_envelope  = tgaussian(start=start,duration=duration,
                               fwhm=fwhm,
                               center=center,
                               order=order)
)

# ----------------------------------------------------------------------------------------
# List of species

Species(
    species_type = "electron",
    initPosition_type = "regular",
    initMomentum_type = "cold",
    n_part_per_cell = 32,
    c_part_max = 1.,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = [-v0, 0.0, 0.0],
    temperature = [0.],
    dynamics_type = pusher,
    radiation_model = radiation,
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 2,
    time_frozen = 29*t0,
    bc_part_type_xmin  = part_cond,
    bc_part_type_xmax  = part_cond,
    isTest = False
)

Species(
    species_type = "positron",
    initPosition_type = "random",
    initMomentum_type = "cold",
    n_part_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_positron,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.],
    dynamics_type = pusher,
    radiation_model = radiation,
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 2,
    time_frozen = 29*t0,
    bc_part_type_xmin  = part_cond,
    bc_part_type_xmax  = part_cond,
    isTest = False
)

Species(
    species_type = "photon",
    initPosition_type = "random",
    initMomentum_type = "cold",
    n_part_per_cell = 0,
    c_part_max = 20.0,
    mass = 0,
    charge = 0.,
    nb_density = n0_photon,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    dynamics_type = "norm",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    bc_part_type_xmin  = part_cond,
    bc_part_type_xmax  = part_cond,
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
    vars=['Uelm','Ukin','Utot',
          'Uexp',
          'Ubal',
          'Urad',
          'Ukin_electron',
          'Ukin_positron',
          'Ukin_photon',
          'Ntot_electron',
          'Ntot_positron',
          'Ntot_photon']
)
