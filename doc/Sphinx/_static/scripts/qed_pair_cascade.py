# _____________________________________________________________________________
#
# QED cascade 3D
# _____________________________________________________________________________

import numpy as np

# ______________________________________________________________________________
# General parameters

c = 299792458
lambdar = 1e-6
reference_angular_frequency_SI = 2*math.pi*c/lambdar

l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle

#Domain size
Lx = 4*l0
Ly = 2.*l0
Lz = 2.*l0

# particle density
n0 = 1.e-5

# initial nb of particles in the simulation
npart0  = 32

# discretization
resx = 48.
resy = 48.
resz = 48.
dx = l0/resx                            # space step
dy = l0/resx
dz = l0/resz
dt  = 0.95 / np.sqrt( 1/dx**2 + 1/dy**2 + 1/dz**2)

# Patch decomposition
number_of_patches = [int(Lx/dx)/8,int(Ly/dy)/8,int(Lz/dz)/8]

# duration of the simulation
Tsim = 2*Lx

pusher = "boris"                         # dynamic type
radiation = "Monte-Carlo"               # Radiation algorithm

time_frozen = 0.45*Lx

particles_per_cell = 1

interpolation_order = 2

radiation_photon_gamma_threshold = 20.

particle_boundary_conditions = [["remove", "remove"],["periodic","periodic"],["periodic","periodic"]]

EM_boundary_conditions = [ ["silver-muller","silver-muller"],["periodic","periodic"],["periodic","periodic"]   ]

# Period for the diags
diag_every = 100

# Laser definition
theta_pol = np.pi/4.    # polarisation parameter (Warning! different than standard definition in SMILEI= pi/4 for CP)
beta_pol  = -1.           # additional polarization parameter (Check! is it left/right handed for x<0 propag.)
                          # -1: left-handed;    +1: right-handed (for x>0 propagation)
tramp = t0/4.               # laser linear ramp time
a0 = 1000.

def laser_env(t):
    if (t<tramp):
        return t/tramp
    else:
        return 1.

def By_p(y,z,t):
    return -beta_pol*np.sin(theta_pol) * a0 * laser_env(t) * np.sin(t)

def Bz_p(y,z,t):
    return a0 * np.cos(theta_pol) * laser_env(t) * np.cos(t)

def By_m(y,z,t):
    return -beta_pol*np.sin(theta_pol) * a0 * laser_env(t) * np.sin(t+np.pi)

def Bz_m(y,z,t):
    return a0 * np.cos(theta_pol) * laser_env(t) * np.cos(t+np.pi)

# Initial list of particles: seed
particle_initialization = np.zeros((4,npart0))
particle_initialization[0,:] = np.multiply(np.ones(npart0),Lx*0.5+0.5*dx)
particle_initialization[1,:] = np.linspace(0.45*Ly,0.55*Ly,npart0)
particle_initialization[2,:] = np.linspace(0.45*Lz,0.55*Lz,npart0)
particle_initialization[3,:] = np.multiply(np.ones(npart0),n0/npart0)

# ______________________________________________________________________________
# Namelists

Main(
    geometry = "3Dcartesian",

    interpolation_order = interpolation_order,

    cell_length = [dx,dy,dz],
    grid_length  = [Lx,Ly,Lz],
    number_of_patches = number_of_patches,

    timestep = dt,
    simulation_time = Tsim,

    EM_boundary_conditions = EM_boundary_conditions,

    print_every = 50,

    random_seed = smilei_mpi_rank,

    clrw = 1,

    reference_angular_frequency_SI = reference_angular_frequency_SI,

    patch_arrangement = "linearized_ZYX",

)

### Laser from the left
Laser(
    box_side = "xmin",
    space_time_profile = [ By_p, Bz_p ]
)
### Laser from the right
Laser(
    box_side = "xmax",
    space_time_profile = [ By_m, Bz_m ]
)

Species(
    name = "electron",
    position_initialization = particle_initialization,
    momentum_initialization = "cold",
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = pusher,
    radiation_model = radiation,
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = radiation_photon_gamma_threshold,
    boundary_conditions = particle_boundary_conditions,
    time_frozen = time_frozen,
    merging_method = "vranic_spherical",
    # Merging properties
    merge_every = 1,
    merge_min_particles_per_cell = 4,
    merge_discretization_scale = "log",
    merge_min_momentum = 1e-2,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [6,6,6],
)

Species(
    name = "positron",
    position_initialization = particle_initialization,
    momentum_initialization = "cold",
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = pusher,
    radiation_model = radiation,
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = radiation_photon_gamma_threshold,
    boundary_conditions = particle_boundary_conditions,
    time_frozen = time_frozen,
    # Merging properties
    merging_method = "vranic_spherical",
    merge_every = 1,
    merge_min_particles_per_cell = 4,
    merge_discretization_scale = "log",
    merge_min_momentum = 1e-2,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [6,6,6],
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 0,
    charge = 0,
    number_density = n0_photon,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "norm",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = particle_boundary_conditions,
    # Mergin properties
    merging_method = "vranic_spherical",
    merge_every = 1,
    merge_min_particles_per_cell = 4,
    merge_discretization_scale = "log",
    merge_min_momentum = 20,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [6,6,6],
)

Vectorization(
    mode = "on",
)

RadiationReaction(
    minimum_chi_discontinuous = 1e-2,
    table_path = "<path to Smilei>/databases"
)

MultiphotonBreitWheeler(
    table_path = "<path to Smilei>/databases"
)

DiagScalar(
    every = 10,
    vars=['Uelm','Ukin','Utot','Uexp','Ubal',
          'Urad',
          'Ukin_bnd',
          'UmBWpairs',
          'Ukin_electron',
          'Ukin_photon',
          'Ukin_positron',
          'Ntot_electron',
          'Ntot_photon',
          'Ntot_positron',
          'Dens_electron',
          'Dens_positron',
          'Dens_photon'],
)

DiagPerformances(
    every = 10,
)

DiagFields(
    every = diag_every,
    time_average = 1,
    fields = ["Ey", "Ez", "By", "Bz"]
)

species_list = ["electron","positron","photon"]

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_every,
		time_average = 1,
		species = [species],
		axes = [
		    ["x", 0, Lx, int(Lx/dx)],
		    ["y", 0, Ly, int(Ly/dy)],
		    ["z", 0, Lz, int(Lz/dz)],
		]
	)

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_every,
		time_average = 1,
		species = [species],
		axes = [
		    ["px", -500, 500, 512],
		    ["py", -500, 500, 512],
		]
	)

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_every,
		time_average = 1,
		species = [species],
		axes = [
		    ["gamma", 0, 500, 256],
		]
	)
