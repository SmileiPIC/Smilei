# _____________________________________________________________________________
#
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
#
# Magnetic shower 3d
#
# _____________________________________________________________________________

import math as math
import numpy as np

# _____________________________________________________________________________
# Main parameters

c = 299792458
lambdar = 1e-6
wr = 2*math.pi*c/lambdar

l0 = 2.0*math.pi              # laser wavelength

Schwinger_E_field= 1.3E18                 # Schwinger electric field
E_norm = 3.2E12                            # Normalization electric field at lambda = 1

chi = 20.
B_field_amplitude = 1000
gamma = chi* Schwinger_E_field/(E_norm*B_field_amplitude)
Rsync = math.sqrt(gamma**2 - 1.)/B_field_amplitude
v = math.sqrt(1.-1./gamma**2)

print("Gamma: {}".format(gamma))
print("Velocity: {}".format(v))

B_field_direction = np.array([0,0,1])
B_field_direction = B_field_direction / np.linalg.norm(B_field_direction)
B_field_vector = B_field_amplitude* B_field_direction

mean_velocity = [0,v,0]

diag_every = 50

Lx = 1.*Rsync
Ly = 1.*Rsync
Lz = 1.*Rsync

n0 = 1e-5                     # particle density

res = 16.                             # nb of cells in one synchrotron radius

particles_per_cell = 8

dt_factor = 0.05
dx = Rsync/res                            # space step
dy = Rsync/res                            # space step
dz = Rsync/res                            # space step
dt  = 1./math.sqrt(1./(dx*dx) + 1./(dy*dy)) # timestep (CFL)
dt *= dt_factor

diag_res = 1024
diag_min = -gamma*0.6
diag_max = gamma*0.6

simulation_time = 70*dt                 # duration of the simulation

merging_method = "vranic_spherical"
pusher = "vay"
radiation = "Monte-Carlo"

# Density profile for inital location of the particles
def n0_photon(x,y,z):
                return 0.


# Density profile for inital location of the particles
def n0_electron(x,y,z):
                return n0


# ______________________________________________________________________________
# Namelists

Main(
    geometry = "3Dcartesian",

    interpolation_order = 2 ,

    cell_length = [dx,dy,dz],
    grid_length  = [Lx,Ly,Lz],

    number_of_patches = [2,2,2],

    time_fields_frozen = simulation_time,

    timestep = dt,
    simulation_time = simulation_time,

    EM_boundary_conditions = [ ["periodic"] ],


    print_every = 10,

    reference_angular_frequency_SI = wr,

)

ExternalField(
    field = "Bx",
    profile = constant(B_field_vector[0])
)

ExternalField(
    field = "By",
    profile = constant(B_field_vector[1])
)

ExternalField(
    field = "Bz",
    profile = constant(B_field_vector[2])
)

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "rectangular",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = mean_velocity,
    temperature = [0.],
    pusher = pusher,
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_gamma_threshold = 10,
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
    # Merging parameters:
    merging_method = merging_method,
    merge_every = 5,
    merge_min_particles_per_cell = 16,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [16,16,16],
    merge_accumulation_correction = False,
)

Species(
    name = "positron",
    position_initialization = "random",
    momentum_initialization = "rectangular",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = n0_electron,
    mean_velocity = mean_velocity,
    temperature = [0.],
    pusher = pusher,
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_gamma_threshold = 10,
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
    # Merging parameters:
    merging_method = merging_method,
    merge_every = 5,
    merge_min_particles_per_cell = 16,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [16,16,16],
    merge_accumulation_correction = False,
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    c_part_max = 1.0,
    mass = 0,
    charge = 0.,
    number_density = n0_photon,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "norm",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
    # Merging parameters:
    merging_method = merging_method,
    merge_every = 5,
    merge_min_particles_per_cell = 8,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [16,16,16],
    merge_accumulation_correction = False,
)

Vectorization(
    mode = "on",
)

RadiationReaction(
    minimum_chi_discontinuous = 1e-2,
    table_path = "/home/mathieu/Documents/Codes/particle_merging/databases"
)

MultiphotonBreitWheeler(
    #table_path = "./"
    table_path = "/home/mathieu/Documents/Codes/particle_merging/databases"
)

DiagScalar(
    every = 1,
    vars=['Uelm','Ukin','Utot','Uexp','Ubal',
          'Urad',
          'UmBWpairs',
          'Ukin_electron',
          'Ukin_photon',
          'Ukin_positron',
          'Ntot_electron',
          'Ntot_photon',
          'Ntot_positron',
          'Dens_electron',
          'Dens_positron',
          'Dens_photon']
)

DiagPerformances(
    every = 1,
)

species_list = ["electron","positron","photon"]

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_every,
		time_average = 1,
		species = [species],
		axes = [
		    ["px", diag_min, diag_max, diag_res],
		    ["py", diag_min, diag_max, diag_res],
		]
	)

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_every,
		time_average = 1,
		species = [species],
		axes = [
		    ["px", diag_min, diag_max, diag_res],
		    ["pz", diag_min, diag_max, diag_res],
		]
	)

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_every,
		time_average = 1,
		species = [species],
		axes = [
		    ["py", diag_min, diag_max, diag_res],
		    ["pz", diag_min, diag_max, diag_res],
		]
	)

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_every,
		time_average = 1,
		species = [species],
		axes = [
		    ["gamma", 0, diag_max, 512],
		]
	)
