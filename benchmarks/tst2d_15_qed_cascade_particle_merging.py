# ----------------------------------------------------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
#
# Input file to test the particle merging in 2D
#
# This scenario simulate the QED cascade induced in the collision of 2
# lineary polarized plane waves.
# ----------------------------------------------------------------------------------------

# ______________________________________________________________________________
# Useful libs

import math

# _____________________________________________________________________________
# Main parameters

# speed of light
c = 299792458
# Reference wavelength used as well for normalization (m)
wavelength = 1e-6
# Angular frequency
angular_frequency = 2*math.pi*c/wavelength

# Normalized laser wavelength
l0 = 2.0*math.pi
# Optical time
t0 = l0

# Domain dimensions
Lx = 5*l0
Ly = 0.5*l0

# Seed density
n0 = 1.e-5

# Discretization of l0
resx = 48.
dx = l0/resx
dy = Ly/24.

# Timestep (0.95 x CFL)
dt  = 0.95/math.sqrt( 1./dx**2 + 1./dy**2 )

# Patch distribution
number_of_patches = [4,4]

# Laser parameters
start = 0                               # Laser start
fwhm = 30*t0                            # Gaussian time fwhm
duration = 90*t0                        # Laser duration
center = duration*0.5                   # Laser profile center
ramp = 30*t0
a0 = 500.

# duration of the simulation
simulation_time = 0.5*Lx + duration

# Frozen time
time_frozen = 0.49*Lx

# Period of large output
diag_time_selection = [2000*dt, 500]

# Few merging parameters
merge_time_selection = [int(time_frozen/dt), 5]
merge_min_particles_per_cell = 8

# Boundary conditions
particle_boundary_conditions = [["remove", "remove"],["periodic", "periodic"]]
EM_boundary_conditions = [["silver-muller","silver-muller"],["periodic","periodic"] ]

# Density profile for inital location of the particles
def n0_electron(x,y):
    if (0.5*(Lx-dx)<x<0.5*(Lx+dx)): #and 0.49*Ly<y<0.51*Ly and 0.49*Lz<z<0.51*Lz
      return n0
    else:
      return 0;

# _____________________________________________________________________________
# Namelist

Main(

    geometry = "2Dcartesian",
    interpolation_order = 2,
    cell_length = [dx,dy],
    grid_length  = [Lx,Ly],
    number_of_patches = number_of_patches,
    timestep = dt,
    simulation_time = simulation_time,
    EM_boundary_conditions = EM_boundary_conditions,
    print_every = 100,
    random_seed = 0,
    clrw = 1,
    patch_arrangement = 'linearized_XY',
    reference_angular_frequency_SI = angular_frequency,

)

LaserPlanar1D(
    box_side         = "xmin",
    a0              = a0,
    omega           = 1.,
#    focus           = [Lx*0.5,Ly*0.5],
#    waist           = 1e9,
#    incidence_angle = [0., 0.],
    polarization_phi = 0.,
    ellipticity     = 0,
    time_envelope  = tgaussian(start=start,duration=duration,fwhm=fwhm,center=center,order=4)
#    time_envelope  = tsin2plateau(start=start,plateau=duration,slope1=ramp,slope2=ramp)
)

LaserPlanar1D(
    box_side         = "xmax",
    a0              = a0,
    omega           = 1.,
#    focus           = [Lx*0.5],
#    waist           = 1e9,
#    incidence_angle = [0., 0.],
    polarization_phi = 0.,
    ellipticity     = 0,
    time_envelope  = tgaussian(start=start,duration=duration,fwhm=fwhm,center=center,order=4)
#    time_envelope  = tsin2plateau(start=start,plateau=duration,slope1=ramp,slope2=ramp)
)

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "boris",
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 20,
    boundary_conditions = particle_boundary_conditions,
    time_frozen = time_frozen,
    
    # merging parameters
    merging_method = "vranic_spherical",
    merge_every = merge_time_selection,
    merge_min_particles_per_cell = merge_min_particles_per_cell,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [16,16,1],
    merge_discretization_scale = "log",
)

Species(
    name = "positron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = 0,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "boris",
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 20.,
    boundary_conditions = particle_boundary_conditions,
    time_frozen = time_frozen,
    
    # merging parameters
    merging_method = "vranic_spherical",
    merge_every = merge_time_selection,
    merge_min_particles_per_cell = merge_min_particles_per_cell,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [16,16,1],
    merge_discretization_scale = "log",
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 0,
    charge = 0,
    number_density = 0,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "norm",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = particle_boundary_conditions,
    
    # Merging parameters
    merging_method = "vranic_spherical",
    merge_every = merge_time_selection,
    merge_min_particles_per_cell = merge_min_particles_per_cell,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [8,8,1],
    merge_accumulation_correction = False,
    merge_discretization_scale = "log",
)

Vectorization(
    mode = "on",
)

RadiationReaction(
    minimum_chi_discontinuous = 1e-2,
)

MultiphotonBreitWheeler(
    #table_path = "./"
)


DiagScalar(
    every = 100,
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

# DiagPerformances(
#    every = 10,
#    flush_every = 100,
#    # patch_information = True,
# )

species_list = ["electron","positron","photon"]


for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_time_selection,
		time_average = 1,
		species = [species],
		axes = [
		    ["gamma", 1., a0, 128, "logscale"],
		]
	)

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = diag_time_selection,
		time_average = 1,
		species = [species],
		axes = [
		    ["chi", 1e-3, 5., 128, "logscale"],
		]
	)
