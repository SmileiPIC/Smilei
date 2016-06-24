# SIMULATION PARAMETERS
# dim: Geometry of the simulation.
#      1d3v = cartesian grid with 1d in space + 3d in velocity.
#      2d3v = cartesian grid with 2d in space + 3d in velocity.
#      3d3v = cartesian grid with 3d in space + 3d in velocity.
#      2drz = cylindrical (r,z) grid with 3d3v particles.
#
import math
dim = "2d3v"

# order of interpolation
interpolation_order = 2

# SIMULATION TIME
# timestep: time step of a single iteration in 1/w0.
# sim_time: duration of the simulation in 1/w0 .

timestep = 0.124
sim_time = 150

# SIMULATION BOX : for all space directions (use vector)
#n_space_win_x: Activates moving window. Number of cell in the total computation domain in the x direction.
nspace_win_x = 896
# cell_length: length of one cell in c/w0.
cell_length  = [0.125, 3.] 
# sim_length: Total length of the simulation in c/w0 .
sim_length = [ cell_length[0]*nspace_win_x,  120.]  

#Moving window parameters
#t_move_win: Starting time of the moving window in 1/w0.
t_move_win = nspace_win_x*cell_length[0]
#vx_win: Velocity of the moving window in c.
vx_win = 0.9997

# ----------------
# Dynamic Load Balancing
# ----------------
#Number of patches must divide number of cells (sim_length/cell_length).
number_of_patches = [128, 8]
#Rebalancing frequency
balancing_freq = 20
coef_cell = 1.
coef_frozen = 0.1
#clrw: width of a cluster in number of cell. Warning: clrw must divide nspace_win_x.
clrw = nspace_win_x/number_of_patches[0]

# RANDOM seed 
# this is used to randomize the random number generator.
random_seed = 0

# DEFINE ALL SPECIES

# species_type: ion, electron, positron, test ...
# initialization_type: regular, cold or (isotrop) Maxwell-Juettner distribution# n_part_per_cell: number of particle-per-cell
# c_part_max: factor on the memory reserved for the total number of particles.
# mass: particle mass in units of the electron mass.
# charge: particle charge in units of e (-e is the electron charge).
# density: species density in units of the normalization density.
# mean_velocity: mean velocity of the species (3D vector) in units of the light velocity.
# temperature: temperature of the species in units of m_e c^2.
# dynamics_type: species type of dynamics = norm or rrLL.
# time_frozen: time during which the particles are frozen in units of the normalization time.
# radiating: boolean, if true incoherent radiation are calculated using the Larmor formula .
#dens_length_x: Length of the plateau - Length of left slope - Length of right slope.
#dens_length_y: Length of the plateau - Length of left slope - Length of right slope.

Species(
    species_type = "proton",
    initPosition_type = "regular",
    initMomentum_type = "cold",
    ionization_model = "none",
    n_part_per_cell = 16, 
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    dynamics_type = "norm",
    time_frozen = 100000.,
    radiating = False,
    bc_part_type_west  = "supp",
    bc_part_type_east  = "supp",
    bc_part_type_south = "supp",
    bc_part_type_north = "supp"
)

Species( 
    species_type = "electron",
    initPosition_type = "regular",
    initMomentum_type = "cold",
    n_part_per_cell = 16,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    dynamics_type = "norm",    
    time_frozen = 0.0,
    radiating = False,
    bc_part_type_west = "supp",
    bc_part_type_east = "supp",
    bc_part_type_south ="stop",
    bc_part_type_north ="stop"
)

# ----------------
# Boundary Conditions for fields
# ----------------
bc_em_type_x = ["silver-muller","silver-muller"]
bc_em_type_y = ["silver-muller","silver-muller"]

# Laser properties
fwhm_field=19.80
LaserGaussian2D(
    boxSide         = "west",
    a0              = 2.,
    focus           = [0., sim_length[1]/2.],
    waist           = 26.16,
    time_envelope   = tgaussian(center=fwhm_field*2**0.5, fwhm=fwhm_field)
)


# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

# print_every (on screen text output) 
print_every = 100

DiagFields(
    every = 100,
    fields = ['Ex','Ey','Rho_electron','Rho_proton','Jx_electron']
)

#Write a restart file every dump_step iterations. 0 for no dump.
dump_step = 0
#Write a restart file every dump_minutes minutes. 0.0 for no dump.
dump_minutes = 0.0

# DIAG ON SCALARS
# every = number of time-steps between each output
#
DiagScalar(every = 100)

exit_after_dump = False
