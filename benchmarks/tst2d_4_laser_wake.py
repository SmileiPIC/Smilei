# SIMULATION PARAMETERS
# dim: Geometry of the simulation.
#      1d3v = cartesian grid with 1d in space + 3d in velocity.
#      2d3v = cartesian grid with 2d in space + 3d in velocity.
#      3d3v = cartesian grid with 3d in space + 3d in velocity.
#      2drz = cylindrical (r,z) grid with 3d3v particles.
#

sim_units = "normalized"
dim = "2d3v"

# order of interpolation
interpolation_order = 2

# SIMULATION TIME
# timestep: time step of a single iteration in 1/w0.
# sim_time: duration of the simulation in 1/w0 .

timestep = 0.124
sim_time = 93

#Moving window
#n_space_win_x: Number of cell in the total computation domain.
#t_move_win: Starting time of the moving window in 1/w0.
#vx_win: Velocity of the moving window in c.
nspace_win_x = 896
t_move_win = 60.
vx_win = 0.9997

# SIMULATION BOX : for all space directions (use vector)
# cell_length: length of one cell in c/w0.
# sim_length: Total length of the simulation in c/w0 .

cell_length  = [0.125, 3.] 
sim_length = [ cell_length[0]*nspace_win_x,  120.]  
bc_em_type_x = ["silver-muller","silver-muller"]
bc_em_type_y = ["silver-muller","silver-muller"]

#Topology:
#number_of_procs: Number of MPI processes in each direction.
#clrw: width of a cluster in number of cell. Warning: clrw must divide nspace_win_x.
number_of_patches = [128, 8] 
balancing_freq = 150
clrw = 1

# PLASMA GEOMETRY
# plasma_geometry: string defining the plasma geometry.
# ***_length: characteristic length for a given geometry.
#


# RANDOM seed 
# this is used to randomize the random number generator.
random_seed = 0

# DEFINE ALL SPECIES

# species_type: ion, electron, positron, test ...
# initialization_type: regular, cold or (isotrop) Maxwell−Juettner distribution# n_part_per_cell: number of particle−per−cell
# c_part_max: factor on the memory reserved for the total number of particles.
# mass: particle mass in units of the electron mass.
# charge: particle charge in units of e (−e is the electron charge).
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
        species_geometry = "trapezoidal",
        vacuum_length   = [0.0,  0.0], 
        dens_length_x = [6000, 20],  
        dens_length_y = [2000, 0], 
	n_part_per_cell = 10, 
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
        species_geometry = "trapezoidal",
        vacuum_length   = [0.0,  0.0], 
        dens_length_x = [6000, 20], 
        dens_length_y = [2000, 0], 
	n_part_per_cell = 10,
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
# LASER PROPERTIES
# ----------------
#
# for each laser define:
# a0: maximum amplitude of the laser electric field (in units of the normalization field).
# angle: angle (in degree) at which the laser enters the simulation box.
# delta: polarization parameter, (0:y) (1:z) (0.707106781:circ).
# time_profile: string defining the time profile.
# double_params: vector of real parameters used by the different time-profiles.
#
Laser(
	a0=2.0,
        boxSide = "west",
	angle=0,
	delta=0.0,              
	time_profile =  "gaussian",
        transv_profile = "gaussian",
	double_params = [7.0, 0.0, 0.0], # gaussian:FWHM in intensity, plateau length, delay (in c/w0 and w0)
        double_params_transv = [60.0, 30.8], # distance from axis, transverse FWHM in intensity (in c/w0)
)

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

# print_every (on screen text output) 
print_every = 100

# every for field dump
fieldsToDump = ['Ex','Ey','Bx','Bz','Rho_electron','Rho_proton','Rho','Jx_electron']
fieldDump_every = 100

# every for averagefield Dump
avgfieldDump_every = 0

#Write a restart file every dump_step iterations. 0 for no dump.
dump_step = 0
#Write a restart file every dump_minutes minutes. 0.0 for no dump.
dump_minutes = 0.0

# DIAG ON SCALARS
# every = number of time-steps between each output
#
DiagScalar(every = 100)

exit_after_dump = False
