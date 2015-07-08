# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField


# sim_units: normalisation units for the input data
#            it is used only in the input data & log file
#            codes outputs are always in "normalised" units
#            'wavelength' : input data are in wavelength-related units
#            'normalized' : input data are put in code (relativistic) units
sim_units = "wavelength"
wavelength_SI = 1.e-6

# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
dim = "1d3v"

# order of interpolation
interpolation_order = 2

# SIMULATION TIME 
# set either the resolution (res_time) or the timestep
# res_time = integer, number of time-steps within one unit of time (`sim_units`)
# timestep = float, time step in units of `sim_units`
# sim_time = float, duration of the simulation  in units of `sim_units`
timestep = 0.002
sim_time  = 0.1


#  optional parameter time_fields_frozen, during which fields are not updated
time_fields_frozen = 100000000000.

# SIMULATION BOX : for all space directions (in 2D & 3D use vector of doubles)
# either use the resolution (res_space) or cell-length (cell_length)
# res_space   = list of integers, number of cells in one unit of space (`sim_units`)
# sim_length  = length of the simulation in units of `sim_units`
# cell_length = cell length  in units of `sim_units`
cell_length = [1]
sim_length  = [40]

# ELECTROMAGNETIC BOUNDARY CONDITIONS
# bc_em_type_long/trans : boundary conditions used for EM fields 
#                         in the longitudinal or transverse directions
#                         'periodic'      : periodic BC (using MPI topology)
#                         'silver-muller' : injecting/absorbing
bc_em_type_long  = "periodic"


# RANDOM seed used to randomize the random number generator
random_seed = 0

# DEFINE ALL SPECIES
# species_type       = string, given name to the species (e.g. ion, electron, positron, test ...)
# initPosition_type  = string, "regular" or "random"
# initMomentum_type  = string "cold", "maxwell-juettner" or "rectangular"
# n_part_per_cell    = integer, number of particles/cell
# c_part_max         = float, factor on the memory reserved for the total number of particles
# mass               = float, particle mass in units of the electron mass
# charge             = float, particle charge in units of the electron charge
# dynamics_type      = string, type of species dynamics = "norm" or "rrLL"
# time_frozen        = float, time during which particles are frozen in units of the normalization time
# radiating          = boolean, if true, incoherent radiation calculated using the Larmor formula 
# vacuum_length      = list of floats, distance from box borders without particles.
# charge_density     = float, species charge density in units of the "critical" density
#     or nb_density for number density
# mean_velocity      = list of floats, mean velocity in units of the speed of light
# temperature        = list of floats, temperature in units of m_e c^2
# SPECIES PROFILES from python function (see doc)
#    Predefined functions: constant, trapezoidal, gaussian, polygonal, cosine
# dens_profile       = python function. Units: n_c
# mvel_[xyz]_profile = python function. Units: c
# temp_[xyz]_profile = python function. Units: m_e c^2
Species(
	species_type = "ghostChargeLayerLeft",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 10,
	mass = 1000.,
	charge =  1.,
	charge_density = 0.00064,
	dens_profile = trapezoidal(xvacuum=5., xplateau=1.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.0001],
	time_frozen = 100000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)
Species(
	species_type = "ghostChargeLayerRight",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 10,
	mass = 1000.,
	charge = -1.,
	charge_density = 0.00064,
	dens_profile = trapezoidal(xvacuum=35., xplateau=1.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.0001],
	time_frozen = 100000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "copper1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 10000,
	mass = 115845.,      # =  mass of Cu atom
	charge = 5.6,
	charge_density = 415.,   # =  density of solid Cu
	dens_profile = trapezoidal(xvacuum=15., xplateau=10.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.00004], # 20 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)
Species(
	species_type = "electron1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 50000,
	mass = 1.0,
	charge = -1.0,
	charge_density = 415.,
	dens_profile = trapezoidal(xvacuum=15., xplateau=10.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.00004], # 20 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "copper2",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 10000,
	mass = 115845.,      # =  mass of Cu atom
	charge = 7.4,
	charge_density = 554.,   # =  density of solid Cu
	dens_profile = trapezoidal(xvacuum=15., xplateau=10.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.0001], # 50 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)
Species(
	species_type = "electron2",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 100000,
	mass = 1.0,
	charge = -1.0,
	charge_density = 554.,
	dens_profile = trapezoidal(xvacuum=15., xplateau=10.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.0001], # 50 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "copper3",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell = 10000,
	mass = 115845.,      # =  mass of Cu atom
	charge = 10.,
	charge_density = 757.,   # =  density of solid Cu
	dens_profile = trapezoidal(xvacuum=15., xplateau=10.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.0002], # 100 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)
Species(
	species_type = "electron3",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 50000,
	mass = 1.0,
	charge = -1.0,
	charge_density = 757.,
	dens_profile = trapezoidal(xvacuum=15., xplateau=10.),
	mean_velocity = [0., 0., 0.],
	temperature = [0.0002], # 100 eV
	time_frozen = 0.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)


# COLLISIONS
# species1    = list of strings, the names of the first species that collide
# species2    = list of strings, the names of the second species that collide
#               (can be the same as species1)
# coulomb_log = float, Coulomb logarithm. If negative or zero, then automatically computed.
Collisions(
	species1 = ["copper1"],
	species2 = ["electron1"],
	coulomb_log = 2.
)
Collisions(
	species1 = ["copper2"],
	species2 = ["electron2"],
	coulomb_log = 2.
)
Collisions(
	species1 = ["copper3"],
	species2 = ["electron3"],
	coulomb_log = 2.
)

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

# print_every (on screen text output) 
print_every = 10

# DIAGNOSTICS ON FIELDS
fieldDump_every    = 5
avgfieldDump_every = 5
ntime_step_avg     = 1


# DIAGNOSTICS ON SCALARS
# every = integer, number of time-steps between each output
DiagScalar(
	every = 1
)


# DIAGNOSTICS ON PARTICLES - project the particles on a N-D arbitrary grid
# ------------------------------------------------------------------------
# output       = string: "density", "charge_density" or "current_density_[xyz]"
#                parameter that describes what quantity is obtained 
# every        = integer > 0: number of time-steps between each output
# time_average = integer > 0: number of time-steps to average
# species      = list of strings, one or several species whose data will be used
# axes         = list of axes
# Each axis is a list: (_type_ _min_ _max_ _nsteps_ ["logscale"] ["edge_inclusive"])
#   _type_ is a string, one of the following options:
#      x, y, z, px, py, pz, p, gamma, ekin, vx, vy, vz, v or charge
#   The data is discretized for _type_ between _min_ and _max_, in _nsteps_ bins
#   The optional "logscale" sets the scale to logarithmic
#   The optional "edge_inclusive" forces the particles that are outside (_min_,_max_)
#     to be counted in the extrema bins
#   Example : axes = ("x", 0, 1, 30)
#   Example : axes = ("px", -1, 1, 100, "edge_inclusive")

DiagParticles(
	output = "density",
	every = 5,
	time_average = 4,
	species = ["electron1"],
	axes = [
		 ["vx",  -0.03,  0.03,    1000]
	]
)
DiagParticles(
	output = "density",
	every = 5,
	time_average = 4,
	species = ["electron2"],
	axes = [
		 ["vx",  -0.04,  0.04,    1000]
	]
)
DiagParticles(
	output = "density",
	every = 5,
	time_average = 4,
	species = ["electron3"],
	axes = [
		 ["vx",  -0.06,  0.06,    1000]
	]
)

