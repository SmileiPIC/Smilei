# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength

referenceAngularFrequency_SI = L0 * 3e8 /1.e-6

Main(
    geometry = "1d3v"

# number_of_patches: list of the number of patches in each dimension
number_of_patches = [ 4 ]

interpolation_order = 2

timestep = 1. * L0
sim_time  = 50 * L0


time_fields_frozen = 100000000000.

# SIMULATION BOX : for all space directions (in 2D & 3D use vector of doubles)
cell_length = [2.5*L0]
sim_length  = [60.*L0]

bc_em_type_x  = ["periodic"]


random_seed = 0

Species(
	species_type = "backgroundelectron",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	ionization_model = "none",
	n_part_per_cell = 10000,
	mass = 1.,
	charge = -1.0,
	charge_density = 10.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.01],
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "electron1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 10000,
	mass = 1.0,
	charge = -1.0,
	charge_density = trapezoidal(0.00001, xplateau=20.*L0),
	mean_velocity = [0.2, 0., 0.],
	temperature = [0.0000001],
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "electron2",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 10000,
	mass = 1.0,
	charge = -1.0,
	charge_density = trapezoidal(0.00001, xvacuum=20.*L0, xplateau=20.*L0),
	mean_velocity = [0.27, 0., 0.],
	temperature = [0.0000001],
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

Species(
	species_type = "electron3",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 10000,
	mass = 1.0,
	charge = -1.0,
	charge_density = trapezoidal(0.00001,xvacuum=40.*L0, xplateau=20.*L0),
	mean_velocity = [0.33, 0., 0.],
	temperature = [0.0000001],
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none"
)

# COLLISIONS
# species1    = list of strings, the names of the first species that collide
# species2    = list of strings, the names of the second species that collide
#               (can be the same as species1)
# coulomb_log = float, Coulomb logarithm. If negative or zero, then automatically computed.
Collisions(
	species1 = ["backgroundelectron"],
	species2 = ["electron1", "electron2", "electron3"],
	coulomb_log = 3.
)


# print_every (on screen text output) 
print_every = 10

# DIAGNOSTICS ON FIELDS
DiagFields(
	every = 1
)


# DIAGNOSTICS ON SCALARS
# every = integer, number of time-steps between each output
DiagScalar(
	every = 1
)


# DIAGNOSTICS ON PARTICLES - project the particles on a N-D arbitrary grid
# ------------------------------------------------------------------------
# output       = string: "density", "charge_density" or "jx_density"
#                parameter that describes what quantity is obtained 
# every        = integer > 0: number of time-steps between each output
# time_average = integer > 0: number of time-steps to average
# species      = list of strings, one or several species whose data will be used
# axes         = list of axes
# Each axis is a list: [_type_,_min_,_max_,_nsteps_,"logscale","edge_inclusive"]
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
	every = 2,
	time_average = 1,
	species = ["electron1", "electron2", "electron3"],
	axes = [
		 ["x",    0*L0,    60.*L0,   3],
		 ["ekin",   0.01,  0.1,   500, "logscale"]
	]
)
