# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

import math
L0 = 2.*math.pi # conversion from normalization length to wavelength

referenceAngularFrequency_SI = L0 * 3e8 /1.e-6

Main(
    geometry = "1d3v"

interpolation_order = 2

# number_of_patches: list of the number of patches in each dimension
number_of_patches = [ 4 ]

timestep = 50 * L0
sim_time  = 2000 * L0


time_fields_frozen = 100000000000.

# SIMULATION BOX : for all space directions (in 2D & 3D use vector of doubles)
cell_length = [20.*L0]
sim_length  = [160.*L0]

bc_em_type_x  = ["periodic"]


random_seed = 0



el = "electron1"
E = 1000. # keV
E /= 511.
vel = math.sqrt(1.-1./(1.+E)**2)
mom = math.sqrt((1.+E)**2-1.)
Species(
	species_type = el,
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 100,
	mass = 1.0,
	charge = -1.0,
	nb_density = 1.,
	mean_velocity = [vel, 0., 0.],
	temperature = [0.0000000001]*3,
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none",
	bc_part_type_south = "none",
	bc_part_type_north = "none",
	c_part_max = 10.
)

Species(
	species_type = "ion1",
	initPosition_type = "regular",
	initMomentum_type = "maxwell-juettner",
	n_part_per_cell= 100,
	mass = 1836.0*27.,
	charge = 0,
	nb_density = 1.,
	mean_velocity = [0., 0., 0.],
	temperature = [0.00000000001]*3,
	time_frozen = 100000000.0,
	bc_part_type_west = "none",
	bc_part_type_east = "none",
	bc_part_type_south = "none",
	bc_part_type_north = "none",
	atomic_number = 13)


# COLLISIONS
# species1    = list of strings, the names of the first species that collide
# species2    = list of strings, the names of the second species that collide
#               (can be the same as species1)
# coulomb_log = float, Coulomb logarithm. If negative or zero, then automatically computed.
Collisions(
	species1 = [el],
	species2 = ["ion1"],
	coulomb_log = 0.00000001,
	ionizing = True
)


# print_every (on screen text output) 
print_every = 10

# DIAGNOSTICS ON FIELDS
DiagFields(
	every = 1000000
)


# DIAGNOSTICS ON SCALARS
# every = integer, number of time-steps between each output
DiagScalar(
	every = 10000000
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
	every = 1,
	species = [el],
	axes = [
		 ["ekin", 0., E, 1000]
	]
)
