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

timestep = 1 * L0
sim_time  = 400 * L0


time_fields_frozen = 100000000000.

# SIMULATION BOX : for all space directions (in 2D & 3D use vector of doubles)
cell_length = [20.*L0]
sim_length  = [160.*L0]

bc_em_type_x  = ["periodic"]


random_seed = 0


Z = 30
A = 65.
density = 10.

electrons = []
ions = []
temperature = []
Tmin = 0.001 # keV
Tmax = 2.
npoints = 10

for i in range(npoints):
	eon = "electron"+str(i)
	ion = "ion"+str(i)
	electrons.append(eon)
	ions.append(ion)
	
	T = math.exp(math.log(Tmin) + float(i)/(npoints-1)*math.log(Tmax/Tmin)) #logscale
	T /= 511.
	temperature.append(T)
	
	Zstar = (Z)*(1.-math.exp(-T/(1./511.)))
	if Zstar<0.01: Zstar=0.01
	
	Species(
		species_type = eon,
		initPosition_type = "regular",
		initMomentum_type = "maxwell-juettner",
		n_part_per_cell= 100,
		mass = 1.0,
		charge = -1.0,
		charge_density = Zstar*density,
		mean_velocity = [0., 0., 0.],
		temperature = [T]*3,
		time_frozen = 100000000.0,
		bc_part_type_west = "none",
		bc_part_type_east = "none",
		bc_part_type_south = "none",
		bc_part_type_north = "none",
		c_part_max = 10.
	)
	
	Species(
		species_type = ion,
		initPosition_type = "regular",
		initMomentum_type = "maxwell-juettner",
		n_part_per_cell= 100,
		mass = 1836.0*A,
		charge = Zstar,
		nb_density = 10.,
		mean_velocity = [0., 0., 0.],
		temperature = [T]*3,
		time_frozen = 100000000.0,
		bc_part_type_west = "none",
		bc_part_type_east = "none",
		bc_part_type_south = "none",
		bc_part_type_north = "none",
		atomic_number = Z
	)
	
	Collisions(
		species1 = [eon],
		species2 = [ion],
		coulomb_log = 0.,
		ionizing = True
	)
	
	DiagParticles(
		output = "ekin_density",
		every = 10,
		species = [eon],
		axes = [ ["x", 0, sim_length[0], 1] ]
	)
	DiagParticles(
		output = "density",
		every = 10,
		species = [eon],
		axes = [ ["x", 0, sim_length[0], 1] ]
	)
	DiagParticles(
		output = "charge_density",
		every = 10,
		species = [ion],
		axes = [ ["x", 0, sim_length[0], 1] ]
	)
	DiagParticles(
		output = "density",
		every = 10,
		species = [ion],
		axes = [ ["x", 0, sim_length[0], 1] ]
	)
	#DiagParticles(
	#	output = "density",
	#	every = 50,
	#	species = [ion],
	#	axes = [ ["charge", -0.5, Z+0.5, Z+1] ]
	#)


# print_every (on screen text output) 
print_every = 100

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



