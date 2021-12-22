# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.*math.pi			# laser wavelength
t0 = l0					# optical cycle
Lsim = [20.*l0,28.*l0]	# length of the simulation
Tsim = 40.*t0			# duration of the simulation
resx = 20.				# nb of cells in on laser wavelength
rest = 30.				# time of timestep in one optical cycle


Main(
	geometry = "2Dcartesian",

	interpolation_order = 2 ,

	cell_length = [l0/resx,l0/resx],
	grid_length  = Lsim,

	number_of_patches = [ 16, 20 ],
	patch_arrangement = "linearized_XY",

	timestep = t0/rest,
	simulation_time = Tsim,

	EM_boundary_conditions = [
		['silver-muller'],
		['reflective'],
	],

)

Vectorization(
    mode = "on",
)

LaserGaussian2D(
	box_side         = "xmin",
	a0              = 0.4,
	focus           = [10.*l0, 25.0*l0],
	waist           = 5.0*l0,
	incidence_angle = 20./180.*math.pi,
	time_envelope   = tgaussian(fwhm=5.*t0, center=10.*t0)
)

Species(
	name = 'ion',
	position_initialization = 'random',
	momentum_initialization = 'cold',
	ionization_model = 'none',
	particles_per_cell = 40,
	c_part_max = 1.0,
	mass = 1836.0,
	charge = 1.0,
	number_density = trapezoidal(2.0,xvacuum=11.*l0,xplateau=2.*l0),
	time_frozen = Tsim,
	boundary_conditions = [
		["reflective", "reflective"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'eon',
	position_initialization = 'random',
	momentum_initialization = 'cold',
	ionization_model = 'none',
	particles_per_cell = 40,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = trapezoidal(2.0,xvacuum=11.*l0,xplateau=2.*l0),
	time_frozen = 0.,
	boundary_conditions = [
		["reflective", "reflective"],
		["periodic", "periodic"],
	],
)


DiagScalar(every=5)


DiagFields(
	every = 15,
	fields = ['Ex','Ey','Rho_eon']
)

DiagPerformances(
    every = 30,
)
