# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi	# laser wavelength
t0 = l0				# optical cicle
Lsim = 10.*l0		# length of the simulation
Tsim = 40.*t0		# duration of the simulation
resx = 500.			# nb of cells in on laser wavelength
rest = resx/0.95	# time of timestep in one optical cycle (0.95 * CFL)

# plasma slab
def f(x):
    if l0 < x < 2.0*l0:
        return 1.0
    else :
        return 0.0

Main(
    geometry = "1Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx],
    grid_length  = [Lsim],
    
    number_of_patches = [ 8 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [ ['silver-muller'] ],
     
    cluster_width = 1
)

Species(
	name = 'ion',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 10,
	mass = 1836.0,
	charge = 1.0,
	number_density = trapezoidal(10.,xvacuum=l0,xplateau=l0),
	temperature = [0.],
	boundary_conditions = [
		["reflective", "reflective"],
	],
)
Species(
	name = 'eon',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	number_density = trapezoidal(10.,xvacuum=l0,xplateau=l0),
	temperature = [0.],
	boundary_conditions = [
		["reflective", "reflective"],
	],
)

LaserPlanar1D(
	box_side = 'xmin',
	a0 = 10.,
    omega = 1.,
    ellipticity = 1.,
    time_envelope = tconstant(),
)


every = int(rest/2.)

DiagFields(
    every = every,
    fields = ['Ex','Ey','Ez','Rho_ion','Rho_eon']
)

DiagScalar(every=every)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = every,
	species = ["ion"],
	axes = [
		["x",  0.,   Lsim, 200],
		["px", -10., 1000., 200]
	]
)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = every,
	species = ["ion"],
	axes = [
		["ekin", 0., 200., 200, "edge_inclusive"]
	]
)


for direction in ["forward", "backward", "both", "canceling"]:
	DiagScreen(
	    shape = "sphere",
	    point = [0.],
	    vector = [Lsim/3.],
	    direction = direction,
	    deposited_quantity = "weight",
	    species = ["eon"],
	    axes = [["ekin", 0., 0.4, 10],],
	    every = 3000
	)
	DiagScreen(
	    shape = "plane",
	    point = [Lsim/3.],
	    vector = [1.],
	    direction = direction,
	    deposited_quantity = "weight",
	    species = ["eon"],
	    axes = [["ekin", 0., 0.4, 10],],
	    every = 3000
	)

