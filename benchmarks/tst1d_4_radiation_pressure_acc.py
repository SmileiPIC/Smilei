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
    geometry = "1d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx],
    sim_length  = [Lsim],
    
    number_of_patches = [ 8 ],
    
    timestep = t0/rest,
    sim_time = Tsim,
     
    bc_em_type_x = ['silver-muller'],
     
    random_seed = 0
)

Species(
	species_type = 'ion',
	initPosition_type = 'random',
	initMomentum_type = 'cold',
	n_part_per_cell = 10,
	mass = 1836.0,
	charge = 1.0,
	nb_density = trapezoidal(10.,xvacuum=l0,xplateau=l0),
	temperature = [0.],
	bc_part_type_xmin = 'refl',
	bc_part_type_xmax = 'refl'
)
Species(
	species_type = 'eon',
	initPosition_type = 'random',
	initMomentum_type = 'cold',
	n_part_per_cell = 10,
	mass = 1.0,
	charge = -1.0,
	nb_density = trapezoidal(10.,xvacuum=l0,xplateau=l0),
	temperature = [0.],
	bc_part_type_xmin = 'refl',
	bc_part_type_xmax = 'refl'
)

LaserPlanar1D(
	boxSide = 'xmin',
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

DiagParticles(
	output = "density",
	every = every,
	species = ["ion"],
	axes = [
		["x",  0.,   Lsim, 200],
		["px", -10., 1000., 200]
	]
)

DiagParticles(
	output = "density",
	every = every,
	species = ["ion"],
	axes = [
		["ekin", 0., 200., 200, "edge_inclusive"]
	]
)


