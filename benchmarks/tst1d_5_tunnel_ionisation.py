# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
l0 = 2.0*math.pi	# wavelength in normalized units
t0 = l0				# optical cycle in normalized units
rest = 200.0		# nb of timestep in 1 optical cycle
resx = 100.0		# nb cells in 1 wavelength
Lsim = l0		# simulation length
Tsim = 5.0*t0		# duration of the simulation


Main(
    geometry = "1d3v",
     
    interpolation_order = 2,
     
    cell_length = [l0/resx],
    sim_length  = [Lsim],
    
    number_of_patches = [ 4 ],
    
    timestep = t0/rest,
    sim_time = Tsim,
     
    bc_em_type_x = ['silver-muller'],
    
    referenceAngularFrequency_SI = 6*math.pi*1e14,
    
    random_seed = 0
)

Species(
	species_type = 'helium',
	ionization_model = 'tunnel',
	ionization_electrons = 'electron',
	atomic_number = 2,
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 100,
	mass = 1836.0,
	charge = 0.0,
	nb_density = trapezoidal(1.0,xvacuum=0.49*l0,xplateau=0.02*l0),
	bc_part_type_west = 'none',
	bc_part_type_east = 'none'
)
Species(
	species_type = 'electron',
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	bc_part_type_west = 'none',
	bc_part_type_east = 'none'
)

LaserPlanar1D(
	boxSide = 'west',
	a0 = 0.1,
    omega = 1.,
    polarizationPhi = math.pi/2.,
    time_envelope = tconstant(),
)
 
DiagScalar(every = 10)
 
DiagFields(
    every = 1
)

DiagParticles(
	output = "density",
	every = 10,
	species = ["electron"],
	axes = [
		["x",  0.45*Lsim, 0.55*Lsim, 200],
		["px", -0.1, 0.1, 200]
	]
)

