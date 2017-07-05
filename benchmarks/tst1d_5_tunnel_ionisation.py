# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
l0 = 2.0*math.pi	# wavelength in normalized units
t0 = l0				# optical cycle in normalized units
rest = 6000.0		# nb of timestep in 1 optical cycle
resx = 4000.0		# nb cells in 1 wavelength
Lsim = 0.01*l0	    # simulation length
Tsim = 0.2*t0		# duration of the simulation


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
    
    random_seed = smilei_mpi_rank
)

Species(
	species_type = 'hydrogen',
	ionization_model = 'tunnel',
	ionization_electrons = 'electron',
	atomic_number = 1,
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 10,
	mass = 1836.0*1000.,
	charge = 0.0,
	nb_density = 0.1,
	bc_part_type_xmin = 'none',
	bc_part_type_xmax = 'none'
)

Species(
	species_type = 'carbon',
	ionization_model = 'tunnel',
	ionization_electrons = 'electron',
	atomic_number = 6,
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 10,
	mass = 1836.0*1000.,
	charge = 0.0,
	nb_density = 0.1,
	bc_part_type_xmin = 'none',
	bc_part_type_xmax = 'none'
)

Species(
	species_type = 'electron',
	initPosition_type = 'regular',
	initMomentum_type = 'cold',
	n_part_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	bc_part_type_xmin = 'none',
	bc_part_type_xmax = 'none',
	track_every = 30
)

LaserPlanar1D(
	boxSide = 'xmin',
	a0 = 0.1,
    omega = 1.,
    polarizationPhi = 0.,
    time_envelope = tconstant(),
)

DiagScalar(every = 20)

DiagFields(
    every = 20,
    time_average = 1,
    fields = ["Ex", "Ey", "Ez"]
)

DiagParticles(
	output = "density",
	every = 20,
	species = ["hydrogen"],
	axes = [
		["charge",  -0.5, 1.5, 2]
	]
)

DiagParticles(
	output = "density",
	every = 20,
	species = ["carbon"],
	axes = [
		["charge",  -0.5, 6.5, 7]
	]
)
