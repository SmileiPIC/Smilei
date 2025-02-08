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
	geometry = "1Dcartesian",
	 
	interpolation_order = 2,
	 
	cell_length = [l0/resx],
	grid_length  = [Lsim],
	
	number_of_patches = [ 4 ],
	
	timestep = t0/rest,
	simulation_time = Tsim,
	 
	EM_boundary_conditions = [ ['silver-muller'] ],
	
	reference_angular_frequency_SI = 6*math.pi*1e14,
	
)

# Tunnel
Species(
	name = 'hydrogen_tunnel',
	ionization_model = 'tunnel',
	ionization_electrons = 'electron_tunnel',
	atomic_number = 1,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 40,
	mass = 1836.0*1000.,
	charge = 0.0,
	number_density = 0.1,
	boundary_conditions = [
		["remove", "remove"],
	],
)

Species(
	name = 'carbon_tunnel',
	ionization_model = 'tunnel',
	ionization_electrons = 'electron_tunnel',
	atomic_number = 6,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 40,
	mass = 1836.0*1000.,
	charge = 0.0,
	number_density = 0.1,
	boundary_conditions = [
		["remove", "remove"],
	],
)

Species(
	name = 'electron_tunnel',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	boundary_conditions = [
		["remove", "remove"],
	],
	keep_interpolated_fields = ["Ex", "Ey", "Ez", "Wx", "Wy", "Wz"],
)

# tunnel_full_PPT
Species(
	name = 'hydrogen_tunnel_full_PPT',
	ionization_model = 'tunnel_full_PPT',
	ionization_electrons = 'electron_tunnel_full_PPT',
	atomic_number = 1,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 40,
	mass = 1836.0*1000.,
	charge = 0.0,
	number_density = 0.1,
	boundary_conditions = [
		["remove", "remove"],
	],
)

Species(
	name = 'carbon_tunnel_full_PPT',
	ionization_model = 'tunnel_full_PPT',
	ionization_electrons = 'electron_tunnel_full_PPT',
	atomic_number = 6,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 40,
	mass = 1836.0*1000.,
	charge = 0.0,
	number_density = 0.1,
	boundary_conditions = [
		["remove", "remove"],
	],
)

Species(
	name = 'electron_tunnel_full_PPT',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	boundary_conditions = [
		["remove", "remove"],
	],
	keep_interpolated_fields = ["Ex", "Ey", "Ez", "Wx", "Wy", "Wz"],
)

# Tong & Lin
Species(
	name = 'hydrogen_tunnel_TL',
	ionization_model = 'tunnel_TL',
	ionization_electrons = 'electron_tunnel_TL',
	atomic_number = 1,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 40,
	mass = 1836.0*1000.,
	charge = 0.0,
	number_density = 0.1,
	boundary_conditions = [
		["remove", "remove"],
	],
)

Species(
	name = 'carbon_tunnel_TL',
	ionization_model = 'tunnel_TL',
	ionization_electrons = 'electron_tunnel_TL',
	atomic_number = 6,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 40,
	mass = 1836.0*1000.,
	charge = 0.0,
	number_density = 0.1,
	boundary_conditions = [
		["remove", "remove"],
	],
)

Species(
	name = 'electron_tunnel_TL',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	boundary_conditions = [
		["remove", "remove"],
	],
	keep_interpolated_fields = ["Ex", "Ey", "Ez", "Wx", "Wy", "Wz"],
)

# BSI
Species(
	name = 'hydrogen_tunnel_BSI',
	ionization_model = 'tunnel_BSI',
	ionization_electrons = 'electron_tunnel_BSI',
	atomic_number = 1,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 40,
	mass = 1836.0*1000.,
	charge = 0.0,
	number_density = 0.1,
	boundary_conditions = [
		["remove", "remove"],
	],
)

Species(
	name = 'carbon_tunnel_BSI',
	ionization_model = 'tunnel_BSI',
	ionization_electrons = 'electron_tunnel_BSI',
	atomic_number = 6,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 40,
	mass = 1836.0*1000.,
	charge = 0.0,
	number_density = 0.1,
	boundary_conditions = [
		["remove", "remove"],
	],
)

Species(
	name = 'electron_tunnel_BSI',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	boundary_conditions = [
		["remove", "remove"],
	],
	keep_interpolated_fields = ["Ex", "Ey", "Ez", "Wx", "Wy", "Wz"],
)

def By(t):
	return 1e-7 * math.sin(t)
def Bz(t):
	return 10 * math.sin(t)

Laser(
	box_side = "xmin",
	space_time_profile = [By, Bz],
)

DiagScalar(every = 20)

DiagFields(
	every = 20,
	time_average = 1,
	fields = ["Ex", "Ey", "Ez"]
)

DiagParticleBinning(
    name = "hydrogen_tunnel",
	deposited_quantity = "weight",
	every = 20,
	species = ["hydrogen_tunnel"],
	axes = [
		["charge",  -0.5, 1.5, 2]
	]
)

DiagParticleBinning(
    name = "hydrogen_tunnel_full_PPT",
	deposited_quantity = "weight",
	every = 20,
	species = ["hydrogen_tunnel_full_PPT"],
	axes = [
		["charge",  -0.5, 1.5, 2]
	]
)

DiagParticleBinning(
    name = "hydrogen_tunnel_TL",
	deposited_quantity = "weight",
	every = 20,
	species = ["hydrogen_tunnel_TL"],
	axes = [
		["charge",  -0.5, 1.5, 2]
	]
)

DiagParticleBinning(
    name = "hydrogen_tunnel_BSI",
	deposited_quantity = "weight",
	every = 20,
	species = ["hydrogen_tunnel_BSI"],
	axes = [
		["charge",  -0.5, 1.5, 2]
	]
)

DiagParticleBinning(
    name = "carbon_tunnel",
	deposited_quantity = "weight",
	every = 20,
	species = ["carbon_tunnel"],
	axes = [
		["charge",  -0.5, 6.5, 7]
	]
)

DiagParticleBinning(
    name = "carbon_tunnel_full_PPT",
	deposited_quantity = "weight",
	every = 20,
	species = ["carbon_tunnel_full_PPT"],
	axes = [
		["charge",  -0.5, 6.5, 7]
	]
)

DiagParticleBinning(
    name = "carbon_tunnel_TL",
	deposited_quantity = "weight",
	every = 20,
	species = ["carbon_tunnel_TL"],
	axes = [
		["charge",  -0.5, 6.5, 7]
	]
)

DiagParticleBinning(
    name = "carbon_tunnel_BSI",
	deposited_quantity = "weight",
	every = 20,
	species = ["carbon_tunnel_BSI"],
	axes = [
		["charge",  -0.5, 6.5, 7]
	]
)

DiagTrackParticles(
	species = "electron_tunnel",
	every = [1,1000,30],
	attributes = ["x","px","py","pz","w","Wy"]
)

DiagTrackParticles(
	species = "electron_tunnel_full_PPT",
	every = [1,1000,30],
	attributes = ["x","px","py","pz","w","Wy"]
)

DiagTrackParticles(
	species = "electron_tunnel_TL",
	every = [1,1000,30],
	attributes = ["x","px","py","pz","w","Wy"]
)

DiagTrackParticles(
	species = "electron_tunnel_BSI",
	every = [1,1000,30],
	attributes = ["x","px","py","pz","w","Wy"]
)

DiagNewParticles(
	species = "electron_tunnel",
	every = 100,
	attributes = ["x","py","w","q"],
)

DiagNewParticles(
	species = "electron_tunnel_full_PPT",
	every = 100,
	attributes = ["x","py","w","q"],
)

DiagNewParticles(
	species = "electron_tunnel_TL",
	every = 100,
	attributes = ["x","py","w","q"],
)

DiagNewParticles(
	species = "electron_tunnel_BSI",
	every = 100,
	attributes = ["x","py","w","q"],
)
