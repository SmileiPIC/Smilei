# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import numpy as np
l0 = 2.0*np.pi	# wavelength in normalized units
t0 = l0				# optical cycle in normalized units

# incoming electromagnetic plane-wave
a0   = 0.1
T    = 4.*t0
U0   = 0.5*a0**2*T

def f(t):
	if (t>0) and (t<2.*T):
		return np.sin(np.pi/2. * t/T)
	else:
		return 0.

def By(t):
	return 0.

def Bz(t):
	if (t>t0):
		return a0 * np.sin(t-t0) * f(t-t0)
	else:
		return 0.		

# hydrogen plasma
D    = 20.*l0
n0   = 0.1*U0 * (1./D) * (5.11e5/13.6)  # chooses slab density so that 20% of the pulse is absorbed during ionization

def nh_(x):
	if (x>2.*T+l0) and (x<D+2.*T+l0):
		return n0
	else:
		return 0.

# numerical parameters
nppc      = 1
dx        = l0/64.
dt        = 0.95*dx
diagEvery = int(t0/dt/4.)
Lsim      = [2.*T+l0+D+2.*T]
Tsim      = Lsim[0] + t0


############################################################################################################
Main(
	geometry = "1Dcartesian",
	 
	interpolation_order = 4,
	 
	cell_length = [dx],
	grid_length = Lsim,
	
	number_of_patches = [ 16 ],
	
	timestep = dt,
	simulation_time = Tsim,
	 
	EM_boundary_conditions = [['silver-muller']],
	
	reference_angular_frequency_SI = 2.*np.pi*3.e8/1.e-6,
	
	random_seed = smilei_mpi_rank
)

Species(
	name = 'hydrogen',
	ionization_model = 'tunnel',
	ionization_electrons = 'electron',
	atomic_number = 1,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = nppc,
	mass = 1836.0 * 1.e6,
	charge = 0.0,
	number_density = nh_,
	boundary_conditions = [["reflective"]],
)

Species(
	name = 'electron',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	boundary_conditions = [["reflective"]],
	time_frozen=2.*Tsim
)

Laser(
	box_side = "xmin",
	space_time_profile = [By, Bz],
)

DiagScalar(every = 1)

DiagFields(
	every = diagEvery,
	time_average = 1,
	fields = ["Ey","Bz_m","Rho_hydrogen","Rho_electron"]
)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = diagEvery,
	species = ["hydrogen"],
	axes = [
		["charge",  -0.5, 1.5, 2]
	]
)

DiagTrackParticles(
	species = "electron",
	every = diagEvery
)