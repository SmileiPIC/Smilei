# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
l0 = 2.0*math.pi	# wavelength in normalized units
t0 = l0				# optical cycle in normalized units
resx = 64.		    # nb cells in 1 wavelength
rest = resx/0.95    # temporal resolution
Lv   = 2.*l0
Lp   = l0/32.
Lx   = 2.*l0+Lp      # simulation length in x direction
Ly   = 2.*l0+Lp      # simulation length in y direction

tmax = 20.*t0
Tsim = Lv+Lp+tmax  # duration of the simulation
nppc = 2048*4        # nb of particle per cells

n0   = 1.e-2         # neutral density

Lmu  = 0.8
I18  = 5.0e-2
aL   = math.sqrt(I18*Lmu**2/1.38)

def n_(x):
	if (Lv<x<Lv+Lp):
		return n0
	else:
		return 0.

Main(
    geometry = "1Dcartesian",
     
    interpolation_order = 2,
     
    cell_length = [l0/resx],
    grid_length  = [Lx],
    
    number_of_patches = [ 1 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [['silver-muller']],
    
    reference_angular_frequency_SI = 2.*math.pi * 3.e8/(Lmu*1.e-6),
    
    random_seed = 0
)

Species(
	name = 'carbon',
	ionization_model = 'tunnel',
	ionization_electrons = 'electron',
	atomic_number = 6,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = nppc,
	mass = 12.*1836.0,
	charge = 0.0,
	number_density = n_,
	boundary_conditions = [['reflective']]
)

Species(
	name = 'electron',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	boundary_conditions = [['reflective']],
#   time_frozen = 2.*Tsim
)

LaserPlanar1D(
    box_side         = "xmin",
    a0              = aL,
    omega           = 1.,
    polarization_phi = 0.,
    ellipticity     = 0.,
    time_envelope   = tgaussian(start=0., duration=None, fwhm=5.*t0, center=Lv+tmax/2., order=2)
)
 
DiagScalar(every = 1)

DiagFields(
    every = 1,
    time_average = 1,
    fields = ["Ex", "Ey", "Ez"]
)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = 1,
	species = ["carbon"],
	axes = [
		["charge",  -0.5, 6.5, 7]
	]
)


