# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------


# DESCRIPTION OF THE SIMULATION
# Two different species (initially neutral) are considered:
# 1) a carbon slab (from Lx/2-1/2 to Lx/2+1/2) for which 2 constant & homogenous ionization rates are given (Z=0->Z=1,Z=1->Z=2)
# 2) a hydrogen slab (over the full box) for which a space dependent [f(x)] ionization rate is given (Z=0->Z=1)
# Lx is the simulation length and all units are plasma related (c/omega_pe0) with omega_pe0 the plasma frequency
# considering all atoms ionized once (only)

import math
import numpy as np

rest = 64	   # nb of timestep in 1/omega_pe0
resx = 32	   # nb cells in c/omega_pe0
Lx   = 10.     # simulation length in c/omega_pe0
Tsim = 20.	   # duration of the simulation 1/omega_pe0

nppc = 500     # number of particles per cells

def carbon_density_profile(x):
    if (Lx/2.-0.5<x) and (x<Lx/2+0.5):
        return 1.
    else:
        return 0.

carbon_rate_charge0 = 0.1   # ionization rate from Z=0 to Z=1
carbon_rate_charge1 = 0.05  # ionization rate from Z=1 to Z=2
def carbon_ionization_rate(particles):
    rate = np.empty_like(particles.x)
    rate[particles.charge==0] = carbon_rate_charge0
    rate[particles.charge==1] = carbon_rate_charge1
    return rate

hydrogen_max_rate = 0.02
def f(x):
    return np.exp(-(x-Lx/2.)**2/2.)

def hydrogen_ionization_rate(particles):
    x = particles.x
    return np.full_like(particles.x,hydrogen_max_rate) * f(x)

Main(
	geometry = "1Dcartesian",

	interpolation_order = 2,

	cell_length = [1./resx],
	grid_length  = [Lx],

	number_of_patches = [ 16 ],

	timestep = 1./rest,
	simulation_time = Tsim,

	EM_boundary_conditions = [ ['silver-muller'] ],

	reference_angular_frequency_SI = 6*math.pi*1e14,

	random_seed = smilei_mpi_rank
)

Species(
	name = 'carbon',
	ionization_model = 'from_rate',
	ionization_electrons = 'electronC',
	ionization_rate = carbon_ionization_rate,
	maximum_charge_state = 2,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = nppc,
	mass = 1836.0*14.,
	charge = 0.0,
	number_density = carbon_density_profile,
	boundary_conditions = [["periodic"]]
)

Species(
	name = 'electronC',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	boundary_conditions = [["periodic"]]
)

Species(
	name = 'hydrogen',
	ionization_model = 'from_rate',
	ionization_electrons = 'electronH',
	ionization_rate = hydrogen_ionization_rate,
	maximum_charge_state = 1,
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = nppc,
	mass = 1836.0,
	charge = 0.0,
	number_density = 1.,
	boundary_conditions = [["periodic"]]
)

Species(
	name = 'electronH',
	position_initialization = 'regular',
	momentum_initialization = 'cold',
	particles_per_cell = 0,
	mass = 1.0,
	charge = -1.0,
	charge_density = 0.0,
	boundary_conditions = [["periodic"]]
)

### DIAGNOSTICS

DiagScalar(every=1)
DiagFields(every=1)

DiagParticleBinning(
	deposited_quantity = "weight",
	every = 1,
	species = ["carbon"],
	axes = [
		["charge",  -0.5, 2.5, 3]
	]
)
