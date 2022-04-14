# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticleBinning,
#           DiagScalar, DiagPhase or ExternalField

import math as m


v0   = 0.
g0   = 1.
Te   = 1.
vth  = m.sqrt(Te)
Lde  = m.sqrt(Te)
Lx   = 64.*Lde
t_sim = Lx/1.

pmax = 5.*m.sqrt(g0*Te)

nppc = 1000
dx   = 0.5*Lde
dt   = 0.95*dx/m.sqrt(1.)  


# MAIN SMILEI INPUT
Main(
     geometry = '1Dcartesian',
     interpolation_order = 2,
     cell_length = [dx],
     grid_length  = [Lx],
     number_of_patches=[1],
     cluster_width=1,
     timestep = dt,
     simulation_time = 0.,
     EM_boundary_conditions = [ ['periodic'] ],
     print_every = int(t_sim/dt/100.)
)

# SPECIES DEFINITION
Species(
	name = 'ion',
	position_initialization = 'random',
	momentum_initialization = 'mj',
	temperature = [Te,Te,Te],
	particles_per_cell = nppc,
	mass = 1.0,
	charge = 1.0,
	number_density = 1.,
	mean_velocity=[0.,0.,v0],
	boundary_conditions = [["periodic"]],
)
Species(
	name = 'eon',
	position_initialization = 'random',
	momentum_initialization = 'mj',
	temperature = [Te,Te,Te],
	particles_per_cell = nppc,
	mass = 1.0,
	charge = -1.0,
	number_density = 1.,
	mean_velocity=[0.,0.,v0],
	boundary_conditions = [["periodic"]],
)



# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

globalEvery = int(t_sim/dt/50.) 

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery ,
    fields = ['Ex','Ey','Ez','By_m','Bz_m','Jz','Rho_eon1','Rho_eon2']
)

DiagParticleBinning(
 	deposited_quantity = "weight",
 	every = globalEvery,
 	time_average = 1,
 	species = ["eon"],
 	axes = [
 		["gamma", 1., 1+15.*Te, 1000]
 	]
)
DiagParticleBinning(
 	deposited_quantity = "weight",
 	every = globalEvery,
 	time_average = 1,
 	species = ["eon"],
 	axes = [
 		["px", -10., 10., 500]
 	]
)
