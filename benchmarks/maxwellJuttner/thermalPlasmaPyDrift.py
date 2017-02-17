# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField

import math as m


g0   = 5.
v0   = m.sqrt(g0**2-1.)/g0
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
     geometry = '1d3v',
     interpolation_order = 2,
     cell_length = [dx],
     sim_length  = [Lx],
     number_of_patches=[1],
     clrw=1,
     timestep = dt,
     sim_time = 0.,
     bc_em_type_x = ['periodic'],
     random_seed = smilei_mpi_rank,
     print_every = int(t_sim/dt/100.)
)

# SPECIES DEFINITION
Species(
	species_type = 'ion',
	initPosition_type = 'random',
	initMomentum_type = 'mj',
	temperature = [Te,Te,Te],
	n_part_per_cell = nppc,
	mass = 1.0,
	charge = 1.0,
	nb_density = 1.,
	mean_velocity=[0.,v0,0.],
	bc_part_type_xmin  = 'none',
	bc_part_type_xmax  = 'none',
	bc_part_type_ymin  = 'none',
	bc_part_type_ymax = 'none'
)
Species(
	species_type = 'eon',
	initPosition_type = 'random',
	initMomentum_type = 'mj',
	temperature = [Te,Te,Te],
	n_part_per_cell = nppc,
	mass = 1.0,
	charge = -1.0,
	nb_density = 1.,
	mean_velocity=[0.,v0,0.],
	bc_part_type_xmin  = 'none',
	bc_part_type_xmax  = 'none',
	bc_part_type_ymin  = 'none',
	bc_part_type_ymax = 'none'
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

DiagParticles(
 	output = "density",
 	every = globalEvery,
 	time_average = 1,
 	species = ["eon"],
 	axes = [
 		["gamma", 1., 1+15.*Te, 1000]
 	]
)
DiagParticles(
 	output = "density",
 	every = globalEvery,
 	time_average = 1,
 	species = ["eon"],
 	axes = [
 		["py", -2., 50., 100]
 	]
)
