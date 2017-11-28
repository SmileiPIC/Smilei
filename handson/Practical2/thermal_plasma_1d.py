# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------


# DEFINING MY OWN VARIABLES
# Simple Python script
# -------------------------

import math 

ne = 1.                    # electron density (code units => reference density)
Ld = 1./64.                # Debye length (in code units => i.e. in skin-depth as ne=1)

Te = ne*Ld**2              # Te normalised in mec^2 (code units)
vth = math.sqrt(Te)        # normalised thermal velocity

dx = 1.*Ld                # spatial resolution
dt = 0.95*dx               # timestep
Lx = 512.                  # simulation length
tsim = 1024.               # duration of the simulation

nppc = 16                  # number of particle-per-cell

diagEvery   = int(64./dt)  # frequency of outputs for DiagField


# DEFINING SMILEI's VARIABLES 
# All in "blocks"
# ---------------------------

Main(
    geometry = "1Dcartesian",
    
    interpolation_order = 2,
    
    timestep = dt,
    simulation_time = tsim,
    
    cell_length = [dx],
    grid_length  = [Lx],
    
    number_of_patches = [ 32 ],
    
    EM_boundary_conditions = [ ['periodic'] ] ,
    
    random_seed = smilei_mpi_rank
)

Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = nppc,
    mass = 1836., 
    charge = 1.0,
    number_density = ne,
    temperature = [Te],
    boundary_conditions = [
    	['periodic'],
    ],
    time_frozen = 2.*tsim
)

Species(
    name = 'eon',
    position_initialization = 'random',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = nppc,
    mass = 1.0,
    charge = -1.0,
    number_density = ne,
    temperature = [Te],
    thermal_boundary_temperature = [Te],
    boundary_conditions = [
    	['periodic'],
    ]
)

LoadBalancing(
    every = 100
)

### DIAGNOSTICS

DiagScalar(every = 2.5/dt)


DiagFields(
    every = diagEvery,
    fields = ['Ex','Ey','Rho_ion','Rho_eon']
)


DiagParticleBinning(
    deposited_quantity = "weight",
    every = diagEvery,
    time_average = 1,
    species = ["eon"],
    axes = [
        ["ekin", 0., 4*Te, 32]
    ]
)


