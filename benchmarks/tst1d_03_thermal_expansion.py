# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math, random

Te_keV = 1.              # electron temperature in keV
Te  = Te_keV/511.        # Te normalised in mec^2 (code units)
vth = math.sqrt(Te)      # normalised thermal velocity
Ld    = vth              # Debye length in normalised units
dx  = Ld/10.             # spatial resolution
Lsim = 40.*Ld            # simulation length
tsim = 50.               # duration of the simulation

mi = 100.0               # ion mass (use reduced one to accelerate computation)
cs = math.sqrt(Te/mi)    # ion acoustic velocity (normalised to c)
Uion = mi*cs*cs          # mean energy used to compute the ion spectrum


Main(
    geometry = "1Dcartesian",
    
    interpolation_order = 2,
    
    timestep = 0.95*dx,
    simulation_time = tsim,
    
    cell_length = [dx],
    grid_length  = [Lsim],
    
    number_of_patches = [ 16 ],
    
    EM_boundary_conditions = [ ['silver-muller','silver-muller'] ] ,
    
    random_seed = int(random.getrandbits(32))+smilei_mpi_rank
)

Species(
    name = 'ion',
    position_initialization = 'random',
    momentum_initialization = 'mj',
    particles_per_cell = 10,
    mass = mi, 
    charge = 1.0,
    number_density = trapezoidal(1., xplateau=20.*Ld),
    temperature = [1.e-6],
    thermal_boundary_temperature = [1.e-6],
    thermal_boundary_velocity = [0.,0.,0.],
    boundary_conditions = [
    	["thermalize", "reflective"],
    ],
)
Species(
    name = 'eon',
    position_initialization = 'random',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = 400,
    mass = 1.0,
    charge = -1.0,
    number_density = trapezoidal(1., xplateau=20.*Ld),
    temperature = [Te],
    thermal_boundary_temperature = [Te],
    thermal_boundary_velocity = [0.,0.,0.],
    boundary_conditions = [
    	["thermalize", "reflective"],
    ],
)

LoadBalancing(
    every = 100
)

every=200

DiagScalar(every = every)#, vars=['Utot','Ubal_norm','Uelm','Ukin','Ukin_ion','Ukin_eon'])    


DiagFields(
    every = every,
    fields = ['Ex','Rho_ion','Rho_eon']
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = every,
    species = ["ion"],
    axes = [
        ["x", 0., Lsim, 50],
        ["px", -3.*cs, 3.*cs, 100]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = every,
    time_average = 1,
    species = ["ion"],
    axes = [
        ["ekin", 0.01*Uion, 10*Uion, 100, "logscale"]
    ]
)


DiagParticleBinning(
    deposited_quantity = "weight",
    every = 10,
    time_average = 10,
    species = ["eon"],
    axes = [
        ["ekin", 0.1*Te, 20*Te, 30, "logscale", "edge_inclusive"]
    ]
)
