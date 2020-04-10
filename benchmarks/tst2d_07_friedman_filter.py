
# ----------------------------------------------------------------------------------------
#                   SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticleBinning,
#           DiagScalar, DiagPhase or ExternalField
#
import math as m
    #Checkpoints(
            #restart_dir = "/ccc/scratch/cont003/gen7678/grassia/Results/Test_S45Laser",
            #        dump_step = 20000,
            #dump_minutes = 240.,
            #dump_deflate = 0,
            #exit_after_dump = True,
            #keep_n_dumps = 2
#)

## CHARACTERISTIC DISTANCES, TEMPERATURE & TIMES
# ---------------------------------------------

# plasma parameters
n0      = 1.0
Mi      = 1.0

Te      = 0.0001

g0      = 10.
v0      = m.sqrt(1.-1./g0**2)
nppc    = 10
sigma   = 0.
B0      = m.sqrt(sigma*2.*g0)
E0      = v0*B0
# simulation box


dx      = 1./4.
dt      = dx*0.5

Lx    = 1216./8.
Ly    = 256/8.
t_sim = 500./8.


def B(x,y):
    return B0
def E(x,y):
    return E0

LoadBalancing(
    every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

Main(
    #dim: Geometry of the simulation
    geometry = '2Dcartesian',
    
    # order of interpolation
    interpolation_order = 2,
    
    # SIMULATION BOX : for all space directions (use vector)
    # cell_length: length of the cell
    # grid_length: length of the simulation in units of the normalization wavelength
    cell_length = [dx,dx],
    grid_length  = [Lx,Ly],
    maxwell_solver = 'Yee',
    
    number_of_patches = [16,16],
    clrw = 1,
    
    # SIMULATION TIME
    # timestep: duration of the timestep
    # simulation_time: duration of the simulation in units of the normalization period
    timestep = dt,
    simulation_time = t_sim,
    
    # ELECTROMAGNETIC BOUNDARY CONDITIONS
    #  periodic = periodic BC (using MPI topology)
    #  silver-muller = injecting/absorbing BC
    #  reflective = consider the ghost-cells as a perfect conductor
    EM_boundary_conditions = [
        ['silver-muller'],
        ['periodic']
    ],
    
    print_every = int(t_sim/dt/50.),
    
    # regular seed
    # this is used to regularize the regular number generator
    random_seed = smilei_mpi_rank
)

CurrentFilter(
    model = "binomial",
    passes = [3,3],
)

FieldFilter(
    model = "Friedman",
    theta = 0.3,
)

Species(
    name = 'pos',
    position_initialization = 'random',
    momentum_initialization = 'mj',
    ionization_model = 'none',
    particles_per_cell = nppc,
    c_part_max = 1.0,
    mass = Mi,
    charge = 1.0,
    number_density = n0,
    mean_velocity = [v0,0.,0.],
    temperature = [Te],
    thermal_boundary_temperature = [0.],
    thermal_boundary_velocity = [0.,0.,0.],
    time_frozen = 0.,
    boundary_conditions = [["reflective"], ["periodic"]],
)

Species(
    name = 'eon',
    position_initialization = 'random',
    momentum_initialization = 'mj',
    ionization_model = 'none',
    particles_per_cell = nppc,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    number_density = n0,
    mean_velocity = [v0,0.,0.],
    temperature = [Te],
    thermal_boundary_temperature = [0.],
    thermal_boundary_velocity = [0.,0.,0.],
    time_frozen = 0,
    boundary_conditions = [["reflective"], ["periodic"]],
)

#ExternalField(
#         field = 'Bz',
#         profile = B
#         )

#ExternalField(
#         field = 'Ey',
#         profile = E
#         )

globalEvery = int(5./dt)

# scalar diagnostics
DiagScalar(every=int(1))

DiagFields(
    every = globalEvery ,
    fields = ['Ey','Bz','Rho_eon','Jy_eon','Jy_pos']
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    time_average = 1,
    species = ["eon"],
    axes = [
        ["x", Lx/2., Lx, 600 ],
        ["gamma", 1, 500, 500]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    time_average = 1,
    species = ["eon"],
    axes = [
        ["x", Lx/2., Lx, 600 ],
        ["px",-75, 75, 400],
        ["py",-75, 75, 400]
    ]
)


