# ----------------------------------------------------------------------------------------
#                   SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField
#
#
#                            THERMAL PLASMA SIMULATION: HYDROGEN PLASMA
#                                    -----------------------
import math as m

### SCALING PARAMETER

SpaceSCALING = 0.5
TimeSCALING  = 1
NbDiagOut    = 32

# PLASMA PARAMETERS
# -----------------
n0      = 1.0                   # plasma density
lde     = 1./8.                 # here I fixe the Debye length (to ensure I have a resolution of 1./8., see below)
T0      = lde**2                # plasma temperature (computed from the Debye length)
nppc    = 64                    # number of particles per cell

# RESOLUTION, BOX SIZE and SIMULATION TIME
# ----------------------------------------

dx = dy = dz = lde             # this is chosen so that numerical heating is avoided and resolution = 1./16.
dt           = dx*0.5           # this is the magic time-step for numerical cherenkov

Lx = Ly = Lz = 16.*SpaceSCALING      # box size
t_sim   = 64.*TimeSCALING       # simulation time

# PATCH SIZE & LOAD BALANCING PARAMETERS
# --------------------------------------
PatchSizex = PatchSizey = PatchSizez = 8    # number of cells per patch (NB: clrw is also fixed to PatchSizex)
NPatchx    = NPatchy    = NPatchz = int( (Lx/dx)/PatchSizex )

#1 less cell per patch along Z so that PatchSizez + 1 can be divided by 4.
Lz = Lz - NPatchz*dz
PatchSizez = 7


LoadBalancingFreq = 320


# SMILEI VARIABLES
# ----------------

Main(
     geometry = '3Dcartesian',
     interpolation_order = 2,
     cell_length = [dx,dy,dz],
     grid_length  = [Lx,Ly,Lz],
     timestep = dt,
     simulation_time = t_sim/8./128.*120.,
     number_of_patches = [NPatchx,NPatchy,NPatchz],
     clrw = PatchSizex,
     maxwell_solver = 'Yee',
     EM_boundary_conditions = [ ['periodic'] ],
     random_seed = smilei_mpi_rank,
     print_every = int(t_sim/dt/100.),
     vecto = "normal"
)

LoadBalancing(
    every = LoadBalancingFreq,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

Species(
        name = 'ion',
        position_initialization = 'regular',
        momentum_initialization = 'mj',
        particles_per_cell = nppc,
        c_part_max = 1.0,
        mass = 1836.0,
        charge = 1.0,
        number_density = 1.0,
        mean_velocity = [0.,0.,0.],
        temperature = [T0],
        boundary_conditions = [
        ["periodic", "periodic"],
        ["periodic", "periodic"],
        ["periodic", "periodic"],
        ]
)

Species(
        name = 'eon',
        position_initialization = 'regular',
        momentum_initialization = 'mj',
        particles_per_cell = nppc,
        c_part_max = 1.0,
        mass = 1.0,
        charge = -1.0,
        number_density = 1.0,
        mean_velocity = [0.,0.,0.],
        temperature = [T0],
        boundary_conditions = [
        ["periodic", "periodic"],
        ["periodic", "periodic"],
        ["periodic", "periodic"],
        ]
)

# DIAGNOSTICS
#DiagScalar(every=int(10))  # at all time steps

#DiagFields(               # managed by NbDiagOut
#         every = int(t_sim/dt/NbDiagOut) ,
#         fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_ion','Rho_eon']
#)
