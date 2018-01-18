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
LoadBalancingFreq = 320


# SMILEI VARIABLES
# ----------------

Main(
     geometry = '3d3v',
     interpolation_order = 2,
     cell_length = [dx,dy,dz],
     sim_length  = [Lx,Ly,Lz],
     timestep = dt,
     sim_time = t_sim/8./128.*120.,
     number_of_patches = [NPatchx,NPatchy,NPatchz],
     clrw = PatchSizex,
     maxwell_sol = 'Yee',
     bc_em_type_x = ['periodic'],
     bc_em_type_y = ['periodic'],
     bc_em_type_z = ['periodic'],
     random_seed = smilei_mpi_rank,
     print_every = int(t_sim/dt/100.)
)

LoadBalancing(
    every = LoadBalancingFreq,
    coef_cell = 1.,
    coef_frozen = 0.1
)

Species(
        species_type = 'ion',
        initPosition_type = 'regular',
        initMomentum_type = 'mj',
        n_part_per_cell = nppc,
        c_part_max = 1.0,
        mass = 1836.0,
        charge = 1.0,
        nb_density = 1.0,
        mean_velocity = [0.,0.,0.],
        temperature = [T0],
        bc_part_type_xmin  = 'none',
        bc_part_type_xmax  = 'none',
        bc_part_type_ymin  = 'none',
        bc_part_type_ymax  = 'none',
        bc_part_type_zmin  = 'none',
        bc_part_type_zmax  = 'none'
)

Species(
        species_type = 'eon',
        initPosition_type = 'regular',
        initMomentum_type = 'mj',
        n_part_per_cell = nppc,
        c_part_max = 1.0,
        mass = 1.0,
        charge = -1.0,
        nb_density = 1.0,
        mean_velocity = [0.,0.,0.],
        temperature = [T0],
        bc_part_type_xmin  = 'none',
        bc_part_type_xmax  = 'none',
        bc_part_type_ymin  = 'none',
        bc_part_type_ymax  = 'none',
        bc_part_type_zmin  = 'none',
        bc_part_type_zmax  = 'none'
)

# DIAGNOSTICS
#DiagScalar(every=int(10))  # at all time steps

#DiagFields(               # managed by NbDiagOut
#         every = int(t_sim/dt/NbDiagOut) ,
#         fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_ion','Rho_eon']
#)





