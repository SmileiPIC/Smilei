# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

# 15 Gbytes GPU Memory test

from math import pi, cos, sin, sqrt
from numpy import s_
import numpy as np

# Physical Params/ Units Library

c                 = 299792458
electron_mass     = 9.10938356e-31
electron_charge   = 1.60217662e-19
lambda_r          = 0.351e-6

mi_over_me        = 1836

wr                = 2*math.pi*c/lambda_r

# ----

l0 = 1.0                # laser wavelength
t0 = l0                 # optical cycle
k0 = 1.0
w0 = 1.0

Lsim = [64.*l0,16.*l0,16.*l0]     # length of the simulation
resx = 8                                    # nb of cells in on laser wavelength
resy,resz = resx,resx

# ----

dx,dy,dz = l0/resx,l0/resy,l0/resz
Lx,Ly,Lz = Lsim[0],Lsim[1],Lsim[2]

solver = 'Bouchard'

if solver == 'Bouchard' :
  rest = resx*2
  custom_oversize = 4
elif solver == 'Yee' :
  rest = resx*np.sqrt(3)/0.96
  custom_oversize = 2

Tsim = 1000*t0/rest              # duration of the simulation

Main(
    geometry               = "3Dcartesian",
    interpolation_order    = 2,
    maxwell_solver         = solver,
    cell_length            = [dx,dy,dz],
    grid_length            = Lsim,
    number_of_patches      = [ 1, 1, 1 ],
    timestep               = t0/rest,
    simulation_time        = Tsim,
    EM_boundary_conditions = [ ['silver-muller','silver-muller'] , ['silver-muller','silver-muller'] , ['silver-muller','silver-muller'] ] ,
    #patch_arrangement      = "linearized_XYZ",
    #patch_arrangement      = "linearized_ZYX",
    gpu_computing          = True,
    solve_poisson          = False,
    custom_oversize        = custom_oversize,
    print_every = 10,
    every_clean_particles_overhead = 200,
)

# Vectorization(
#     mode = "on",
#     reconfigure_every = 20,
#     initial_mode = "on"
# )

# LoadBalancing(
#     initial_balance = True,
#     every = 250,
#     cell_load = 1.,
#     frozen_particle_load = 0.1
# )

Species(
	name                    = "eon",
	particles_per_cell      = 3,
	position_initialization = "random",
	momentum_initialization = "mj",
	mass                    = 1.0,
	charge                  = -1.0,
        temperature             = [1/511],
	charge_density          = 0.05,
	mean_velocity           = [0.,0.,0.],
	boundary_conditions = [
            ["thermalize","thermalize"],["thermalize","thermalize"],["thermalize","thermalize"]
	],
        thermal_boundary_temperature = [1/511],
	is_test = False
)

Species(
        name                    = "ion",
        particles_per_cell      = 3,
        position_initialization = "eon",
        momentum_initialization = "mj",
        mass                    = 1.0,
        charge                  = +1.0,
        temperature             = [0.3/511],
        charge_density          = 0.05,
        mean_velocity           = [0.,0.,0.],
        boundary_conditions = [
                ["thermalize","thermalize"],["thermalize","thermalize"],["thermalize","thermalize"]
        ],
        thermal_boundary_temperature = [0.3/511],
        is_test = False
)

Species(
        name                    = "eon2",
        particles_per_cell      = 1,
        position_initialization = "random",
        momentum_initialization = "mj",
        mass                    = 1.0,
        charge                  = -1.0,
        temperature             = [1/511],
        charge_density          = 0.05,
        mean_velocity           = [0.,0.,0.],
        boundary_conditions = [
            ["thermalize","thermalize"],["thermalize","thermalize"],["thermalize","thermalize"]
        ],
        thermal_boundary_temperature = [1/511],
        is_test = False
)

Species(
    name                    = "ion2",
    particles_per_cell      = 1,
        position_initialization = "eon2",
        momentum_initialization = "mj",
    mass                    = 1.0,
    charge                  = +1.0,
    temperature             = [0.3/511],
    charge_density          = 0.05,
    mean_velocity           = [0.,0.,0.],
    boundary_conditions = [
        ["thermalize","thermalize"],["thermalize","thermalize"],["thermalize","thermalize"]
    ],
    thermal_boundary_temperature = [0.3/511],
    is_test = False
)

globalEvery = int(2)

DiagProbe(
    #name = "my probe diag",
    datatype = "float",
    every = 100,
    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz"],
    origin   = [0., 0., 0.],
    vectors = [
        [0,0,Lsim[2]],
        [0,Lsim[1],0],
        [Lsim[0],0,0]
    ],
    number = [
        int(Lsim[2]/l0*resz),
        int(Lsim[1]/l0*resy),
        int(Lsim[0]/l0*resx)
    ]
)

DiagScalar(
    every=globalEvery
)
