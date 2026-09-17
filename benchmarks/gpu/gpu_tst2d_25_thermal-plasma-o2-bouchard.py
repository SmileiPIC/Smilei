# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

from math import pi, cos, sin, sqrt

l0 = 2.0*pi                # laser wavelength
t0 = l0                    # optical cycle
Lsim = [32.*l0,32.*l0]     # length of the simulation
resx,resy = 16,16                  # nb of cells in on laser wavelength

solver = 'Bouchard'

if solver == 'Bouchard' :
  rest = resx*2
  custom_oversize = 4
elif solver == 'Yee' :
  rest = resx*sqrt(2)/0.96
  custom_oversize = 2

Tsim = 1000.*t0/rest

Main(
    geometry = "2Dcartesian",
    maxwell_solver = solver,
    interpolation_order = 2,
    cell_length = [l0/resx,l0/resy],
    grid_length  = Lsim,
    number_of_patches = [ 1, 1 ],
    timestep = t0/rest,
    simulation_time = Tsim,
    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
    ],
    gpu_computing          = True,
    solve_poisson          = False,
    custom_oversize        = custom_oversize,
    print_every = 10,
    every_clean_particles_overhead = 100,   
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
	name = "eon",
	position_initialization = "random",
	momentum_initialization = "mj",
	particles_per_cell = 55,
	mass = 1.0,
	charge = -1.0,
        temperature = [1/511],
	charge_density = 1.,
	mean_velocity = [0.,0.,0.],
	boundary_conditions = [
		["thermalize", "thermalize"],["thermalize", "thermalize"]
	],
        thermal_boundary_temperature = [1/511],
	is_test = False
)

Species(
        name = "ion",
        position_initialization = "eon",
        momentum_initialization = "mj",
        particles_per_cell = 55,
        mass = 1836.0,
        charge = +1.0,
        temperature = [0.3/511],
        charge_density = 1.,
        mean_velocity = [0.,0.,0.],
        boundary_conditions = [
                ["thermalize", "thermalize"],["thermalize", "thermalize"]
        ],
        thermal_boundary_temperature = [0.3/511],
        is_test = False
)

Species(
        name = "eon2",
        position_initialization = "random",
        momentum_initialization = "mj",
        particles_per_cell = 9,
        mass = 1.0,
        charge = -1.0,
        temperature = [1/511],
        charge_density = 1.,
        mean_velocity = [0.,0.,0.],
        boundary_conditions = [
                ["thermalize", "thermalize"],["thermalize", "thermalize"]
        ],
        thermal_boundary_temperature = [1/511],
        is_test = False
)

Species(
        name = "ion2",
        position_initialization = "eon2",
        momentum_initialization = "mj",
        particles_per_cell = 9,
        mass = 1836.0,
        charge = +1.0,
        temperature = [0.3/511],
        charge_density = 1.,
        mean_velocity = [0.,0.,0.],
        boundary_conditions = [
                ["thermalize", "thermalize"],["thermalize", "thermalize"]
        ],
        thermal_boundary_temperature = [0.3/511],
        is_test = False
)


globalEvery = int(rest)

DiagScalar(
    every=globalEvery
)

DiagProbe(
    #name = "my probe diag",
    datatype = "float",
    every = 100,
    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz"],
    #origin   = [0., 0., 0.],
    origin   = [0., 0.],
    vectors = [
        #[0,0,Lsim[2]],
        #[0,Lsim[1],0],
        #[Lsim[0],0,0]
        [0,Lsim[1]],
        [Lsim[0],0]
    ],
    number = [
        #int(Lsim[2]/l0*resz),
        int(Lsim[1]/l0*resy),
        int(Lsim[0]/l0*resx)
    ]
)
