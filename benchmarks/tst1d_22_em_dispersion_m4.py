import numpy as np

Lx = 200
nx = 64

nt = 1024
dt_multiplier = 0.2

dx = Lx / nx

def f(x):
    if dx < np.mod(x,Lx) <= 2.0*dx:  return 1
    else: return 0

Main(
    geometry = "1Dcartesian",
    maxwell_solver = "M4",
    timestep_over_CFL = dt_multiplier, 
    number_of_timesteps = nt,
    cell_length = [dx],
    grid_length  = [Lx],
    number_of_patches = [8], 
    solve_poisson = False,
    print_every = 100,
)

DiagFields(
    every = [4],
    fields = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez'],
)

ExternalField(
    field = "Ey",
    profile = f,
)

ExternalField(
    field = "Ez",
    profile = f,
)
