import numpy as np

Lx = 200
Ly = 150
nx = 64
ny = 32

nt = 1024
dt_multiplier = 0.2

dx = Lx / nx
dy = Ly / ny

def f(x, y):
    if dx < np.mod(x,Lx) <= 2.0*dx and dy < np.mod(y,Ly) <= 2.0*dy : 
        return 1
    else: 
        return 0

Main(
    geometry = "2Dcartesian",
    maxwell_solver = "M4",
    timestep_over_CFL = dt_multiplier, 
    number_of_timesteps = nt,
    cell_length = [dx, dy],
    grid_length  = [Lx, Ly],
    number_of_patches = [8, 4], 
    solve_poisson = False,
    print_every = 100,
)

DiagFields(
    every = [4],
    fields = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez'],
)

ExternalField(
    field = "Ez",
    profile = f,
)
