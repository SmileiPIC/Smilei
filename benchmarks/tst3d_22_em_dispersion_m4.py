import numpy as np

Lx = 100
Ly = 150
Lz = 200
nx = 16
ny = 16
nz = 16

# Lx = 200
# Ly = 150
# Lz = 100
# nx = 256
# ny = 128
# nz = 192

nt = 1024
dt_multiplier = 0.2

dx = Lx / nx
dy = Ly / ny
dz = Lz / nz

def fx(x, y, z):
    if (dy < np.mod(y,Ly) <= 2.0*dy and dz < np.mod(z,Lz) <= 2.0*dz):
        return 1
    else:
        return 0

def fy(x, y, z):
    if (dz < np.mod(z,Lz) <= 2.0*dz and dx < np.mod(x,Lx) <= 2.0*dx):
        return 1
    else:
        return 0

def fz(x, y, z):
    if (dx < np.mod(x,Lx) <= 2.0*dx and dy < np.mod(y,Ly) <= 2.0*dy):
        return 1
    else:
        return 0

Main(
    geometry = "3Dcartesian",
    maxwell_solver = "M4",
    timestep_over_CFL = dt_multiplier, 
    number_of_timesteps = nt,
    cell_length = [dx, dy, dz],
    grid_length  = [Lx, Ly, Lz],
    number_of_patches = [2, 2, 2], 
    solve_poisson = False,
    print_every = 100,
)

DiagFields(
    every = [4],
    fields = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez'],
)

ExternalField(
    field = "Ex",
    profile = fx,
)

ExternalField(
    field = "Ey",
    profile = fy,
)

ExternalField(
    field = "Ez",
    profile = fz,
)
