# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Tsim = t0           # duration of the simulation
laser_fwhm = 19.80
dx = 0.2
dy = 0.5
dz = dy
nx = 800
ny = 100
nz = ny
Lsim = [dx*nx,ny*dy, nz*dz] # length of the simulation
dt = 0.18

Main(
    geometry = "3Dcartesian",
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [dx, dy, dz],
    grid_length  = Lsim,
    number_of_patches = [ 4, 4, 4 ],
    timestep = dt,
    simulation_time = 2.*nx*dt,
     
    EM_boundary_conditions = [["silver-muller"]],
    
    #print_every = 1,
)

MovingWindow(
    time_start = Main.grid_length[0]-50*dx,
    velocity_x = 0.9997
)

LaserGaussian3D(
    a0              = 1.,
    omega           = 1.,
    #focus           = [Lsim[0]/2., 0.],
    focus           = [0. , 0., 0.],
    waist           = 8.,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)
#LaserGaussian2D(
#    box_side         = "xmin",
#    a0              = 2.,
#    focus           = [0., 0.],
#    waist           = 26.16,
#    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
#)

Species( 
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 16,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.000001],
    pusher = "boris",    
    time_frozen = 0.0,
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
        ["remove", "remove"],
    ],
)

globalEvery = int(1)

list_fields = ['Ex','Ey','Ez','Rho','Jx', 'Jy', 'Jz']
#DiagScalar(every=globalEvery)

DiagFields(
    every = 100,
    fields = list_fields
)

#DiagProbe(
#    every = 10,
#    origin = [1., 10., 0.],
#    fields = []
#)
#DiagProbe(
#    every = 10,
#    origin = [0., 10., 0.],
#    corners = [[Lsim[0], 10., 0.]],
#    number=[100],
#    fields = []
#)
#DiagProbe(
#    every = 10,
#    origin = [0., -10., 0.],
#    corners = [[Lsim[0], -10., 0.]],
#    number=[100],
#    fields = []
#)


