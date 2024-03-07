# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

from numpy import pi, linspace, random, meshgrid, cos, sin, vstack
dx = 0.251327
dtrans = 1.96349
dt = 0.96 * dx
nx =  960
ntrans = 256
Lx = nx * dx
Ltrans = ntrans * dtrans
npatch_x = 16
npatch_trans = 8
Nit = 2000

ne = 0.0045
begin_upramp = 10.  #Density is 0 before that and up ramp starts.
Lupramp = 100. #Length of the upramp 
Lplateau = 2000.  #Length of the plateau 
Ldownramp = 2356.19 #Length of the down ramp
xplateau = begin_upramp + Lupramp # Start of the plateau
begin_downramp = xplateau + Lplateau # Beginning of the output ramp. 
finish = begin_downramp + Ldownramp # End of plasma

# hdf file to contain all particles
def preprocess():
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    if comm.Get_rank() > 0: return
    nppc_x = 2
    nppc_trans = 1
    W = ne * 2.*pi * dx * dtrans / (nppc_x*nppc_trans)
    ncell_x = int((begin_downramp - xplateau) / dx)
    x = linspace( xplateau, begin_downramp, ncell_x*nppc_x )
    r = linspace( 0.1*dtrans, Ltrans-0.1*dtrans, ntrans*nppc_trans )
    x, r = meshgrid( x, r )
    x = x.flatten()
    r = r.flatten()
    theta = 2.*pi*random.random(len(x))
    y = r * cos(theta)
    z = r * sin(theta)
    w = W * r
    from h5py import File
    with File("plasma.h5","w") as f:
        f["position/x"] = x
        f["position/y"] = y
        f["position/z"] = z
        f["weight"] = w
    with File("momentum.h5","w") as f:
        f["momentum/x"] = x*0.
        f["momentum/y"] = y*0.
        f["momentum/z"] = z*0.

Main(
    geometry = "AMcylindrical",
    number_of_AM=2,
    interpolation_order = 2,
    timestep = dt,
    simulation_time = dt*Nit,
    cell_length  = [dx, dtrans],
    grid_length = [ Lx,  Ltrans],
    number_of_patches = [npatch_x, npatch_trans],
    cluster_width = 5,
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    solve_poisson = False,
    print_every = 100
)

MovingWindow(
    time_start = Main.grid_length[0] - 50*dx, #Leaves 2 patches untouched, in front of the laser.
    velocity_x = 0.996995486
)

Species(
    name = "electron",
    position_initialization = "plasma.h5",
    momentum_initialization = "momentum.h5",
    ionization_model = "none",
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    # mean_velocity = [0., 0., 0.],
    time_frozen = 0.0,
    boundary_conditions = [
    	["remove", "remove"],
    	["reflective", "remove"],
    ],
)

laser_fwhm = 82. 
LaserGaussianAM(
    box_side         = "xmin",
    a0              = 2.,
    focus           = [10.],  
    waist           = 120.,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)

DiagProbe(
	every = 100,
	origin = [0., -Ltrans, 0.],
	corners = [
              [Lx, -Ltrans, 0.],
              [0.,  Ltrans, 0.]
                  ],
	number = [nx, ntrans],
)

DiagParticleBinning(
    every = 100,
    species = ["electron"],
    deposited_quantity = "weight",
    axes = [
        ["moving_x",0.,Lx,10],
        ["y",-Ltrans,Ltrans,6],
        ["z",-0.2*Ltrans, 0.2*Ltrans, 1]
    ]
)

DiagPerformances(
	every = 1000,
)

