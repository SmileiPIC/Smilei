# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
dx = 0.251327
dtrans = 1.96349
dt = 0.96 * dx
nx =  960
ntrans = 256
Lx = nx * dx
Ltrans = ntrans * dtrans
npatch_x = 64
npatch_trans =32
Nit = 2000


Main(
    geometry = "AMcylindrical",
    number_of_AM=2,
    interpolation_order = 2,
    timestep = dt,
    simulation_time = dt*Nit,
    cell_length  = [dx, dtrans],
    grid_length = [ Lx,  Ltrans],
    number_of_patches = [npatch_x, npatch_trans],
    uncoupled_grids=True,
    clrw = 5,
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    random_seed = smilei_mpi_rank,
    solve_poisson = False,
    print_every = 100
)

MovingWindow(
    time_start = Main.grid_length[0] - 50*dx, #Leaves 2 patches untouched, in front of the laser.
    velocity_x = 0.996995486
)

ne = 0.0045
begin_upramp = 10.  #Density is 0 before that and up ramp starts.
Lupramp = 100. #Length of the upramp 
Lplateau = 15707.  #Length of the plateau 
Ldownramp = 2356.19 #Length of the down ramp
xplateau = begin_upramp + Lupramp # Start of the plateau
begin_downramp = xplateau + Lplateau # Beginning of the output ramp. 
finish = begin_downramp + Ldownramp # End of plasma

g = polygonal(xpoints=[begin_upramp, xplateau, begin_downramp, finish], xvalues=[0, ne, ne, 0.])

def my_profile(x,y):
    return g(x,y)

Species( 
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    ionization_model = "none",
    particles_per_cell = 30,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = my_profile,  # Here absolute value of the charge is 1 so charge_density = nb_density
    mean_velocity = [0., 0., 0.],
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
    focus           = [10.,0.],  
    waist           = 120.,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)

DiagProbe(
	every = 1000,
	origin = [0., 2*dtrans, 0.],
	corners = [
              [Main.grid_length[0], 2*dtrans, 0.]
                  ],
	number = [nx],
)

#DiagProbe(
#	every = 1000,
#	origin = [0., -Main.grid_length[1], 0.],
#	corners =  [
#           [Main.grid_length[0], -Main.grid_length[1], 0.],
#	   [0., Main.grid_length[1], 0.],
#                   ],
#	number = [nx,2*ntrans],
#)

