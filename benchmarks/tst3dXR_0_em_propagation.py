# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi        # laser wavelength
t0 = l0                 # optical cycle
Lsim = [20.,50.]  # length of the simulation
Tsim = 1.*t0           # duration of the simulation
resx = 28.              # nb of cells in on laser wavelength
rest = 10#40.              # time of timestep in one optical cycle 

Main(
    geometry = "3drz",
    Nmode = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    number_of_patches = [ 1, 1 ],
    timestep = t0/rest,
    simulation_time = 100*t0/rest,
     
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["Buneman","Buneman"],
    ],
    
    random_seed = smilei_mpi_rank
)

LaserGaussian2D(
    a0              = 1.,
    omega           = 1.,
    focus           = [Lsim[0], Lsim[1]/2.],
    waist           = 8.,
    time_envelope   = tgaussian()
)


globalEvery = int(1)


#DiagScalar(every=globalEvery)

DiagFields(
    every = 100,
    fields = ['Ex_mode_1','Er_mode_1','Et_mode_1','Bx_mode_0','Br_mode_0','Bt_mode_0']
)

#DiagProbe(
#    every = 100,
#    number = [100, 100],
#    pos = [0., 10.*l0],
#    pos_first = [20.*l0, 0.*l0],
#    pos_second = [3.*l0 , 40.*l0],
#    fields = []
#)
#
#DiagProbe(
#    every = 10,
#    pos = [0.1*Lsim[0], 0.5*Lsim[1]],
#    fields = []
#)

