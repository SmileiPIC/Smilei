# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi		# laser wavelength
t0 = l0					# optical cicle
Lsim = [50.*l0,50.*l0]	# length of the simulation
Tsim = 400.*t0			# duration of the simulation
resx = 8.				# nb of cells in on laser wavelength
rest = 16.				# time of timestep in one optical cycle 

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 4 ,
    global_factor=[4,4], 
    norder = [2,2],
    is_spectral = True,
    is_pxr = True, 

    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4, 4  ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
    print_every=10,
     
    EM_boundary_conditions = [
        ['silver-muller','silver-muller'],
        ['silver-muller','silver-muller'],
    ],
    random_seed = smilei_mpi_rank
)

globalEvery = 20


Antenna(
    field='Jz',
    time_profile= lambda t: math.sin(2.*t/t0)*math.exp(-((t-10*t0)/(8*t0))**2),
    space_profile=gaussian(100., xfwhm=l0/resx, yfwhm=4*l0, xcenter=Main.grid_length[0]*0.2, ycenter=Main.grid_length[1]*0.5)
)


DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ez','Jz','Rho_electron','Rho_species1']
)

