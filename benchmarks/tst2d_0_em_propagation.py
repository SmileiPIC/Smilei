# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi		# laser wavelength
t0 = l0					# optical cicle
Lsim = [20.*l0,50.*l0]	# length of the simulation
Tsim = 50.*t0			# duration of the simulation
resx = 10.				# nb of cells in on laser wavelength
rest = 30.				# time of timestep in one optical cycle 

Main(
    geometry = "2d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    sim_length  = Lsim,
    
    number_of_patches = [ 8, 4 ],
    
    timestep = t0/rest,
    sim_time = Tsim,
     
    bc_em_type_x = ['silver-muller'],
    bc_em_type_y = ['silver-muller'],
    
    random_seed = 0
)

LaserGaussian2D(
    a0              = 1.,
    omega           = 1.,
    focus           = [Lsim[0], Lsim[1]/2.],
    waist           = 8.,
    incidence_angle = 0.5,
    time_envelope   = tgaussian()
)


globalEvery = int(rest)

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez']
)

