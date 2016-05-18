# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math
l0 = 2.0*math.pi	# wavelength in normalized units
t0 = l0				# optical cycle in normalized units
rest = 102.0		# nb of timestep in 1 optical cycle
resx = 100.0		# nb cells in 1 wavelength

Main(
    geometry = "1d3v",
    interpolation_order = 2,
    
    cell_length = [l0/resx],
    sim_length  = [4.0*l0],
    
    number_of_patches = [ 4 ],
    
    timestep = t0/rest,
    sim_time = 4.0*t0,
    
    bc_em_type_x = ['silver-muller'],
    
    random_seed = 0,
    
    print_every = int(rest/2.0)
)

Laser(
    omega          = 1.,
    chirp_profile  = tpolynomial(order2=0.005),
    time_envelope  = tgaussian(),
    space_envelope = [1., 0.],
)

DiagScalar(
    every = 1
)

DiagFields(
    every = int(rest/2.0),
    fields = ['Ex','Ey','Ez','By_m','Bz_m']
)

DiagProbe(
    every = 1, 
    pos = [0.0]   
)

DiagProbe(
    every = 5,
    pos = [0.0],
    pos_first = Main.sim_length,
    number = [1000]
)
