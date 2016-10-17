# ----------------------------------------------------------------------------------------
#                     SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi    # laser wavelength
t0 = l0            # optical cicle
Lsim = [10.*l0,10.*l0]    # length of the simulation
Tsim = 10.*t0        # duration of the simulation
resx = 100.        # nb of cells in on laser wavelength
rest = 150.        # time of timestep in one optical cycle 

Main(
    geometry = "2d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx],
    sim_length  = Lsim,
    
    number_of_patches = [ 8, 8 ],
    
    timestep = t0/rest,
    sim_time = Tsim,
     
    bc_em_type_x = ['silver-muller'],
    bc_em_type_y = ['periodic'],
    
    random_seed = 0
)

LaserGaussian2D(
    boxSide         = "west",
    a0              = 150.,
    focus           = [10.*l0, 5.*l0],
    waist           = 2.0*l0,
    ellipticity     = 1.,
    time_envelope   = ttrapezoidal(slope1=t0)
)


Species(
    species_type = 'ion',
    initPosition_type = 'regular',
    initMomentum_type = 'cold',
    n_part_per_cell = 4,
    mass = 1836.0,
    charge = 1.0,
    nb_density = trapezoidal(100.0,xvacuum=l0,xplateau=0.44*l0),
    bc_part_type_xmin = 'refl',
    bc_part_type_xmax = 'refl',
    bc_part_type_ymin = 'none',
    bc_part_type_ymax = 'none'
)
Species(
    species_type = 'eon',
    initPosition_type = 'regular',
    initMomentum_type = 'mj',
    n_part_per_cell = 4,
    mass = 1.0,
    charge = -1.0,
    nb_density = trapezoidal(100.0,xvacuum=l0,xplateau=0.44*l0),
    temperature = [0.001],
    bc_part_type_xmin = 'refl',
    bc_part_type_xmax = 'refl',
    bc_part_type_ymin = 'none',
    bc_part_type_ymax = 'none'
)



globalEvery = int(rest/2.)

DiagScalar(every=globalEvery)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho_ion','Rho_eon']
)
