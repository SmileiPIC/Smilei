def my_function(x) :
    return x**2


###################################################

mysim=smilei()

myspec1=species(mysim)
myspec1.species_type = 'ion'
myspec1.mass   = 100.0
myspec1.charge = 1
myspec1.n_part_per_cell = 2
myspec1.c_part_max = 1.0
myspec1.species_geometry = 'harris'
myspec1.density = 0.2
myspec1.dens_dbl_params = 0.2
myspec1.dens_length_y = (5.0, 150.0, 450.0)
myspec1.mvel_z_profile = 'harris'
myspec1.mean_velocity = (0.0, 0.0, 1.0)
myspec1.mvel_z_int_params = 0
myspec1.mvel_z_dbl_params = (0.2, 0.25, 0.025, 0.2)
myspec1.mvel_z_length_x = (10.0, 210.0, 430.0)
myspec1.mvel_z_length_y = (5.0,  150.0, 450.0)
myspec1.initialization_type  = 'mj'
myspec1.temperature = 0.02604166667
myspec1.bc_part_type_west  = 'none'
myspec1.bc_part_type_east  = 'none'
myspec1.bc_part_type_south = 'none'
myspec1.bc_part_type_north = 'none'
myspec1.initPosition_type = 'regular'
myspec1.initMomentum_type = 'mj'
myspec1.dynamics_type='norm'

myspec2=myspec1

smilei.species.append(myspec2)
myspec2

mysim.interpolation_order=2

mysim.sim_units='normalized'
mysim.dim='1d3v'


mysim.cell_length = 6.28/10000.
mysim.res_space   = 3.2
mysim.sim_length   = 3.2

mysim.bc_em_type_long  = 'periodic'
mysim.bc_em_type_trans = 'periodic'


###################################


mysim.profile=my_function

