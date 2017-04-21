import math

L0 = 2.*math.pi # Wavelength in PIC units

def exp(t):
	return math.exp(-t)


Main(
	geometry = "3d3v",
	
	interpolation_order = 2,
	
	timestep = 0.005 * L0,
	sim_time  = .2 * L0,
	
	cell_length = [0.01 * L0]*3,
	sim_length  = [1. * L0]*3,
	
	number_of_patches = [ 4 ]*3,
	
	bc_em_type_x = ["periodic"],
	bc_em_type_y = ["periodic"],
	bc_em_type_z = ["periodic"],
	print_every = 10
)


nportion = 5
portion_width = Main.sim_length[0]/nportion

iportion = 0
poly0 = polygonal( xpoints=[(iportion+i/3.)*portion_width for i in range(4)], xvalues=[0., 2., -1, 0.] )
Antenna(
	field = "Jx",
	space_profile = lambda x,y,z : poly0(x,y,z) * math.exp(-0.1*(y-0.5*L0)**2) * math.exp(-0.1*(z-0.5*L0)**2),
	time_profile = exp
)
iportion = 1
poly1 = polygonal( xpoints=[(iportion+i/3.)*portion_width for i in range(4)], xvalues=[0., 2., -1, 0.] )
Antenna(
	field = "Jy",
	space_profile = lambda x,y,z : poly1(x,y,z) * math.exp(-0.1*(y-0.5*L0)**2) * math.exp(-0.1*(z-0.5*L0)**2),
	time_profile = exp
)
iportion = 2
poly2 = polygonal( xpoints=[(iportion+i/3.)*portion_width for i in range(4)], xvalues=[0., 2., -1, 0.] )
Antenna(
	field = "Jz",
	space_profile = lambda x,y,z : poly2(x,y,z) * math.exp(-0.1*(y-0.5*L0)**2) * math.exp(-0.1*(z-0.5*L0)**2),
	time_profile = exp
)

iportion = 3
poly3 = polygonal( xpoints=[(iportion+i/3.)*portion_width for i in range(4)], xvalues=[0., 2., -1, 0.] )
Species( 
	species_type = "test0",
	initPosition_type = "random",
	initMomentum_type = "cold",
	n_part_per_cell = 10,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	nb_density = lambda x,y,z : poly3(x,y,z) * math.exp(-(y-0.5*L0)**2) * math.exp(-(z-0.5*L0)**2),
	dynamics_type = "norm",
	time_frozen = 1000000.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none",
	bc_part_type_ymin = "none",
	bc_part_type_ymax = "none",
	bc_part_type_zmin = "none",
	bc_part_type_zmax = "none"
)
iportion = 4
poly4 = polygonal( xpoints=[(iportion+i/3.)*portion_width for i in range(4)], xvalues=[0., 2., -1, 0.] )
Species( 
	species_type = "test1",
	initPosition_type = "random",
	initMomentum_type = "cold",
	n_part_per_cell = 10,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	nb_density = lambda x,y,z : poly4(x,y,z) * math.exp(-(y-0.5*L0)**2) * math.exp(-(z-0.5*L0)**2),
	dynamics_type = "norm",
	time_frozen = 1000000.0,
	bc_part_type_xmin = "none",
	bc_part_type_xmax = "none",
	bc_part_type_ymin = "none",
	bc_part_type_ymax = "none",
	bc_part_type_zmin = "none",
	bc_part_type_zmax = "none"
)

DiagScalar(
    every = 10, 
)
DiagFields(
    every = 10, 
)
