import math

L0 = 2.*math.pi # Wavelength in PIC units

def exp(t):
	return math.exp(-t)


Main(
	geometry = "3d3v",
	
	interpolation_order = 2,
	
	timestep = 0.005 * L0,
	sim_time  = .24 * L0,
	
	cell_length = [0.01 * L0]*3,
	sim_length  = [1. * L0]*3,
	
	number_of_patches = [ 4 ]*3,
	
	EM_boundary_conditions = [ ["periodic"] ],
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
	name = "test0",
	position_initialization = "random",
	momentum_initialization = "cold",
	n_part_per_cell = 10,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	nb_density = lambda x,y,z : poly3(x,y,z) * math.exp(-(y-0.5*L0)**2) * math.exp(-(z-0.5*L0)**2),
	mean_velocity = [0.00001, 0.00001, 0.00001],
	dynamics_type = "norm",
	time_frozen = 1.30, # Move only after timestep 40
	boundary_conditions = [
		["periodic", "periodic"],
		["periodic", "periodic"],
		["periodic", "periodic"],
	],
)
iportion = 4
poly4 = polygonal( xpoints=[(iportion+i/3.)*portion_width for i in range(4)], xvalues=[0., 2., -1, 0.] )
Species( 
	name = "test1",
	position_initialization = "random",
	momentum_initialization = "cold",
	n_part_per_cell = 10,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	nb_density = lambda x,y,z : poly4(x,y,z) * math.exp(-(y-0.5*L0)**2) * math.exp(-(z-0.5*L0)**2),
	mean_velocity = [0.00001, 0.00001, 0.00001],
	dynamics_type = "norm",
	time_frozen = 1000000.0,
	boundary_conditions = [
		["periodic", "periodic"],
		["periodic", "periodic"],
		["periodic", "periodic"],
	],
)

DiagScalar(
    every = 10, 
)
DiagFields(
    every = 45, 
)

DiagProbe(
    every = 40,
    number = [10],
    origin = [0., L0/2, L0/2],
    corners = [
        [L0, L0/2, L0/2]
    ],
    fields = []
)

