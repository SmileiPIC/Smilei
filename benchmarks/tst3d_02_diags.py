import math

L0 = 2.*math.pi # Wavelength in PIC units

def exp(t):
	return math.exp(-t)


Main(
	geometry = "3Dcartesian",
	
	interpolation_order = 2,
	
	timestep = 0.005 * L0,
	simulation_time  = .24 * L0,
	
	cell_length = [0.01 * L0]*3,
	grid_length  = [1. * L0]*3,
	
	number_of_patches = [ 4 ]*3,
	
	EM_boundary_conditions = [ ["periodic"] ],
	print_every = 10,
	
#)


nportion = 5
portion_width = Main.grid_length[0]/nportion

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
	particles_per_cell = 10,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = lambda x,y,z : poly3(x,y,z) * math.exp(-(y-0.5*L0)**2) * math.exp(-(z-0.5*L0)**2),
	mean_velocity = [0.00001, 0.00001, 0.00001],
	pusher = "boris",
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
	particles_per_cell = 10,
	c_part_max = 1.0,
	mass = 1.0,
	charge = 1.0,
	number_density = lambda x,y,z : poly4(x,y,z) * math.exp(-(y-0.5*L0)**2) * math.exp(-(z-0.5*L0)**2),
	mean_velocity = [0.00001, 0.00001, 0.00001],
	pusher = "boris",
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

axes = [
	["x"        , 0., Main.grid_length[0], 10],
	["y"        , 0., Main.grid_length[1], 10],
	["z"        , 0., Main.grid_length[2], 10],
	["px"       , -0.0001, 0.0001, 30],
	["py"       , -0.0001, 0.0001, 30],
	["pz"       , -0.0001, 0.0001, 30],
	["p"        , 0., 0.0001, 30],
	["gamma"    , 1., 1.+1e-5, 30], # precision insufficient
	["ekin"     , 0., 1e-9, 30],
	["vx"       , -0.0001, 0.0001, 30],
	["vy"       , -0.0001, 0.0001, 30],
	["vz"       , -0.0001, 0.0001, 30],
	["v"        , 0., 0.0001, 30],
	["vperp2"   , 0., 1e-9, 30],
	["charge"   , -5., 5., 10],
	[lambda particles: particles.x, 0., Main.grid_length[0], 10],
]

for axis in axes:
	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 10,
		species = ["test0"],
		axes = [axis]
	)

quantities = [
	"weight"          ,
	"weight_charge"   ,
	"weight_charge_vx",
	"weight_charge_vy",
	"weight_charge_vz",
	"weight_ekin"     ,
	"weight_p"        ,
	"weight_px"       ,
	"weight_py"       ,
	"weight_pz"       ,
	"weight_vx_px"    ,
	"weight_vy_py"    ,
	"weight_vz_pz"    ,
	"weight_vx_py"    ,
	"weight_vx_pz"    ,
	"weight_vy_pz"    ,
	"weight_ekin_vx"  ,
	lambda particles: particles.weight
]

for quantity in quantities:
	DiagParticleBinning(
		deposited_quantity = quantity,
		every = 10,
		species = ["test0"],
		axes = [["x" , 0., Main.grid_length[0], 10]]
	)
