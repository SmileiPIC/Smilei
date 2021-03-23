import math

L0 = 2.*math.pi # Wavelength in PIC units
Lcell = [0.01*L0, 0.01*L0]
Lsim = [1.*L0, 1.*L0]

# Make a profile in a file
def preprocess():
	if smilei_mpi_rank == 0:
		from numpy import linspace, meshgrid, ceil, exp
		import h5py
		x = linspace(0., Lsim[0], int(ceil(Lsim[0]/Lcell[0])) + 1 )
		y = linspace(0., Lsim[1], int(ceil(Lsim[1]/Lcell[1])) + 1 )
		xx, yy = meshgrid(x, y)
		with h5py.File("test_profile.h5","w") as f:
			g = f.create_group("some_group")
			g.create_dataset("the_profile", data = exp(-((xx-0.4*L0)/(0.1*L0))**2 -((yy-0.6*L0)/(0.1*L0))**2) )


Main(
	geometry = "2Dcartesian",
	
	interpolation_order = 2,
	
	timestep = 0.005 * L0,
	simulation_time  = 0.01 * L0,
	
	cell_length = Lcell,
	grid_length  = Lsim,
	
	number_of_patches = [ 4 ]*2,
	
	time_fields_frozen = 10000000.,
	
	EM_boundary_conditions = [
		["periodic"],
		["periodic"],
	], 
	print_every = 10,
	solve_poisson = False,
	
    random_seed = smilei_mpi_rank
)

def custom(x, y):
	if y<L0/2.: return 0.
	else: return math.exp(-y)

profiles = {
"constant"   :constant   (1.),
"trapezoidal":trapezoidal(1., xvacuum=0.1*L0, xplateau=0.4*L0, xslope1=0.1*L0, xslope2=0.1*L0, yvacuum=0.3*L0, yplateau=0.1*L0, yslope1=0.2*L0, yslope2=0.2*L0),
"gaussian"   :gaussian   (1., xvacuum=0.1*L0, xlength =0.5*L0, xfwhm=0.2*L0, xcenter=0.2*L0, xorder=2, yvacuum=0.2*L0, ylength =0.6*L0, yfwhm=0.4*L0, ycenter=0.5*L0, yorder=4),
"polygonal"  :polygonal  (xpoints=[0.1*L0, 0.2*L0, 0.4*L0, 0.8*L0], xvalues=[1.,0.5,0.8, 0.1]),
"cosine"     :cosine     (1., xamplitude=0.4, xvacuum=0.3*L0, xlength=0.4*L0, xphi=0.1*L0, xnumber=3, yamplitude=0.2, yvacuum=0.2*L0, ylength=0.6*L0, yphi=0.3*L0, ynumber=10),
"polynomial" :polynomial (x0=0.4*L0, y0=0.5*L0, order0=1., order1=[-1./L0,-0.1/L0], order2=[(1./L0)**2,(0.1/L0)**2,(0.1/L0)**2]),
"custom"     :custom,
"file"       :"test_profile.h5/some_group/the_profile"
}

for name, profile in profiles.items():
	Species(
		name = name,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= 100,
		mass = 1.0,
		charge = 1.0,
		number_density = profile,
		time_frozen = 10000.0,
		boundary_conditions = [
			["periodic", "periodic"],
			["periodic", "periodic"],
		],
	)


# NON-RELATIVISTIC MAXWELL-JUTTNER INITIALIZATION
Te = 0.01
Species(
	name = "eon",
	position_initialization = 'random',
	momentum_initialization = 'mj',
	temperature = [Te,Te,Te],
	particles_per_cell = 300,
	mass = 1.0,
	charge = -1.0,
	number_density = 1.,
	mean_velocity=[0., 0., 0.],
	time_frozen = 10000.,
	boundary_conditions = [
		["periodic", "periodic"],
		["periodic", "periodic"],
	],
	is_test = True
)


for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
	ExternalField(
		field = field,
		profile = gaussian(0.1)
	)


DiagParticleBinning(
 	deposited_quantity = "weight",
 	every = 1000.,
 	species = ["eon"],
 	axes = [
 		["px", -0.4, 0.4, 100]
 	]
)



DiagFields(
	every = 5,
)


