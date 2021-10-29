import math

L0 = 2.*math.pi # Wavelength in PIC units

dx = 0.01 * L0
Lsim = 1. * L0

# Make a profile in a file
def preprocess():
	if smilei_mpi_rank == 0:
		from numpy import linspace, ceil, exp
		import h5py
		x = linspace(0., Lsim, int(ceil(Lsim/dx)) + 1 )
		with h5py.File("test_profile.h5","w") as f:
			g = f.create_group("some_group")
			g.create_dataset("the_profile", data = exp(-((x-0.4*L0)/(0.1*L0))**2) )

Main(
	geometry = "1Dcartesian",
	
	interpolation_order = 2,
	
	timestep = 0.005 * L0,
	simulation_time  = 1. * L0,
	
	cell_length = [dx],
	grid_length  = [Lsim],
	
	number_of_patches = [ 4 ],
	
	time_fields_frozen = 10000000.,
	
	EM_boundary_conditions = [ ["periodic"] ],
	
	print_every = 10,
	solve_poisson = False,
	
)


# SPATIAL PROFILES
def custom(x):
	if x<L0/2.: return 0.
	else: return math.exp(-x)
profiles = {
"constant"   :constant   (1.),
"trapezoidal":trapezoidal(1., xvacuum=0.1*L0, xplateau=0.4*L0, xslope1=0.1*L0, xslope2=0.3*L0),
"gaussian"   :gaussian   (1., xvacuum=0.1*L0, xlength =0.5*L0, xfwhm=0.2*L0, xcenter=0.2*L0, xorder=2),
"polygonal"  :polygonal  (xpoints=[0.1*L0, 0.2*L0, 0.4*L0, 0.8*L0], xvalues=[1.,0.5,0.8, 0.1]),
"cosine"     :cosine     (1., xamplitude=0.4, xvacuum=0.3*L0, xlength=0.4*L0, xphi=0.1*L0, xnumber=3),
"polynomial" :polynomial (x0=0.4*L0, order0=1., order1=-1./L0, order2=(1./L0)**2),
"custom"     :custom,
"file"       :"test_profile.h5/some_group/the_profile"
}
for name, profile in profiles.items():
	Species(
		name = name,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= 1000,
		mass = 1.0,
		charge = 1.0,
		number_density = profile,
		time_frozen = 10000.0,
		boundary_conditions = [
			["periodic", "periodic"],
		],
		temperature = [0.1]
	)

# TEMPORAL PROFILES
def tcustom(t):
	if t<L0/2.: return 0.
	else: return math.exp(-t)
tprofiles = [
["tconstant"   ,tconstant   ()],
["ttrapezoidal",ttrapezoidal(start=0.2*L0, plateau=0.4*L0, slope1=0.2*L0, slope2=0.1*L0)],
["tgaussian"   ,tgaussian   (start=0.2*L0, duration=0.7*L0, fwhm=0.4*L0, center=0.4*L0, order=2)],
["tpolygonal"  ,tpolygonal  (points=[0.1*L0, 0.2*L0, 0.4*L0, 0.8*L0], values=[1.,0.5,0.8, 0.1])],
["tcosine"     ,tcosine     (base=0.2, amplitude=1., start=0.2*L0, duration=0.5*L0, phi=1., freq=4.)],
["tpolynomial" ,tpolynomial (t0=0.4*L0, order0=-1., order1=-1./L0, order2=(3./L0)**2)],
["tcustom"     ,tcustom]
]
antenna_width = Main.grid_length[0]/len(tprofiles)
for i, (name, tprofile) in enumerate(tprofiles):
	Antenna(
		field = "Jz",
		space_profile = trapezoidal(1., xvacuum=i*antenna_width, xplateau=antenna_width),
		time_profile = tprofile
	)

# RELATIVISTIC MAXWELL-JUTTNER INITIALIZATION
Te = 1.
g0 = 5.
mj_species = ['eon_nodrift', 'eon_xdrift', 'eon_ydrift', 'eon_zdrift']
for eon in mj_species:
	vmean = [0., 0., 0.]
	direction = eon[4]
	if direction in "xyz":
		vmean[ "xyz".index(direction) ] = math.sqrt(1.-g0**-2)
		pmax = 100.
	else:
		direction = "x"
		pmax = 10.
	
	Species(
		name = eon,
		position_initialization = 'random',
		momentum_initialization = 'mj',
		temperature = [Te,Te,Te],
		particles_per_cell = 100000,
		mass = 1.0,
		charge = -1.0,
		number_density = 1.,
		mean_velocity=vmean,
		time_frozen = 10000.,
		boundary_conditions = [["periodic"]],
		is_test = True
	)
	
	DiagParticleBinning(
	 	deposited_quantity = "weight",
	 	every = 1000.,
	 	species = [eon],
	 	axes = [
	 		["p"+direction, -10., pmax, 500]
	 	]
	)

for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
	ExternalField(
		field = field,
		profile = gaussian(0.1)
	)

DiagFields(
	every = 2,
	fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jz"] + ["Rho_"+name for name in profiles.keys()],
)
