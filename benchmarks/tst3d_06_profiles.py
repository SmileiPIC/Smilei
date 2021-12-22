import math

L0 = 2.*math.pi # Wavelength in PIC units

Main(
	geometry = "3Dcartesian",
	
	interpolation_order = 2,
	
	timestep = 0.005 * L0,
	simulation_time  = 0.02 * L0,
	
	cell_length = [0.01 * L0]*3,
	grid_length  = [1. * L0]*3,
	
	number_of_patches = [ 4 ]*3,
	
	time_fields_frozen = 0.01,
	
	EM_boundary_conditions = [ ["periodic"] ],
	print_every = 1,
	solve_poisson = False,
	
)

def custom(x, y, z):
	if z<L0/2.: return 0.
	else: return math.exp(-z)

profiles = {
"constant"   :constant   (1.),
"trapezoidal":trapezoidal(1.,
							xvacuum=0.1*L0, xplateau=0.4*L0, xslope1=0.1*L0, xslope2=0.1*L0,
							yvacuum=0.3*L0, yplateau=0.1*L0, yslope1=0.2*L0, yslope2=0.2*L0,
							zvacuum=0.2*L0, zplateau=0.6*L0, zslope1=0.1*L0, zslope2=0.1*L0),
"gaussian"   :gaussian   (1.,
							xvacuum=0.1 *L0, xlength=0.5*L0, xfwhm=0.2*L0, xcenter=0.2*L0, xorder=2,
							yvacuum=0.2 *L0, ylength=0.6*L0, yfwhm=0.4*L0, ycenter=0.5*L0, yorder=4,
							zvacuum=0.05*L0, zlength=0.7*L0, zfwhm=0.6*L0, zcenter=0.4*L0, zorder=2),
"polygonal"  :polygonal  (xpoints=[0.1*L0, 0.2*L0, 0.4*L0, 0.8*L0], xvalues=[1.,0.5,0.8, 0.1]),
"cosine"     :cosine     (1.,
							xamplitude=0.4, xvacuum=0.3*L0, xlength=0.4*L0, xphi=0.1*L0, xnumber=5,
							yamplitude=0.2, yvacuum=0.2*L0, ylength=0.6*L0, yphi=0.3*L0, ynumber=10,
							zamplitude=0.1, zvacuum=0.1*L0, zlength=0.6*L0, zphi=0. *L0, znumber=2),
"polynomial" :polynomial (x0=0.4*L0, y0=0.5*L0, z0=0.6*L0,
							order0=1.,
							order1=[-1./L0,-0.1/L0, 2./L0],
							order2=[(0.1/L0)**2]*6),
"custom"     :custom
}

for name, profile in profiles.items():
	Species(
		name = name,
		position_initialization = "regular",
		momentum_initialization = "maxwell-juettner",
		particles_per_cell= 8,
		mass = 1.0,
		charge = 1.0,
		number_density = profile,
		time_frozen = 10000.0,
		boundary_conditions = [
			["periodic", "periodic"],
			["periodic", "periodic"],
			["periodic", "periodic"],
		],
	)


for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
	ExternalField(
		field = field,
		profile = gaussian(10.)
	)
	PrescribedField(
		field = field,
		profile = lambda x,y,z,t: 10.*np.cos(x+y+z)*np.sin(math.pi/2*t/Main.simulation_time)
	)


DiagFields(
	every = 4,
)


