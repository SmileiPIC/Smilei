import math

L0 = 2.*math.pi # Wavelength in PIC units
T = 1. * L0;
Main(
	geometry = "1d3v",
	
	interpolation_order = 2,
	
	timestep = 0.01 * L0,
	sim_time  = T,
	
	cell_length = [0.1 * L0],
	sim_length  = [1. * L0],
	
	number_of_patches = [ 4 ],
	
	time_fields_frozen = 10000000.,
	
	bc_em_type_x  = ["periodic"],
	
	print_every = 10
)

profiles = [
	tconstant   (),
	ttrapezoidal(start=0.2*T, plateau=0.4*T, slope1=0.2*T, slope2=0.1*T),
	tgaussian   (start=0.2*T, duration=0.7*T, fwhm=0.4*T, center=0.4*T, order=2),
	tpolygonal  (points=[0.1*T, 0.2*T, 0.4*T, 0.8*T], values=[1.,0.5,0.8, 0.1]),
	tcosine     (base=0.2, amplitude=1., start=0.2*T, duration=0.5*T, phi=1., freq=4.),
	tpolynomial (t0=0.4*T, order0=-1., order1=-1./T, order2=(3./T)**2)
]

n_profiles = len(profiles)
Main.sim_length[0] *= n_profiles

for i, profile in enumerate(profiles):
	Antenna(
		field = "Jz",
		space_profile = trapezoidal(1., xvacuum=i*L0, xplateau=L0),
		time_profile = profile
	)

DiagFields(
	every = 2,
	fields = ["Jz"]
)


