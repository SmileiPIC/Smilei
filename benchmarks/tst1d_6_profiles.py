import math

L0 = 2.*math.pi # Wavelength in PIC units

Main(
	geometry = "1d3v",
	
	interpolation_order = 2,
	
	timestep = 0.005 * L0,
	sim_time  = 0.01 * L0,
	
	cell_length = [0.01 * L0],
	sim_length  = [1. * L0],
	
	number_of_patches = [ 4 ],
	
	time_fields_frozen = 10000000.,
	
	bc_em_type_x  = ["periodic"],
	
	print_every = 10
)

profiles = {
"constant"   :constant   (1.),
"trapezoidal":trapezoidal(1., xvacuum=0.1*L0, xplateau=0.4*L0, xslope1=0.1*L0, xslope2=0.1*L0),
"gaussian"   :gaussian   (1., xvacuum=0.1*L0, xlength =0.5*L0, xfwhm=0.2*L0, xcenter=0.2*L0, xorder=2),
"polygonal"  :polygonal  (xpoints=[0.1*L0, 0.2*L0, 0.4*L0, 0.8*L0], xvalues=[1.,0.5,0.8, 0.1]),
"cosine"     :cosine     (1., xamplitude=0.4, xvacuum=0.3*L0, xlength=0.4*L0, xphi=0.1*L0, xnumber=3),
"polynomial" :polynomial (x0=0.4*L0, order0=1., order1=-1./L0, order2=(1./L0)**2)
}

for name, profile in profiles.items():
	Species(
		species_type = name,
		initPosition_type = "regular",
		initMomentum_type = "maxwell-juettner",
		n_part_per_cell= 1000,
		mass = 1.0,
		charge = 1.0,
		nb_density = profile,
		time_frozen = 10000.0,
		bc_part_type_xmin = "none",
		bc_part_type_xmax = "none"
	)


DiagFields(
	every = 5,
)


