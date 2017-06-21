# _____________________________________________________________________________
#
# Electron trajectory in a plane wave 
# with a Gaussian temporal profile.
#
# Validation in the relativist regime
# 
# _____________________________________________________________________________

import math

# _____________________________________________________________________________
# Main parameters

l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle
Lx = 50*l0

n0 = 1e-8                     # particle density

Tsim = 120.*t0                 # duration of the simulation
resx = 64.                    # nb of cells in one laser wavelength

dx = l0/resx                            # space step
dt  = 0.95 * dx                 		# timestep (0.95 x CFL)

a0 = 5
start = 0                               # Laser start
fwhm = 10*t0                            # Gaussian time fwhm
duration = 90*t0                        # Laser duration
center = duration*0.5                   # Laser profile center

pusher_list = ["norm","vay","higueracary"]  # dynamic type

# Density profile for inital location of the particles
def n0_(x):
        if (dx<x<2*dx):
                return n0
        else:
                return 0.

# ______________________________________________________________________________
# Namelists

Main(
    geometry = "1d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [dx],
    sim_length  = [Lx],
    
    number_of_patches = [32],
    
    timestep = dt,
    sim_time = Tsim,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    
    random_seed = 0
)

LaserPlanar1D(
    box_side         = "xmin",
    a0              = a0,
    omega           = 1.,
    polarization_phi = 0.,
    ellipticity     = 1,
    time_envelope  = tgaussian(start=start,duration=duration,fwhm=fwhm,center=center,order=2)
)

for ipusher,pusher in enumerate(pusher_list):
    Species(
        species_type = "electron_" + pusher,
        position_initialization = "centered",
        momentum_initialization = "cold",
        n_part_per_cell = 10,
        c_part_max = 1.0,
        mass = 1.0,
        charge = -1.0,
        charge_density = n0_,
        mean_velocity = [0., 0.0, 0.0],
        temperature = [0.],
        dynamics_type = pusher,
	    boundary_conditions = [["periodic"]],
        track_every = 10,
        track_flush_every = 100,
        is_test = True
    )


