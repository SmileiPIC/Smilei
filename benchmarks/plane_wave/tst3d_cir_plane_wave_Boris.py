# _____________________________________________________________________________
#
# Electron trajectory in a plane wave with Gaussian 
# temporal profile.
#
# _____________________________________________________________________________

import math

# _____________________________________________________________________________
# Main parameters

l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle
Lx = 12*l0
Ly = 2.*l0
Lz = 2.*l0

n0 = 1e-8                     # particle density

Tsim = 90.*t0                 # duration of the simulation
resx = 40.                    # nb of cells in one laser wavelength

dx = l0/resx                            # space step
dt  = 0.95 * dx/math.sqrt(3.)		# timestep (0.95 x CFL)

start = 0                               # Laser start
fwhm = 10*t0                            # Gaussian time fwhm
duration = 90*t0                        # Laser duration
center = duration*0.5                   # Laser profile center

pusher = "boris"                         # dynamic type

# Density profile for inital location of the particles
def n0_(x,y,z):
        if (dx<x<2*dx) and (0.5*Ly-dx<y<0.5*Ly+dx) and (0.5*Lz-dx<z<0.5*Lz+dx):
                return n0
        else:
                return 0.

# ______________________________________________________________________________
# Namelists

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [dx,dx,dx],
    grid_length  = [Lx,Ly,Lz],
    
    number_of_patches = [ 4,4,4 ],
    
    timestep = dt,
    simulation_time = Tsim,
    
    EM_boundary_conditions = [
        ['silver-muller'],
        ['periodic'],
        ['periodic'],
    ],
    
    random_seed = 0
)

#Laser(
#    box_side        = "xmin",
#    omega          = 1.,
#    chirp_profile  = tconstant(),
#    time_envelope  = tgaussian(start=0,duration=90*t0,fwhm=30*t0,center=45*t0,order=2),
#    space_envelope = [ 1. , 1. ],
#    phase          = [ 0., 0. ]
#)

LaserGaussian3D(
    box_side         = "xmin",
    a0              = 2.,
    omega           = 1.,
    focus           = [0., 0.5*Ly, 0.5*Lz],
    waist           = 1e9,
    incidence_angle = [0., 0.],
    polarization_phi = 0.,
    ellipticity     = 1,
    time_envelope  = tgaussian(start=start,duration=duration,fwhm=fwhm,center=center,order=2)
)

Species(
    name = "electron",
    position_initialization = "centered",
    momentum_initialization = "cold",
    particles_per_cell = 10,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_,
    mean_velocity = [0., 0.0, 0.0],
    temperature = [0.],
    pusher = pusher,
    boundary_conditions = [
    	["remove", "remove"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
    is_test = False
)

DiagScalar(
    every = 10,
    vars=['Uelm','Ukin','Utot','Uexp','Ubal']
)

DiagFields(
    every = 100,
    fields = ['Ex','Ey','Ez','By','Bz']
)

DiagTrackParticles(
    species = "electron",
    every = 2,
    flush_every = 100,
)
