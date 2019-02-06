import math
import cmath
from numpy import exp, sqrt, arctan, vectorize, real, sin, cos, arctan, zeros_like, arange
from scipy import integrate,zeros_like 
from math import log
import h5py

dx = 0.393
dtrans = 4.712
dt = 0.96 * dx
nx =  600 
ntrans = 100*2 
Lx = nx * dx
Ltrans = ntrans * dtrans
npatch_x = 8 
npatch_trans = 8
Nit =  1



Main(
    geometry = "3Dcartesian",
    interpolation_order = 2 ,
    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],
    number_of_patches = [ npatch_x, npatch_trans, npatch_trans ],
    timestep = dt,
    simulation_time = dt*Nit,
    clrw = 5,
    EM_boundary_conditions = [ ['silver-muller'] ],
    random_seed = smilei_mpi_rank
)

#LaserGaussian3D(
#    a0              = 1.,
#    omega           = 1.,
#    focus           = [2.*10, 5.*l0, 5.*l0],
#    waist           = 10,
#    incidence_angle = [0., 0.],
#    time_envelope   = tgaussian(center=2*10., fwhm=10.)
#)

################ Laser gaussian pulse, defined through external fields ###################

# Electromagnetic fields of a gaussian beam (fundamental mode), linearly polarized in the y direction
# formulas from B. Quesnel, P. Mora, PHYSICAL REVIEW E 58, no. 3, 1998 
# (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.58.3719)

a0    = 0.01
laser_FWHM_E = 58.81
focus = [ Main.grid_length[0]/2., Main.grid_length[1]/2.,Main.grid_length[2]/2.]
laser_initial_position = Main.grid_length[0]/2.

c_vacuum = 1.
waist       = 157.
omega       = 1.
Zr          = omega * waist**2/2.  # Rayleigh length



# time gaussian function
def time_gaussian(fwhm, center, order=2):
    import math
    sigma = (0.5*fwhm)**order/log(2.0)
    def f(t):
        return exp( -(t-center)**order / sigma )
    return f

time_envelope_t              = time_gaussian(center=laser_initial_position                  , fwhm=laser_FWHM_E)
time_envelope_t_plus_half_dt = time_gaussian(center=(laser_initial_position+c_vacuum*0.5*dt), fwhm=laser_FWHM_E)

# laser waist function
def w(x):
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        return w

def space_envelope(x,y,z):
        coeff = omega * (x-focus[0]) * w(x)**2 / (2.*Zr**2)
        invWaist2 = (w(x)/waist)**2
        spatial_amplitude = w(x) * exp( -invWaist2*(  (y-focus[1])**2 + (z-focus[2])**2 )  )
        phase = coeff * ( (y-focus[1])**2 + (z-focus[2])**2 )
        return a0 * spatial_amplitude * exp( 1j*( phase-arctan((x-focus[0])/ Zr)  )  )

def complex_exponential_comoving(x,t):
        csi = x-c_vacuum*t-laser_initial_position # comoving coordinate
        return exp(1j*csi)

### Electromagnetic field
# Electric field        
#def Ex(x,y,z):
#        invWaist2 = (w(x)/waist)**2
#        complexEx = 2.* (y-focus[1]) * invWaist2 * space_envelope(x,y,z) * complex_exponential_comoving(x,0.)
#        return real(complexEx)*time_envelope_t(x)

# Define Ex as solution of Poisson
def Ex(x,y,z):
    A = zeros_like(x)
    nx, ny, nz = x.shape
    ix = round((x[0,0,0]+5*dx/2.)/dx) 
    iy = round((y[0,0,0]+2*dtrans)/dtrans )
    iz = round((z[0,0,0]+2*dtrans)/dtrans )
    h5f = h5py.File('tabulated.h5','r') 
    B = h5f['dataset_1'][ix:ix+nx,iy:iy+ny,iz:iz+nz]
    print(ix, nx, iy, ny, iz, nz)
    h5f.close()
    return A


def Ey(x,y,z):
        complexEy  = 1j * space_envelope(x,y,z) * complex_exponential_comoving(x,0)
        return real(complexEy)*time_envelope_t(x)


def Ez(x,y,z):
        return 0.*x

# Magnetic field
def Bx(x,y,z):
        invWaist2 = (w(x)/waist)**2
        complexBx = 2.* (z-focus[2]) * invWaist2 * space_envelope(x,y,z) * complex_exponential_comoving(x,dt/2.)
        return real(complexBx)*time_envelope_t_plus_half_dt(x)

def By(x,y,z):
        return 0.*x

def Bz(x,y,z):
        complexBz = 1j * space_envelope(x,y,z) * complex_exponential_comoving(x,dt/2.)
        return real(complexBz)*time_envelope_t_plus_half_dt(x)

field_profile = {'Ex': Ex, 'Ey': Ey, 'Ez': Ez, 'Bx': Bx, 'By': By, 'Bz': Bz}

for field in ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz']:
        ExternalField(
                field = field,
                profile = field_profile[field],
        )


##########################################################################################

globalEvery = int(100)

DiagScalar(
    every=globalEvery
)

#DiagFields(
#    every = globalEvery,
#    fields = ['Ex','Ey','Ez']
#)
#from numpy import s_
#DiagFields(
#    every = globalEvery,
#    fields = ['Ex','Ey','Ez'],
#    subgrid = s_[4:100:3, 5:400:10, 6:300:80]
#)

DiagProbe(
    every = 10,
    origin = [0.1*Main.grid_length[0], 0.5*Main.grid_length[1], 0.5*Main.grid_length[2]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [nx],
    origin = [0.1*Main.grid_length[0], 0.5*Main.grid_length[1], 0.5*Main.grid_length[2]],
    corners = [[0.9*Main.grid_length[0], 0.5*Main.grid_length[1], 0.5*Main.grid_length[2]]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [nx, ntrans],
    origin = [0.*Main.grid_length[0], 0.*Main.grid_length[1], 0.5*Main.grid_length[2]],
    corners = [
        [1.*Main.grid_length[0], 0. *Main.grid_length[1], 0.5*Main.grid_length[2]],
        [0.*Main.grid_length[0], 1.*Main.grid_length[1], 0.5*Main.grid_length[2]],
    ],
    fields = []
)

