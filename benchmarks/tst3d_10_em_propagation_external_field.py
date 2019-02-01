import math
import cmath
from numpy import exp, sqrt, arctan, vectorize, real, sin, cos, arctan
from scipy import integrate, zeros_like
from math import log

l0 = 2.0*math.pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [7.*l0,10.*l0,10.*l0]  # length of the simulation
Tsim = 90.                # duration of the simulation
resx = 16.                    # nb of cells in one laser wavelength
rest = 30.                    # nb of timesteps in one optical cycle 

dt = t0/rest

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4,4,4 ],
    
    timestep = dt,
    simulation_time = Tsim,
    
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

a0    = 1.
laser_FWHM_E = 10.
focus = [ 2.*laser_FWHM_E, 5.*l0, 5.*l0 ]
laser_initial_position = 2.*laser_FWHM_E

c_vacuum = 1.
dx = Main.cell_length[0]
dy = Main.cell_length[1]
dz = Main.cell_length[2]
dt = Main.timestep

waist       = 10.
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
    integration_constant = 0.
    nx = x.shape[0]
    ny = x.shape[1]
    nz = x.shape[2]
    A = zeros_like(x)
    y2d = y[0,:,:]-focus[1]
    z2d = z[0,:,:]-focus[2]
    xpatch = x[:,0,0]
    rsquare2d = y2d**2+z2d**2
    dx = Main.cell_length[0]
    xprior = scipy.arange(dx/2.,xpatch[0],dx/2.) #Ex is dual along x

    def r2overRC(xp):
        if xp != 0. :
            return -0.5*rsquared2d/(xp + Zr**2/xp)
        else:
            return 0.
    def yoverRC(xp):
        if xp != 0. :
            return y2d/(xp + Zr**2/xp)
        else:
            return 0.
    def invWaist2(xp):
        return  (w(xp)/waist)**2
    def spatial_amplitude(xp):
        return   w(xp) * exp( -invWaist2(xp)*rsquared2d)
    def modified_phase(xp):
        return xp + r2overRC(xp) + arctan(xp/Zr)
    def integrand_hyperslab(xp):
        return a0*time_envelope_t(xp)*spatial_amplitude(xp)*( 2*y2d*invWaist2(xp)*sin(modified_phase(xp)) + yoverRC(xp) * cos(modified_phase(xp) ))

    for xp in xprior:
        A[0,:,:] += integrand_hyperslab(xp)
    A[0,:,:] = A[0,:,:]*dx  #This is the trapezoid method with homogeneous dx.  

    for ix in range(nx-1):
        slab = integrand_hyperslab(xpatch[ix])*dx/2.
        A[ix,:,:] += slab
        A[ix+1,:,:] = A[ix,:,:] + slab
    A[-1,:,:] += integrand_hyperslab(xpatch[nx-1])*dx/2.
    
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

globalEvery = int(rest)

DiagScalar(
    every=globalEvery
)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez']
)
from numpy import s_
DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez'],
    subgrid = s_[4:100:3, 5:400:10, 6:300:80]
)

DiagProbe(
    every = 10,
    origin = [0.1*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [30],
    origin = [0.1*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    corners = [[0.9*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [10, 10],
    origin = [0.1*Lsim[0], 0.*Lsim[1], 0.5*Lsim[2]],
    corners = [
        [0.9*Lsim[0], 0. *Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0.9*Lsim[1], 0.5*Lsim[2]],
    ],
    fields = []
)

DiagProbe(
    every = 100,
    number = [4, 4, 4],
    origin = [0.1*Lsim[0], 0.*Lsim[1], 0.5*Lsim[2]],
    corners = [
        [0.9*Lsim[0], 0. *Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0.9*Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0. *Lsim[1], 0.9*Lsim[2]],
    ],
    fields = []
)
