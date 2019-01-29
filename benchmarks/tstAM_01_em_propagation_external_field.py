# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
import math
import cmath
from numpy import exp, sqrt, arctan, vectorize, real
from math import log

dx = 0.2
dr = 1.5
nx = 800
nr = 100
Lsim = [dx*nx,nr*dr] # length of the simulation
dt = 0.18

Main(
    geometry = "AMcylindrical",
    number_of_AM = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [dx, dr],
    grid_length  = Lsim,
    number_of_patches = [ 32, 4 ],
    timestep = dt,
    simulation_time = 2.5*nx*dt,
     
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    
    random_seed = smilei_mpi_rank,
    print_every = 100,
)

#MovingWindow(
#    time_start = Main.grid_length[0]-50*dx,
#    velocity_x = 0.9997
#)




################ Laser gaussian pulse, defined through external fields ###################

# Electromagnetic fields of a gaussian beam (fundamental mode), linearly polarized in the y direction
# formulas from B. Quesnel, P. Mora, PHYSICAL REVIEW E 58, no. 3, 1998 
# (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.58.3719)
# version AM created following also A.F.Lifschitz, Journ. Comput. Phys. 228 (2009) 1803–1814
# (doi:10.1016/j.jcp.2008.11.017)

laser_FWHM_E = 19.80
waist      = 25.
laser_initial_position = 2**0.5*laser_FWHM_E
a0         = 2.
focus      = [laser_initial_position]


c_vacuum = 1.
dx = Main.cell_length[0]
dr = Main.cell_length[1]
dt = Main.timestep


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

def space_envelope(x,r):
        coeff = omega * (x-focus[0]) * w(x)**2 / (2.*Zr**2)
        invWaist2 = (w(x)/waist)**2
        spatial_amplitude = w(x) * exp( -invWaist2*( r**2  ) )
        phase = coeff * (r**2) 
        return a0 * spatial_amplitude * exp( 1j*( phase-arctan((x-focus[0])/ Zr)  )  )

def complex_exponential_comoving(x,t):
        csi = x-c_vacuum*t-laser_initial_position # comoving coordinate
        return exp(1j*csi)

### Electromagnetic field
# Electric field        
def Ey(x,r):
        complexEy = 1j * space_envelope(x,r) * complex_exponential_comoving(x,0)
        return real(complexEy)*time_envelope_t(x)
        
def Er_mode_1(x,r):
        return Ey(x,r)
        
def Et_mode_1(x,r):
        return -1j*Ey(x,r)
        
def El_mode_1(x,r):
        invWaist2 = (w(x)/waist)**2
        complexEx = 2.* r * invWaist2 * space_envelope(x,r) * complex_exponential_comoving(x,0.)
        # to obtain Ex, this Ex will be multiplied by cos(theta) in the AM reconstruction [y = r*cos(theta)]
        Ex = real(complexEx)*time_envelope_t(x)
        return Ex

# Magnetic field
def Bz(x,r):
        complexBz = 1j * space_envelope(x,r) * complex_exponential_comoving(x,dt/2.)
        return real(complexBz)*time_envelope_t_plus_half_dt(x)
        
def Br_mode_1(x,r):
        return 1j*Bz(x,r)
        
def Bt_mode_1(x,r):
        return Bz(x,r)

def Bl_mode_1(x,r):
        invWaist2 = (w(x)/waist)**2
        complexBx = 2.* r * invWaist2 * space_envelope(x,r) * complex_exponential_comoving(x,dt/2.)
        # to obtain Bx, this Bx will be multiplied by sin(theta) in the AM reconstruction [z = r*sin(theta)]
        Bx = real(complexBx)*time_envelope_t_plus_half_dt(x)
        return 1j*Bx
                
field_profile = {'El_mode_1': El_mode_1, 'Er_mode_1': Er_mode_1, 'Et_mode_1': Et_mode_1, 'Bl_mode_1': Bl_mode_1, 'Br_mode_1': Br_mode_1, 'Bt_mode_1': Bt_mode_1}

for field in ['El_mode_1', 'Er_mode_1', 'Et_mode_1', 'Bl_mode_1', 'Br_mode_1', 'Bt_mode_1']:
        ExternalField(
                field = field,
                profile = field_profile[field],
        )

DiagFields(
    every = 100,
    fields = ["Br_m_mode_0", "Br_m_mode_1","Bl_m_mode_0","Bl_m_mode_1","Bt_m_mode_0","Bt_m_mode_1","Bt_mode_0","Bt_mode_1","Bl_mode_0","Bl_mode_1","Br_mode_0","Br_mode_1","Er_mode_0","Er_mode_1","Et_mode_0","Et_mode_1","El_mode_0","El_mode_1" ]
)

#DiagProbe(
#    every = 10,
#    origin = [1., 10., 0.],
#    fields = []
#)
#DiagProbe(
#    every = 10,
#    origin = [0., 10., 0.],
#    corners = [[Lsim[0], 10., 0.]],
#    number=[100],
#    fields = []
#)
#DiagProbe(
#    every = 10,
#    origin = [0., -10., 0.],
#    corners = [[Lsim[0], -10., 0.]],
#    number=[100],
#    fields = []
#)


