# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
import math
import cmath
from numpy import exp, sqrt, arctan, vectorize, real
from math import log

dx = 0.125
dr = 1.5
nx = 960*2
nr = 200
Lsim = [dx*nx,nr*dr] # length of the simulation
dt = dx # spectral solver feature

Main(
    geometry = "AMcylindrical",
    number_of_AM = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [dx, dr],
    grid_length  = Lsim,
    number_of_patches = [ 64, 4 ],
    timestep = dt,
    simulation_time = 1*dt,
     
    EM_boundary_conditions = [
        ["zero","zero"],
        #["periodic","periodic"],
        #["silver-muller","silver-muller"],
        ["silver-muller","buneman"],
    ],
    
    random_seed = smilei_mpi_rank,
    print_every = 10,
    is_spectral=True,
    uncoupled_grids = True,
    is_pxr = True,
    norder = [32,0],
    pseudo_spectral_guardells = 64,
    apply_rotational_cleaning = False,

)

#MovingWindow(
#    time_start = 0.,
#    velocity_x = 1.
#)

ne = 0.0045
begin_upramp = Main.grid_length[0]/2.  #Density is 0 before that and up ramp starts.
Lupramp = 0.2*Main.grid_length[0]  #Length of the plateau 
Lplateau = Lupramp  #Length of the plateau 
Ldownramp = 0. #Length of the down ramp
xplateau = begin_upramp + Lupramp # Start of the plateau
begin_downramp = xplateau + Lplateau # Beginning of the output ramp. 
finish = begin_downramp + Ldownramp # End of plasma

g = polygonal(xpoints=[begin_upramp, xplateau, begin_downramp, finish], xvalues=[0, ne, ne, 0.])

def my_profile(x,y):
    if( y < 0.75*Main.grid_length[1]):
        limiter = 1.
    else:
        limiter = 0.
    return limiter*g(x,y)

Species( 
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    ionization_model = "none",
    particles_per_cell = 4*4*16,
    regular_number = [4,4,16],
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = my_profile,  # Here absolute value of the charge is 1 so charge_density = nb_density
    mean_velocity = [0., 0., 0.],
    time_frozen = 0.0,
    boundary_conditions = [
    	["remove","remove"],
    	["remove","remove"],
    ],
)

Species( 
    name = "ion",
    position_initialization = "electron",
    momentum_initialization = "cold",
    ionization_model = "none",
    particles_per_cell = 4*4*16,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = my_profile,  # Here absolute value of the charge is 1 so charge_density = nb_density
    mean_velocity = [0., 0., 0.],
    time_frozen = 10000000.0,
    boundary_conditions = [
    	["remove","remove"],
    	["remove","remove"],
    ],
)



################ Laser gaussian pulse, defined through external fields ###################

# Electromagnetic fields of a gaussian beam (fundamental mode), linearly polarized in the y direction
# formulas from B. Quesnel, P. Mora, PHYSICAL REVIEW E 58, no. 3, 1998 
# (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.58.3719)
# version AM created following also A.F.Lifschitz, Journ. Comput. Phys. 228 (2009) 1803–1814
# (doi:10.1016/j.jcp.2008.11.017)

laser_FWHM_E = 19.80*2.
waist      = 25.
laser_initial_position = Main.grid_length[0]/4. #2**0.5*laser_FWHM_E
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

def time_sin2(fwhm, center):
    import scipy
    import math
    slope1=fwhm
    slope2=fwhm
    start=center-fwhm
    def f(t):
        w = scipy.zeros_like(t)
        w[t < start+slope1+slope2] = scipy.cos(0.5*math.pi*(t[t < start+slope1+slope2]-start-slope1)/slope2)**2
        w[t < start+slope1]        = scipy.sin(0.5*math.pi*(t[t < start+slope1]-start)/slope1)**2 
        w[t<start] = 0.
        return w
    return f

time_envelope_t              = time_sin2(center=laser_initial_position                  , fwhm=laser_FWHM_E)

# laser waist function
def w(x):
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        #w[x>Main.grid_length[0]/2.] = 0.
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
        #complexBz = 1j * space_envelope(x,r) * complex_exponential_comoving(x,dt/2.)
        complexBz = 1j * space_envelope(x,r) * complex_exponential_comoving(x, 0.)
        return real(complexBz)*time_envelope_t(x)
        
def Br_mode_1(x,r):
        return 1j*Bz(x,r)
        
def Bt_mode_1(x,r):
        return Bz(x,r)

def Bl_mode_1(x,r):
        invWaist2 = (w(x)/waist)**2
        #complexBx = 2.* r * invWaist2 * space_envelope(x,r) * complex_exponential_comoving(x,dt/2.)
        complexBx = 2.* r * invWaist2 * space_envelope(x,r) * complex_exponential_comoving(x,0.)
        # to obtain Bx, this Bx will be multiplied by sin(theta) in the AM reconstruction [z = r*sin(theta)]
        Bx = real(complexBx)*time_envelope_t(x)
        return 1j*Bx
                
field_profile = {'El_mode_1': El_mode_1, 'Er_mode_1': Er_mode_1, 'Et_mode_1': Et_mode_1, 'Bl_mode_1': Bl_mode_1, 'Br_mode_1': Br_mode_1, 'Bt_mode_1': Bt_mode_1}

for field in ['El_mode_1', 'Er_mode_1', 'Et_mode_1', 'Bl_mode_1', 'Br_mode_1', 'Bt_mode_1']:
        ExternalField(
                field = field,
                profile = field_profile[field],
        )

DiagFields(
    every = 50,
    fields = ["Br","Br_m","Bl_m","Bt_m","Bl","Bt","Et","El", "Jr_electron","Jt_electron","Rho_electron","Jl_electron","Er","Rho","Rho_ion","Jr_ion","Jl_ion"]
    #fields = ["Br","Br_m","Bl_m","Bt_m","Bl","Bt","Er","Et","El"]
    #fields = ["Br_m_mode_0", "Br_m_mode_1","Bl_m_mode_0","Bl_m_mode_1","Bt_m_mode_0","Bt_m_mode_1","Bt","Bl_mode_0","Bl_mode_1","Br_mode_0","Br_mode_1","Er_mode_0","Er_mode_1","Et_mode_0","Et_mode_1","El_mode_0","El_mode_1","Rho_mode_0", "Rho_mode_1", "Jr_mode_0","Jr_mode_1","Jt_mode_0","Jt_mode_1","Jl_mode_0","Jl_mode_1","Rho_electron_mode_0","Rho_electron_mode_1"]
)

#DiagProbe(
#    every = 10,
#    origin = [1., 10., 0.],
#    fields = []
#)
#DiagProbe(
#    every = 200,
#    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx", "Jy", "Jz", "Rho", "Rho_electron_mode_0","Jr_electron_mode_0","Jr_electron_mode_1","Jr_ion_mode_0","Jr_ion_mode_1"],
#    origin = [0., -Lsim[1], 0.],
#    corners = [[Lsim[0], -Lsim[1], 0.], [0., Lsim[1], 0.]],
#    number=[nx, 2*nr],
#)
#DiagProbe(
#    every = 10,
#    origin = [0., -10., 0.],
#    corners = [[Lsim[0], -10., 0.]],
#    number=[100],
#    fields = []
#)
#DiagTrackParticles(
#    species = "ion",
#    every = 100,
#    attributes = ["x", "px", "py", "weight", "q"]
#)
#DiagTrackParticles(
#    species = "electron",
#    every = 100,
#    attributes = ["x", "px", "py", "weight", "q"]
#)

