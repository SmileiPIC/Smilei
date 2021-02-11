# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
import math
import cmath
from numpy import exp, sqrt, arctan, vectorize, real, zeros_like
from math import log

dx = 0.251327
dt = dx # spectral solver feature
dtrans = 1.96349
nx =  960
ntrans = 256
Lx = nx * dx
Ltrans = ntrans * dtrans
npatch_x = 64
npatch_trans =32
Nit = 500

Main(
    geometry = "AMcylindrical",
    number_of_AM=2,
    interpolation_order = 1,
    timestep = dt,
    simulation_time = dt*Nit,
    cell_length  = [dx, dtrans],
    grid_length = [ Lx,  Ltrans],
    number_of_patches = [npatch_x, npatch_trans],
    #clrw = 5,
    EM_boundary_conditions = [
        ["zero","zero"],
        ["silver-muller","buneman"],
    ],
    random_seed = smilei_mpi_rank,
    solve_poisson = False,
    print_every = 100,
    is_spectral=True,
    uncoupled_grids = True,
    is_pxr = True,
    norder = [32,0],
    pseudo_spectral_guardells = 64,
    apply_rotational_cleaning = True,
    number_of_damping_cells = [44],

)

MovingWindow(
    time_start = 0., #Leaves 2 patches untouched, in front of the laser.
    velocity_x = 0.996995486
)

ne = 0.0045
begin_upramp = Main.grid_length[0] + 10.  #Density is 0 before that and up ramp starts.
Lupramp = 100. #Length of the upramp 
Lplateau = 15707.  #Length of the plateau 
Ldownramp = 2356.19 #Length of the down ramp
xplateau = begin_upramp + Lupramp # Start of the plateau
begin_downramp = xplateau + Lplateau # Beginning of the output ramp. 
finish = begin_downramp + Ldownramp # End of plasma

g = polygonal(xpoints=[begin_upramp, xplateau, begin_downramp, finish], xvalues=[0, ne, ne, 0.])

def my_profile(x,y):
    return g(x,y)

Species( 
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    ionization_model = "none",
    particles_per_cell = 30,
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
    particles_per_cell = 30,
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

laser_FWHM_E = 41.
waist      = 120.
laser_initial_position = Main.grid_length[0]/2. 
a0         = 2.
focus      = [begin_upramp]


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

time_envelope_t = time_gaussian(center=laser_initial_position, fwhm=laser_FWHM_E)

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
                
field_profile = {'Er_mode_1': Er_mode_1, 'Et_mode_1': Et_mode_1}

for field in ['Er_mode_1', 'Et_mode_1']:
        ExternalField(
                field = field,
                profile = field_profile[field],
        )

DiagProbe(
	every = [500,500],
	origin = [0., -Ltrans/6., Ltrans/6.],
	corners = [
              [Main.grid_length[0], -Ltrans/6., Ltrans/6.],
                  ],
	number = [nx],
        fields = ["Ey","Jy"],
)

DiagFields(
        every = 500,
)

DiagPerformances(
	every = 1000,
)
