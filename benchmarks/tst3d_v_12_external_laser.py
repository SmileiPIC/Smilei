# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
import numpy
import math
import cmath

lambda0 = 0.8e-6               # m
c = 299792458                  # m/s
omega0 = 2*math.pi*c/lambda0   # rad/s
c_ov_omega0 = c/omega0         # m/rad
eps0=8.854187e-12              # Vacuum permittivity, F/m
e=1.60217646e-19               # Elementary charge, C
me=9.10938215e-31              # Electron mass, kg
ncrit = eps0*omega0**2*me/e**2 # Plasma critical number density [m-3]
normalized_electron_charge = -1 # For electrons


dx = lambda0 / 32. / c_ov_omega0
dtrans = lambda0   / c_ov_omega0
dt = dx * 0.95
ncells_x = 1792
Lx = ncells_x * dx
ntrans = 128
Ltrans = ntrans * dtrans
npatch_x = 128
npatch_trans = 16
Nit = 4000

Main(
    geometry = "3Dcartesian",
    interpolation_order = 2,
    timestep = dt,
    simulation_time = dt*Nit,
    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],
    number_of_patches = [npatch_x, npatch_trans, npatch_trans],
    clrw = 14,
    EM_boundary_conditions = [ ["silver-muller"] ],
    EM_boundary_conditions_k = [ [1., 0., 0.],[-1., 0., 0.],[1., 0.005, 0.],[1., -0.005, 0.],[1., 0., 0.005],[1., 0., -0.005] ],
    random_seed = 0,
    solve_poisson = False,
    solve_relativistic_poisson = True,
    print_every = 100,
    save_magnectic_fields_for_SM = False,
)

Vectorization(
    mode = "adaptive",
    reconfigure_every = 20,
    initial_mode = "off"
)

MovingWindow(
    time_start = 0, #Laser and bunch are correctly positionned right from the start
    velocity_x = 0.9999569
)

LoadBalancing(
    initial_balance = False,
    every = [int(ncells_x/dt), 20], # Start balancing only when the box is full of plasma
    cell_load = 2.,
)


# LASER -------------------------------------------------------------------------
laser_FWHM_I             = 30e-15 * omega0 
laser_FWHM_E             = laser_FWHM_I * math.sqrt(2) 
laser_focus_x          = Main.grid_length[0]
laser_initial_position = Main.grid_length[0] - 2 * laser_FWHM_E
waist                  = 20.0e-6/c_ov_omega0
a0                     = math.sqrt(2) 
##  BACKGROUND ELECTRON ------------------------------------------
#
n0 = 1.5e23 / ncrit 

Species( 
    name = "electron_background",
    position_initialization = "random",
    momentum_initialization = "cold",
    ionization_model = "none",
    particles_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0,  # Here absolute value of the charge is 1 so charge_density = nb_density
    mean_velocity = [0.0, 0.0, 0.0],
    pusher = 'boris',
    time_frozen = 0.0,
    boundary_conditions = [
    	["remove", "remove"],
    	["remove", "remove"],
    	["remove", "remove"],
    ],
)


################################# Laser field, from external fields

focus = [laser_focus_x, Main.grid_length[1]/2., Main.grid_length[2]/2.]

c_vacuum = 1.
dx = Main.cell_length[0]
dy = Main.cell_length[1]
dz = Main.cell_length[2]
dt = Main.timestep

omega       = 1.
Zr          = omega * waist**2/2.  # Rayleigh length

# time gaussian function
def time_gaussian(fwhm, center, order=2):
    import math
    import numpy as np
    sigma = (0.5*fwhm)**order/math.log(2.0)

    def f(t):
        return np.exp( -(t-center)**order / sigma )

    return f


time_envelope_t              = time_gaussian(center=laser_initial_position                  , fwhm=laser_FWHM_E)
time_envelope_t_plus_half_dt = time_gaussian(center=(laser_initial_position+c_vacuum*0.5*dt), fwhm=laser_FWHM_E)

# laser waist function
def w(x):
        import numpy as np
        w  = np.sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        return w

def coeff(x):
        import numpy as np
        coeff = omega * (x-focus[0]) * w(x)**2 / (2.*Zr**2)
        return coeff

def spatial_amplitude(x,y,z):
        import numpy as np
        invWaist2 = (w(x)/waist)**2
        return w(x) * np.exp( -invWaist2*(  (y-focus[1])**2 + (z-focus[2])**2 )  )

# laser phase   
def phase(x,y,z):
        import numpy as np

# laser waist function
def w(x):
        import numpy as np
        w  = np.sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        return w

def coeff(x):
        import numpy as np
        coeff = omega * (x-focus[0]) * w(x)**2 / (2.*Zr**2)
        return coeff

def spatial_amplitude(x,y,z):
        import numpy as np
        invWaist2 = (w(x)/waist)**2
        return w(x) * np.exp( -invWaist2*(  (y-focus[1])**2 + (z-focus[2])**2 )  )

# laser phase   
def phase(x,y,z):
        import numpy as np
        return coeff(x) * ( (y-focus[1])**2 + (z-focus[2])**2 )

def Gouy_phase(x):
        import numpy as np
        return np.arctan(   (x-focus[0]) / Zr     )

def space_envelope(x,y,z):
        import numpy as np
        return a0 * spatial_amplitude(x,y,z) * np.exp(1j*phase(x,y,z)) * np.exp(-1j*Gouy_phase(x))

def complex_exponential_comoving(x,t):
        import numpy as np
        csi = x-c_vacuum*t-laser_initial_position # comoving coordinate
        return np.exp(1j*csi)

### Electromagnetic field
# Electric field        
def Ex(x,y,z):
        import numpy as np
        invWaist2 = (w(x)/waist)**2
        complexEx = 2.* (y-focus[1]) * invWaist2 * space_envelope(x,y,z) * complex_exponential_comoving(x,0.)
        return np.multiply(np.real(complexEx),time_envelope_t(x))

def Ey(x,y,z):
        import numpy as np
        complexEy  = 1j * space_envelope(x,y,z) * complex_exponential_comoving(x,0)
        return np.multiply(np.real(complexEy),time_envelope_t(x))


def Ez(x,y,z):
        import numpy as np
        return np.zeros(shape=np.shape(x))

# Magnetic field
def Bx(x,y,z):
        import numpy as np
        invWaist2 = (w(x)/waist)**2
        complexBx = 2.* (z-focus[2]) * invWaist2 * space_envelope(x,y,z) * complex_exponential_comoving(x,dt/2.)
        return np.multiply(np.real(complexBx),time_envelope_t_plus_half_dt(x))

def By(x,y,z):
        import numpy as np
        return np.zeros(shape=np.shape(x))

def Bz(x,y,z):
        import numpy as np
        complexBz = 1j * space_envelope(x,y,z) * complex_exponential_comoving(x,dt/2.)
        return np.multiply(np.real(complexBz),time_envelope_t_plus_half_dt(x))



field_profile = {'Ex': Ex, 'Ey': Ey, 'Ez': Ez, 'Bx': Bx, 'By': By, 'Bz': Bz}

for field in ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz']:
        ExternalField(
                field = field,
                profile = field_profile[field],
        )


###############################################################################


DiagProbe(
	every = 200,
	origin = [0., Main.grid_length[1]/2., Main.grid_length[2]/2.],
	corners = [
              [Main.grid_length[0], Main.grid_length[1]/2., Main.grid_length[2]/2.]
                  ],
	number = [ncells_x],
	fields = ['Ex','Ey','Rho']
)

DiagProbe(
	every = 1000,
	origin = [0., Main.grid_length[1]/4., Main.grid_length[2]/2.],
	corners =  [
           [Main.grid_length[0], Main.grid_length[1]/4., Main.grid_length[2]/2.],
	   [0., 3.*Main.grid_length[1]/4., Main.grid_length[2]/2.],
                   ],
	number = [ncells_x/4,ntrans/2],
	fields = ['Ex','Ey','Ez','Bx','By','Bz','Rho']
)

DiagProbe(
	every = 5000,
	origin = [0., 0., Main.grid_length[2]/2.],
	corners = [
           [Main.grid_length[0], 0., Main.grid_length[2]/2.],
	   [0., Main.grid_length[1], Main.grid_length[2]/2.],
                  ],
	number = [ncells_x,ntrans],
	fields = ['Ex','Ey','Ez','Rho','Jx']
)

DiagScalar(every = 200, vars=['Uelm','ExMax','ExMaxCell','EyMax','EyMaxCell', 'RhoMax', 'RhoMaxCell'])

DiagPerformances(
    every = 500,
)


