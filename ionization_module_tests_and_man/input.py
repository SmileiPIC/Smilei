import numpy as np
# import sys
# sys.path.append('/home/ent/Smilei_ionization/simulations/')


l0 = 2.0*np.pi                # laser wavelength
t0 = l0                       # optical cicle
Lsim = [12.*l0]               # length of the simulation
Tsim = 24.*t0                 # duration of the simulation               
resx = 128.                    # nb of cells in one laser wavelength
dx = dy = dz = l0/resx
cell_length  = [dx]
dt  = 0.95*dx/np.sqrt(1.)
rest = t0/dt
timesteps = int(Tsim/dt)

Z = 18
mass = 39.948


# laser pulse input parameters
a0 = 100.
omega = 1.
N_cycles = 10.
xf = Lsim[0]/2.
tcenter = N_cycles*l0

Lmu  = 0.8                              # mu-m

c = 299792458.                          # m/sec
omega_ref = 2*np.pi*c/(Lmu*1.e-6)       # 1/sec
eps0 = 8.8541878128e-12                 # F⋅m^−1
charge_SI = 1.60217663e-19              # C
m_e_SI = 9.1093837e-31                  # kg
N_r = eps0*m_e_SI*omega_ref**2/charge_SI**2
print('Reference density = ',N_r*10**(-6), ' cm-3')

I_S = 4.65e29                #W/cm^2
l_C =  3.8615901e-13         #m
I_ref = 0.5*(a0*omega_ref*l_C/c)**2.*4.65e29
print('Reference intensity = ',I_ref, ' W/cm^2')


# n0 = 5.e13/(N_r*1.e-6)  # initial density 5.e13 cm^(-3)
# nppc = 16                               # nb of particle per cells


n0 = 5.e13/(N_r*1.e-6)
thickness = 0.25*l0
position = Lsim[0]/2.-thickness/2.
def nf(x):
    return n0*np.heaviside(x-position, 1.)*np.heaviside(position+thickness-x, 1.)
nppc = 128

Main(
    geometry = "1Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = cell_length,
    grid_length  = Lsim,
    
    number_of_patches = [1],
    
    timestep = dt,
    simulation_time = Tsim,

    reference_angular_frequency_SI = omega_ref,
    
    EM_boundary_conditions = [ ["silver-muller", "silver-muller"] ],
    
)

def sin2_envelope(t, tcenter, N):
    if (-np.pi*N <= t-tcenter) and (t-tcenter <= np.pi*N):
        return np.sin((t-tcenter+np.pi*N)/2./N)**2
    else:
        return 0.



Laser(
    box_side = "xmin", 
    omega = omega,
    space_time_profile = [ lambda t: 0., 
                           lambda t: a0*sin2_envelope(t, tcenter, N_cycles)\
                                       *np.cos(omega*((t-tcenter)-xf)) ],  
)


Species(
    name = 'atom_tunnel',
    ionization_model = 'tunnel',
    ionization_electrons = 'electron',
    atomic_number = Z,
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = mass*1836.0,
    charge = 0.,
    # number_density = n0,
    number_density = nf,
    boundary_conditions = [['remove','remove']]
)


Species(
    name = 'atom_TL',
    ionization_model = 'tunnel_TL',
    ionization_tl_parameter = 6,
    ionization_electrons = 'electron',
    atomic_number = Z,
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = mass*1836.0,
    charge = 0.,
    # number_density = n0,
    number_density = nf,
    boundary_conditions = [['remove','remove']]
)

Species(
    name = 'atom_BSI',
    ionization_model = 'tunnel_BSI',
    ionization_electrons = 'electron',
    atomic_number = Z,
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = mass*1836.0,
    charge = 0.,
    # number_density = n0,
    number_density = nf,
    boundary_conditions = [['remove','remove']]
)

Species(
    name = 'atom_full_PPT',
    ionization_model = 'tunnel_full_PPT',
    ionization_electrons = 'electron',
    atomic_number = Z,
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = mass*1836.0,
    charge = 0.,
    # number_density = n0,
    number_density = nf,
    boundary_conditions = [['remove','remove']]
)

Species(
    name = 'electron',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = 0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.0,
    boundary_conditions = [['remove','remove']],
#   time_frozen = 2.*Tsim
)


# DiagScalar(
#     every = rest
# )


DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    species = ["atom_tunnel"],
    axes = [
        ["charge", -0.5, Z+0.5, Z+1]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    species = ["atom_TL"],
    axes = [
        ["charge", -0.5, Z+0.5, Z+1]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    species = ["atom_BSI"],
    axes = [
        ["charge", -0.5, Z+0.5, Z+1]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 1,
    species = ["atom_full_PPT"],
    axes = [
        ["charge", -0.5, Z+0.5, Z+1]
    ]
)

DiagFields(
    every = 100,
    fields = ['Ex','Ey','Ez','Bx','By','Bz']
)