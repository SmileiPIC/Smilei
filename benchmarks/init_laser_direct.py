# ----------------------------------------------------------------------------------------
#					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
import math
import cmath
import numpy as np
from math import log
import scipy.constants as scc
nx =1600
nr = 100
dx= 0.2    #0.785
dr= 1.5    # 0.785
Lsim = [dx*nx,nr*dr] 
dt =0.18   #0.787366579847066      #3.3440009543674392E-016 #0.18
#laser_FWHM_E = 19.80
W0= 25.  # 31.41592653589793    # 25 #  31.41592653589793     # 4.e-6
#laser_initial_position = 2**0.5*laser_FWHM_E
#a0         = 2.
E0=2.
#focus      = [laser_initial_position]
lambda0= 6.283185307179586 #0.8e-6
theta_pol=0. 
cep_phase=0. 
phi2_chirp=0.
prop_dir=1.
ctau= 19.981310667784744 #39.269908169872416   #5.e-6
z0=Lsim[0]/2.
zf=Lsim[0]/2.
c = scc.c
#zmin= -157.    #-20.e-6
#zmax=157. #20.e-6
#rmin=0.
#rmax=157. # 20.e-6
#z=np.linspace(zmin, zmax,nz,endpoint=False)
#dz = z[1]-z[0]
#r=np.linspace(rmin, rmax, nr, endpoint=False)
#dr = r[1]-r[0]
#r += dr/2
zr=(np.pi*W0**2)/lambda0
inv_zr= 1./zr
k0 = 2*np.pi/lambda0
inv_ctau2 = 1./(ctau)**2
stretch_factor = 1 - 2j * phi2_chirp * c**2 * inv_ctau2
nstep = 201
norder_l=0
norder_r=0
Main(
    geometry = "AMcylindrical",
    number_of_AM = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [dx, dr],
    grid_length  = Lsim,
    number_of_patches = [2, 1 ],
    timestep = dt,
    simulation_time = nstep*dt,
     
    EM_boundary_conditions = [
        ["periodic","periodic"],
        ["buneman","buneman"],
    ],
    
    random_seed = smilei_mpi_rank,
    print_every = 10,
    is_spectral=True,
    uncoupled_grids = True,
    is_pxr = True,
    norder = [0,0],
    pseudo_spectral_guardells = 20
)

################ Laser gaussian pulse, defined through external fields ###################

# Electromagnetic fields of a gaussian beam (fundamental mode), linearly polarized in the y direction
# formulas from B. Quesnel, P. Mora, PHYSICAL REVIEW E 58, no. 3, 1998 
# (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.58.3719)
# version AM created following also A.F.Lifschitz, Journ. Comput. Phys. 228 (2009) 1803–1814
# (doi:10.1016/j.jcp.2008.11.017)


def diffract_factor(x):
    diff=1. + 1j * prop_dir*(x - zf) * inv_zr
    return diff

#( prop_dir*(x- z0)-c*t ) is the right expression but here we define for t=0
def exp_argument(x,r):
    exp = - 1j*cep_phase\
         + 1j*k0*(prop_dir*(x- z0))\
         - (r**2) / (W0**2 * diffract_factor(x)) \
         - 1./stretch_factor*inv_ctau2 * \
         ( prop_dir*(x - z0)  )**2
    return exp

def profil_gauss(x,r):
    profil= np.exp(exp_argument(x,r)) /(diffract_factor(x) * stretch_factor**0.5)
    return profil


"""Class that calculates a Gaussian laser pulse."""
def GaussianLaser(x, r ):
    # Diffraction and stretch_factor
    # Calculate the argument of the complex exponential
    #Get the transverse profile
    # Get the projection along x and y, with the correct polarization
    Er = E0 * profil_gauss(x,r)*np.exp(1j*theta_pol)
    #Et = - 1j *E0 * profil_gauss(x,r)*np.exp(1j*theta_pol)
    return Er

print("stope here")

def Er(x,r):
    return GaussianLaser(x,r )
def Er_mode_1(x,r):
    return np.real(Er(x,r))
def Et_mode_1(x,r):
    return  -1j*Er_mode_1(x,r)

def Br_mode_1(x,r):
    #return  -Et_mode_1(x,r)
    return  1j*Er_mode_1(x,r)
def Bt_mode_1(x,r):
    return  Er_mode_1(x,r)


field_profile = { 'Er_mode_1': Er_mode_1, 'Et_mode_1': Et_mode_1, 'Br_mode_1': Br_mode_1, 'Bt_mode_1': Bt_mode_1}

for field in ['Er_mode_1', 'Et_mode_1', 'Br_mode_1', 'Bt_mode_1']:
        ExternalField(
                field = field,
                profile = field_profile[field],
        )

DiagFields(
    every = 10,
    fields = ["Br_m_mode_0", "Br_m_mode_1","Bl_m_mode_0","Bl_m_mode_1","Bt_m_mode_0","Bt_m_mode_1","Bt_mode_0","Bt_mode_1","Bl_mode_0","Bl_mode_1","Br_mode_0","Br_mode_1","Er_mode_0","Er_mode_1","Et_mode_0","Et_mode_1","El_mode_0","El_mode_1" ]
)
