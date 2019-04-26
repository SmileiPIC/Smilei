# ----------------------------------------------------------------------------------------
#					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
import math
import cmath
from numpy import exp, sqrt, arctan, vectorize, real
from math import log
import scipy.constants as scc
nz = 400
nr = 200
dt = 3.3440009543674392E-016 #0.18
laser_FWHM_E = 19.80
waist      = 4.e-6
laser_initial_position = 2**0.5*laser_FWHM_E
#a0         = 2.
focus      = [laser_initial_position]

ctau=5.e-6
z0=0.
zf=0.
c = scc.c
zmin=-20.e-6
zmax=20.e-6
rmin=0.
rmax=20.e-6
z=np.linspace(zmin, zmax,Nz,endpoint=False)
dz = z[1]-z[0]
r=np.linspace(rmin, rmax, Nr, endpoint=False)
dr = r[1]-r[0]
r += dr/2
Lsim = [dz*nz,nr*dr] # length of the simulation
zr=(pi*W0**2)/lambda0
nstep = 50
Main(
    geometry = "AMcylindrical",
    number_of_AM = 2,
    interpolation_order = 2 ,
    solve_poisson = False,
    cell_length = [dz, dr],
    grid_length  = Lsim,
    number_of_patches = [ 32, 4 ],
    timestep = dt,
    simulation_time = nstep*dt,
     
    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],
    
    random_seed = smilei_mpi_rank,
    print_every = 100,
)

################ Laser gaussian pulse, defined through external fields ###################

# Electromagnetic fields of a gaussian beam (fundamental mode), linearly polarized in the y direction
# formulas from B. Quesnel, P. Mora, PHYSICAL REVIEW E 58, no. 3, 1998 
# (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.58.3719)
# version AM created following also A.F.Lifschitz, Journ. Comput. Phys. 228 (2009) 1803–1814
# (doi:10.1016/j.jcp.2008.11.017)


"""Class that calculates a Gaussian laser pulse."""

    def GaussianLaser(  a0, waist, ctau, z0, zf=None, theta_pol=0.,
    lambda0=0.8e-6, cep_phase=0., phi2_chirp=0.,
    prop_dir=1, r, z ):
        """
        Define a linearly-polarized Gaussian laser profile.
	More precisely, the electric field **near the focal plane**
	is given by:
	.. math::
	E(\\boldsymbol{x},t) = a_0\\times E_0\,
	\exp\left( -\\frac{r^2}{w_0^2} - \\frac{(z-z_0-ct)^2}{c^2\\tau^2} \\right)
	\cos[ k_0( z - z_0 - ct ) - \phi_{cep} ]
	where :math:`k_0 = 2\pi/\\lambda_0` is the wavevector and where
	:math:`E_0 = m_e c^2 k_0 / q_e` is the field amplitude for :math:`a_0=1`.
	.. note::
	The additional terms that arise **far from the focal plane**
	(Gouy phase, wavefront curvature, ...) are not included in the above
	formula for simplicity, but are of course taken into account by
	the code, when initializing the laser pulse away from the focal plane.
	Parameters
	----------
	a0: float (dimensionless)
	The peak normalized vector potential at the focal plane, defined
	as :math:`a_0` in the above formula.
	waist: float (in meter)
	Laser waist at the focal plane, defined as :math:`w_0` in the
	above formula.
	tau: float (in second)
	The duration of the laser (in the lab frame),
	defined as :math:`\\tau` in the above formula.
	z0: float (in meter)
	The initial position of the centroid of the laser
	(in the lab frame), defined as :math:`z_0` in the above formula.
	zf: float (in meter), optional
	The position of the focal plane (in the lab frame).

	"""

	k0 = 2*np.pi/lambda0
	E0 = 1. #a0*m_e*c**2*k0/e
	zr = 0.5*k0*waist**2

	# If no focal plane position is given, use z0
	if zf is None:
	zf = z0

	# Store the parameters
	inv_zr = 1./zr
	w0 = waist
	inv_ctau2 = 1./(ctau)**2

	# Note: this formula is expressed with complex numbers for compactness
	# and simplicity, but only the real part is used in the end
	# (see final return statement)
	# The formula for the laser (in complex numbers) is obtained by
	# multiplying the Fourier transform of the laser at focus
	# E(k_x,k_y,\omega) = exp( -(\omega-\omega_0)^2(\tau^2/4 + \phi^(2)/2)
	# - (k_x^2 + k_y^2)w_0^2/4 ) by the paraxial propagator
	# e^(-i(\omega/c - (k_x^2 +k_y^2)/2k0)(z-z_foc))
	# and then by taking the inverse Fourier transform in x, y, and t

	# Diffraction and stretch_factor
	diffract_factor = 1. + 1j * prop_dir*(z - zf) * inv_zr
	stretch_factor = 1 - 2j * phi2_chirp * c**2 * inv_ctau2
	# Calculate the argument of the complex exponential
	exp_argument = - 1j*cep_phase \
	+ 1j*k0*( prop_dir*(z - z0) - c*t ) \
	- (r**2) / (w0**2 * diffract_factor) \
	- 1./stretch_factor*inv_ctau2 * \
	( prop_dir*(z - z0) - c*t )**2
	# Get the transverse profile
	profile = np.exp(exp_argument) /(diffract_factor * stretch_factor**0.5)


	# Get the projection along x and y, with the correct polarization
	Er = E0 * profile*exp(1j*theta_pol)
	Et = - 1j *E0 * profile*exp(1j*theta_pol)

        return( Er, Et )


Er_mode_1=Er
Et_mode_1= Et

Br_mode_1= -Et_mode_1/clight
Bt_mode_1= Er_mode_1/clight


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
