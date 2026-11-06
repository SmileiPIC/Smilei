# Vacuum propagation of a Gaussian laser envelope injected from xmin in AMcylindrical geometry,
# derived from the benchmark tstAM_00_em_propagation.py.
#
# import_from_lasy = False : the envelope is injected with LaserEnvelopeGaussianAM
# import_from_lasy = True  : the same laser is first written to a file in the LASY format
#                            (openPMD, thetaMode geometry), without importing lasy,
#                            then injected with LaserEnvelopeFromLasyAM.
#                            The file contains, depending on lasy_envelope_field,
#                            "normalized_vector_potential" : the complex envelope of e*A/(m_e*c), used directly
#                            "electric_field"              : the complex envelope of E in V/m, from which
#                                                            LaserEnvelopeFromLasyAM computes the envelope of A
# The runs should give the same envelope (with the electric field, within the accuracy of the
# conversion from the envelope of E to the envelope of A, second order in the timestep of the LASY file),
# so this namelist can be used as a benchmark of the LASY interface for the envelope model.

import math
import numpy as np
import scipy.constants as cst

import_from_lasy        = True
lasy_envelope_field     = "normalized_vector_potential" # or "electric_field"
lasy_filename           = "lasy_gaussian_envelope_00000.h5"

# Reference angular frequency (rad/s), needed to convert the LASY file to/from SI units
laser_wavelength_SI     = 0.8e-6
reference_frequency_SI  = 2. * math.pi * cst.c / laser_wavelength_SI
# Laser angular frequency (rad/s) and the same in Smilei units
omega0_SI               = reference_frequency_SI
omega                   = omega0_SI / reference_frequency_SI

# Grid (same as tstAM_00_em_propagation.py)
cell_length_x           = 0.2
cell_length_r           = 1.5
number_of_cells_x       = 800
number_of_cells_r       = 100
grid_length             = [cell_length_x * number_of_cells_x, cell_length_r * number_of_cells_r]
timestep                = 0.18

Main(
    geometry                       = "AMcylindrical",
    # the envelope model in AMcylindrical geometry uses only the azimuthal mode 0
    number_of_AM                   = 1,
    interpolation_order            = 2 ,
    solve_poisson                  = False,
    cell_length                    = [cell_length_x, cell_length_r],
    grid_length                    = grid_length,
    number_of_patches              = [ 32, 4 ],
    timestep                       = timestep,
    simulation_time                = 2.5 * number_of_cells_x * timestep,
    reference_angular_frequency_SI = reference_frequency_SI,

    EM_boundary_conditions = [
        ["silver-muller","silver-muller"],
        ["buneman","buneman"],
    ],

    print_every = 100,
)

MovingWindow(
    time_start = Main.grid_length[0] - 50 * cell_length_x,
    velocity_x = 0.9997
)

# Laser parameters (same as tstAM_00_em_propagation.py), in Smilei units
a0                      = 2.
waist                   = 25.
focus_x                 = 0.
laser_fwhm              = 19.80
laser_center            = 2**0.5 * laser_fwhm

envelope_solver               = "explicit"
envelope_boundary_conditions  = [ ["reflective", "reflective"], ["reflective", "reflective"], ]

if not import_from_lasy:

    LaserEnvelopeGaussianAM(
        a0                           = a0,
        omega                        = omega,
        focus                        = [focus_x],
        waist                        = waist,
        time_envelope                = tgaussian(center=laser_center, fwhm=laser_fwhm),
        envelope_solver              = envelope_solver,
        Envelope_boundary_conditions = envelope_boundary_conditions,
        box_side                     = "xmin",
        polarization_phi             = 0.,
        ellipticity                  = 0.
    )

else:

    # Time window written in the LASY file: the laser field is ~1e-11 of its peak after it
    lasy_time_end       = laser_center + 3. * laser_fwhm
    # The LASY grid is finer than the Smilei one, and contains the Smilei sampling points
    # (t = (i+1/2)*timestep, r = j*cell_length_r/2) of LaserEnvelopeFromLasyAM
    lasy_dt             = timestep / 2.
    lasy_dr             = cell_length_r / 4.
    lasy_number_of_t    = int(math.ceil(lasy_time_end / lasy_dt)) + 1
    lasy_number_of_r    = int(math.ceil((grid_length[1] + 3. * cell_length_r) / lasy_dr)) + 1
    # Number of timesteps during which the envelope is injected
    number_of_injection_timesteps = int(math.ceil(lasy_time_end / timestep))

    def write_lasy_file(filename):
        # Writes the laser of LaserEnvelopeGaussianAM (focus, waist, a0, tgaussian time envelope,
        # box_side = "xmin") in the LASY format: E = Re[ E_envelope * exp(-i*omega0*t) ] * polarization,
        # with the envelope stored as azimuthal modes on a (t, r) grid, in SI units.
        # Only the m = 0 mode is non zero, since the envelope has azimuthal symmetry.
        import h5py

        # Normalized time and radius of the LASY grid
        lasy_t                    = np.arange(lasy_number_of_t) * lasy_dt
        lasy_r                    = np.arange(lasy_number_of_r) * lasy_dr
        t_grid, r_grid            = np.meshgrid(lasy_t, lasy_r, indexing="ij")

        # Envelope of the vector potential at x = 0, same formulas as LaserEnvelopeGaussianAM
        rayleigh_length           = omega * waist**2 / 2.
        waist_ratio               = math.sqrt(1. / (1. + (focus_x / rayleigh_length)**2))
        one_ov_waist_at_plane_2   = (waist_ratio / waist)**2
        curvature_phase_coeff     = -omega * focus_x * waist_ratio**2 / (2. * rayleigh_length**2)
        gouy_phase                = math.atan(-focus_x / rayleigh_length)
        space_envelope_A          = (a0 * omega * waist_ratio * np.exp(-one_ov_waist_at_plane_2 * r_grid**2)
                                     * np.exp(1j * (curvature_phase_coeff * r_grid**2 - gouy_phase)))

        # Same time envelope as tgaussian, and its time derivative
        tgaussian_sigma           = (0.5 * laser_fwhm)**2 / math.log(2.)
        time_envelope             = np.exp(-(t_grid - laser_center)**2 / tgaussian_sigma)
        time_envelope_derivative  = -2. * (t_grid - laser_center) / tgaussian_sigma * time_envelope

        if lasy_envelope_field == "normalized_vector_potential":

            # Envelope of the normalized vector potential, dimensionless
            envelope_modes          = (space_envelope_A * time_envelope)[np.newaxis, :, :].astype(np.complex128)
            unit_dimension          = np.zeros(7)

        else:

            # Envelope of E from the one of A: E_envelope = -(d/dt - i*omega) A_envelope
            reference_electric_field_SI = cst.m_e * cst.c * reference_frequency_SI / cst.e
            envelope_E_SI           = (space_envelope_A * (1j * omega * time_envelope - time_envelope_derivative)
                                       * reference_electric_field_SI)
            envelope_modes          = envelope_E_SI[np.newaxis, :, :].astype(np.complex128)
            # Electric field: kg m s^-3 A^-1
            unit_dimension          = np.array([1., 1., -3., -1., 0., 0., 0.])

        with h5py.File(filename, "w") as f:
            # openPMD attributes, as written by LASY (and checked by FBPIC's reader)
            f.attrs["openPMD"]            = np.bytes_("1.1.0")
            f.attrs["openPMDextension"]   = np.uint32(0)
            f.attrs["basePath"]           = np.bytes_("/data/%T/")
            f.attrs["meshesPath"]         = np.bytes_("meshes/")
            f.attrs["iterationEncoding"]  = np.bytes_("fileBased")
            f.attrs["iterationFormat"]    = np.bytes_("lasy_gaussian_envelope_%T.h5")
            f.attrs["software"]           = np.bytes_("lasy")
            f.attrs["softwareVersion"]    = np.bytes_("0.6.0")
            f.attrs["comment"]            = np.bytes_("Written by a Smilei namelist in the LASY format, without lasy")

            iteration_group               = f.create_group("data/0")
            iteration_group.attrs["dt"]         = 1.
            iteration_group.attrs["time"]       = 0.
            iteration_group.attrs["timeUnitSI"] = 1.

            meshes_group                  = f.create_group("data/0/meshes")
            dset                          = meshes_group.create_dataset("laserEnvelope", data=envelope_modes)
            dset.attrs["geometry"]           = np.bytes_("thetaMode")
            dset.attrs["geometryParameters"] = np.bytes_("m=1;imag=+")
            dset.attrs["axisLabels"]         = np.array([b"t", b"r"])
            dset.attrs["dataOrder"]          = np.bytes_("C")
            dset.attrs["gridSpacing"]        = np.array([lasy_dt / reference_frequency_SI,
                                                         lasy_dr * cst.c / reference_frequency_SI])
            dset.attrs["gridGlobalOffset"]   = np.array([0., 0.])
            dset.attrs["gridUnitSI"]         = 1.
            dset.attrs["position"]           = np.array([0., 0.])
            dset.attrs["unitSI"]             = 1.
            dset.attrs["unitDimension"]      = unit_dimension
            dset.attrs["envelopeField"]      = np.bytes_(lasy_envelope_field)
            dset.attrs["timeOffset"]         = 0.
            dset.attrs["angularFrequency"]   = omega0_SI
            # E along y in Smilei, i.e. along x in LASY (linear polarization, polarization_phi = 0)
            dset.attrs["polarization"]       = np.array([1. + 0j, 0. + 0j])

    # The namelist is executed by every MPI process: only one writes the file,
    # and all wait for it before reading it in LaserEnvelopeFromLasyAM
    try:
        from mpi4py import MPI
        mpi_rank                  = MPI.COMM_WORLD.Get_rank()
        mpi_barrier               = MPI.COMM_WORLD.Barrier
    except ImportError:
        mpi_rank                  = 0
        mpi_barrier               = lambda: None

    if mpi_rank == 0:
        write_lasy_file(lasy_filename)
    mpi_barrier()

    # The reference frequency is the one of the Main block, omega0_SI is checked against the file
    LaserEnvelopeFromLasyAM(
        lasy_filename,
        number_of_injection_timesteps,
        omega0_SI                    = omega0_SI,
        envelope_solver              = envelope_solver,
        Envelope_boundary_conditions = envelope_boundary_conditions
    )

DiagFields(
    every = 100,
    fields = ["Env_A_abs_mode_0", "Env_E_abs_mode_0"]
)

DiagProbe(
    every   = 100,
    origin  = [0., 2. * cell_length_r, 0.],
    corners = [ [Main.grid_length[0], 2. * cell_length_r, 0.] ],
    number  = [number_of_cells_x],
    fields  = ["Env_A_abs", "Env_E_abs"]
)