# Vacuum propagation of a Gaussian laser in AMcylindrical geometry,
# derived from the benchmark tstAM_00_em_propagation.py.
#
# import_from_lasy = False : the laser is injected with LaserGaussianAM, as in tstAM_00
# import_from_lasy = True  : the same laser is first written to a file in the LASY format
#                            (openPMD, thetaMode geometry), without importing lasy,
#                            then injected with LaserFromLasyAM.
#                            The file contains, depending on lasy_envelope_field,
#                            "electric_field"              : the complex envelope of E in V/m
#                            "normalized_vector_potential" : the complex envelope of e*A/(m_e*c),
#                                                            converted to E by LaserFromLasyAM
# The runs should give the same fields, so this namelist can be used as a benchmark
# of the LASY interface.

import math
import numpy as np
import scipy.constants as cst

import_from_lasy          = True
lasy_envelope_field       = "electric_field" # or "normalized_vector_potential"
lasy_filename             = "lasy_gaussian_laser_00000.h5"
# write_benchmark_lasy_file = True  : the namelist writes the benchmark Gaussian laser to lasy_filename
# write_benchmark_lasy_file = False : lasy_filename is an existing LASY file, which is only read
# To use your own LASY file, set lasy_filename to its path AND write_benchmark_lasy_file = False.
# The namelist never overwrites a file it did not write: it stops with an error instead.
write_benchmark_lasy_file = True

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
    number_of_AM                   = 2,
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

if not import_from_lasy:

    LaserGaussianAM(
        box_side        = "xmin",
        a0              = a0,
        omega           = omega,
        focus           = [focus_x],
        waist           = waist,
        time_envelope   = tgaussian(center=laser_center, fwhm=laser_fwhm)
    )

else:

    # Time window written in the LASY file: the laser field is ~1e-11 of its peak after it
    lasy_time_end       = laser_center + 3. * laser_fwhm
    # The LASY grid is finer than the Smilei one, and contains the Smilei sampling points
    # (t = (i+1/2)*timestep, r = j*cell_length_r/2) of LaserFromLasyAM
    lasy_dt             = timestep / 2.
    lasy_dr             = cell_length_r / 4.
    lasy_number_of_t    = int(math.ceil(lasy_time_end / lasy_dt)) + 1
    lasy_number_of_r    = int(math.ceil((grid_length[1] + 3. * cell_length_r) / lasy_dr)) + 1
    benchmark_lasy_file_comment = "Written by a Smilei namelist in the LASY format, without lasy"

    def write_lasy_file(filename):
        # Writes the laser of LaserGaussianAM (focus, waist, a0, tgaussian time envelope)
        # in the LASY format: E = Re[ envelope * exp(-i*omega0*t) ] * polarization,
        # with the envelope stored as azimuthal modes on a (t, r) grid, in SI units.
        # Only the m = 0 mode is non zero, since the envelope has azimuthal symmetry.
        import h5py

        # Normalized time and radius of the LASY grid
        lasy_t                    = np.arange(lasy_number_of_t) * lasy_dt
        lasy_r                    = np.arange(lasy_number_of_r) * lasy_dr
        t_grid, r_grid            = np.meshgrid(lasy_t, lasy_r, indexing="ij")

        # Gaussian beam at the injection plane, same formulas as LaserGaussianAM
        rayleigh_length           = omega * waist**2 / 2.
        waist_ratio               = math.sqrt(1. / (1. + (focus_x / rayleigh_length)**2))
        one_ov_waist_at_plane_2   = (waist_ratio / waist)**2
        curvature_phase_coeff     = -omega * focus_x * waist_ratio**2 / (2. * rayleigh_length**2)
        space_envelope            = a0 * omega * waist_ratio * np.exp(-one_ov_waist_at_plane_2 * r_grid**2)
        phase                     = curvature_phase_coeff * r_grid**2

        # Same time envelope as tgaussian, delayed by phase/omega as in Smilei's separable laser
        tgaussian_sigma           = (0.5 * laser_fwhm)**2 / math.log(2.)
        time_envelope             = np.exp(-(t_grid - phase / omega - laser_center)**2 / tgaussian_sigma)

        # Smilei's separable laser gives E_y = space_envelope * time_envelope * sin(omega*t - phase),
        # which is Re[ 1j * exp(1j*phase) * space_envelope * time_envelope * exp(-1j*omega*t) ]
        reference_electric_field_SI = cst.m_e * cst.c * reference_frequency_SI / cst.e

        if lasy_envelope_field == "electric_field":

            envelope_SI           = (1j * np.exp(1j * phase) * space_envelope * time_envelope
                                     * reference_electric_field_SI)
            envelope_modes        = envelope_SI[np.newaxis, :, :].astype(np.complex128)
            # Electric field: kg m s^-3 A^-1
            unit_dimension        = np.array([1., 1., -3., -1., 0., 0., 0.])

        else:

            # Envelope of the normalized vector potential giving the same electric field,
            # integrating in time dA_envelope/dt = i*omega*A_envelope - E_envelope
            # from A_envelope = 0 at t = 0 (the reference laser is zero before), on a time grid
            # refinement_factor times finer than the LASY one. Each fine step is integrated
            # exactly for E_envelope linear inside the step (exponential integrator).
            phase_along_r          = curvature_phase_coeff * lasy_r**2
            space_envelope_along_r = a0 * omega * waist_ratio * np.exp(-one_ov_waist_at_plane_2 * lasy_r**2)
            def normalized_E_envelope(t):
                return (1j * np.exp(1j * phase_along_r) * space_envelope_along_r
                        * np.exp(-(t - phase_along_r / omega - laser_center)**2 / tgaussian_sigma))

            refinement_factor      = 8
            fine_dt                = lasy_dt / refinement_factor
            exp_i_omega_fine_dt    = np.exp(1j * omega * fine_dt)
            # integral over the fine step of exp(i*omega*(fine_dt-s)) ds
            integral_constant_part = (exp_i_omega_fine_dt - 1.) / (1j * omega)
            # integral over the fine step of exp(i*omega*(fine_dt-s)) * s/fine_dt ds
            integral_linear_part   = (integral_constant_part - exp_i_omega_fine_dt / (1j * omega)
                                      + (exp_i_omega_fine_dt - 1.) / ((1j * omega)**2 * fine_dt))

            A_envelope             = np.zeros((lasy_number_of_t, lasy_number_of_r), dtype=np.complex128)
            A_envelope_fine        = np.zeros(lasy_number_of_r, dtype=np.complex128)
            E_envelope_fine_old    = normalized_E_envelope(0.)
            for it in range(lasy_number_of_t - 1):
                for i_fine in range(refinement_factor):
                    E_envelope_fine_new = normalized_E_envelope((it * refinement_factor + i_fine + 1) * fine_dt)
                    A_envelope_fine     = (exp_i_omega_fine_dt * A_envelope_fine
                                           - E_envelope_fine_old * integral_constant_part
                                           - (E_envelope_fine_new - E_envelope_fine_old) * integral_linear_part)
                    E_envelope_fine_old = E_envelope_fine_new
                A_envelope[it+1, :]     = A_envelope_fine

            envelope_modes        = A_envelope[np.newaxis, :, :]
            # Normalized vector potential: dimensionless
            unit_dimension        = np.zeros(7)

        with h5py.File(filename, "w") as f:
            # openPMD attributes, as written by LASY (and checked by FBPIC's reader)
            f.attrs["openPMD"]            = np.bytes_("1.1.0")
            f.attrs["openPMDextension"]   = np.uint32(0)
            f.attrs["basePath"]           = np.bytes_("/data/%T/")
            f.attrs["meshesPath"]         = np.bytes_("meshes/")
            f.attrs["iterationEncoding"]  = np.bytes_("fileBased")
            f.attrs["iterationFormat"]    = np.bytes_("lasy_gaussian_laser_%T.h5")
            f.attrs["software"]           = np.bytes_("lasy")
            f.attrs["softwareVersion"]    = np.bytes_("0.6.0")
            f.attrs["comment"]            = np.bytes_(benchmark_lasy_file_comment)

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
            # E along y in Smilei, i.e. along x in LASY
            dset.attrs["polarization"]       = np.array([1. + 0j, 0. + 0j])

    # The namelist is executed by every MPI process: only one writes the file,
    # and all wait for it before reading it in LaserFromLasyAM
    try:
        from mpi4py import MPI
        mpi_rank                  = MPI.COMM_WORLD.Get_rank()
        mpi_barrier               = MPI.COMM_WORLD.Barrier
    except ImportError:
        mpi_rank                  = 0
        mpi_barrier               = lambda: None

    def check_benchmark_lasy_file_can_be_written(filename):
        # The file can be (over)written only if it does not exist or if it was written by this benchmark
        import os
        import h5py
        if os.path.exists(filename):
            try:
                with h5py.File(filename, "r") as f:
                    file_comment  = f.attrs.get("comment", b"")
            except OSError:
                file_comment      = b""
            if isinstance(file_comment, (bytes, np.bytes_)):
                file_comment      = file_comment.decode()
            if file_comment != benchmark_lasy_file_comment:
                raise Exception("The file %s exists and was not written by this benchmark, so it is not overwritten. "
                                "To read it, set write_benchmark_lasy_file = False; "
                                "to run the benchmark, choose another lasy_filename." % filename)

    if write_benchmark_lasy_file:
        # The check is done by all processes (so that all stop in case of error),
        # and all of them finish it before the file is written
        check_benchmark_lasy_file_can_be_written(lasy_filename)
        mpi_barrier()
        if mpi_rank == 0:
            write_lasy_file(lasy_filename)
        mpi_barrier()

    # The reference frequency is the one of the Main block, omega0_SI is checked against the file.
    # N_time is not given: the laser is injected during the whole time axis of the LASY file.
    LaserFromLasyAM(
        lasy_filename,
        omega0_SI = omega0_SI
    )

DiagFields(
    every = 100,
    fields = ["Br_m_mode_0", "Br_m_mode_1","Bl_m_mode_0","Bl_m_mode_1","Bt_m_mode_0","Bt_m_mode_1","Bt_mode_0","Bt_mode_1","Bl_mode_0","Bl_mode_1","Br_mode_0","Br_mode_1","Er_mode_0","Er_mode_1","Et_mode_0","Et_mode_1","El_mode_0","El_mode_1" ]
)