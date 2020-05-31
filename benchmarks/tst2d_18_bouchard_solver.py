# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------------------------------

#!/usr/bin/python
import math  as m
import numpy as np

# Space and time properties for the simulation (w/ respect of CFL)
# ---------------------------------------------------------------

l0 = 2.*m.pi                 # laser wavelength
t0 = l0                      # optical cycle

resx  = 8.                  # nb of cells in one laser wavelength x
resy  = 8.                  # nb of cells in one laser wavelength y

solver = 'Bouchard'
order=2

# One 2 pass current filter in y-direction, perpendicular to the x-velocity of the plasma

# My FIR Kernel
my_kernelFIR = [
    0.,
    0.01364035272420037,
    -3.5059674461847047e-17,
    -0.023941402143973,
    -1.8698493046318425e-16,
    0.046588211176174774,
    -1.2855213969343918e-16,
    -0.0951163438409052,
    -2.921639538487254e-16,
    0.3145154967908793,
    0.4999999999999999,
    0.3145154967908793,
    -2.921639538487254e-16,
    -0.0951163438409052,
    -1.2855213969343918e-16,
    0.046588211176174774,
    -1.8698493046318425e-16,
    -0.023941402143973,
    -3.5059674461847047e-17,
    0.01364035272420037,
    0.
]

CurrentFilter(
    model = "customFIR",
    passes = [2,2],
    kernelFIR = my_kernelFIR
)

fromcflfactor = 1.00 # have to be less than 1.
rest = 2*resx/fromcflfactor

T0 = int((15.+10.)*rest)
T2 = int((40.+10.)*rest)
T1 = int((30.+10.)*rest)

Lsim = [120*l0,60.*l0]               # dimension of the simulation
Tsim = 110*t0                      # duration of the simulation

# Conversion formula
# ------------------

T_keV   = 1./511.                # factor in order to express temperature in keV (T0 = 1/511 eV = 1 keV)

Th = 1e-12*T_keV
Tc = 1e-12*T_keV

# ---------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------

mfac = 1
dfac = 1

Main(
        geometry                       = "2Dcartesian",
        interpolation_order            = order,
        custom_oversize                = 10,
        timestep                       = t0/rest,
        simulation_time                = Tsim,
        cell_length                    = [l0/resx,l0/resx],
        grid_length                    = Lsim,
        number_of_patches              = [8,8],
        clrw                           = 1,
        maxwell_solver                 = solver,
        EM_boundary_conditions         = [
                                          ["periodic","periodic"],
                                          ["periodic","periodic"]
                                         ],
        random_seed                    = smilei_mpi_rank,
        solve_poisson                  = True,
        reference_angular_frequency_SI = 2.*m.pi*3e14/0.8,
        print_every                    = 10
)

LoadBalancing(
        initial_balance = True,
        every           = 1e6
)

# ---------------------------
# DEFINING SPECIES PROPERTIES
# ---------------------------

Species(
    name                           = "pon_c",
    particles_per_cell             = 1,
    # atomic_number                  = 1,
    mass                           = +1.0,
    charge                         = +1.0,
    number_density                 = constant(1,0,0), 
    position_initialization        = "centered",
    momentum_initialization        = "mj",
    time_frozen                    = 0,
    temperature                    = [1e-2],
    mean_velocity                  = [0.999999995, 0, 0],
    # thermal_boundary_temperature = [1e-3*T_keV],
    # thermal_boundary_velocity    = [0., 0., -beta],
    boundary_conditions            = [
                                      ['periodic','periodic'],
                                      ['periodic','periodic']
                                     ]
)

Species(
    name                           = "eon_c",
    particles_per_cell             = 1,
    # atomic_number                = 1,
    mass                           = +1.0,
    charge                         = -1.0,
    number_density                 = constant(1,0,0),
    position_initialization        = "pon_c",
    momentum_initialization        = "mj",
    time_frozen                    = 0,
    temperature                    = [1e-2],
    mean_velocity                  = [0.999999995, 0, 0],
    # thermal_boundary_temperature   = [10.*T_keV],
    # thermal_boundary_velocity      = [0., 0., -beta],
    boundary_conditions            = [
                                      ['periodic','periodic'],
                                      ['periodic','periodic']
                                     ]
)

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

# Diagnostics on Scalar and Fields in order to have an overall vizualisation
# --------------------------------------------------------------------------

DiagScalar(
  every = 2,
  vars = [
    "Utot",
    "Ukin_eon_c",
    "Ukin_pon_c",
    "Uelm",
    "Uelm_Bz_m",
    "Uelm_By_m",
    "Uelm_Bx_m",
    "Uelm_Ez",
    "Uelm_Ey",
    "Uelm_Ex",
    "Ubal",
    "Ubal_norm",
    "Ukin_bnd",
    "Uelm_bnd"]
)

# Fields Images ------------------------------------------------------------

dfac = 1.
Tsim_over_t0 = 110.
Lsim_over_l0 = [120.,60.]

DiagFields(
    every        = [0.*rest,Tsim_over_t0*rest,100],
    flush_every  = [0.*rest,Tsim_over_t0*rest,2000],
    subgrid      = np.s_[int(40*resy-10*resy):int(40*resy+20*resy), int(30*resy-10*resy):int(30*resy+10*resy)],
    fields       = ["Bz","Ey"]
)
