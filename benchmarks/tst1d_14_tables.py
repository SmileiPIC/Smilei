# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Descrition:
# Computation of the radiation and multiphoton Breit-Wheeler tables
#
# Purpose:
# To test that the generation and output process work well by generating small tables
#
# Validation:
# - external tables
#
# ----------------------------------------------------------------------------------------

import math
from math import sqrt

# _____________________________________________________________________________
# Main parameters

c = 299792458
lambdar = 1e-6
wr = 2*math.pi*c/lambdar

# ______________________________________________________________________________
# Namelists

Main(
    geometry = "1Dcartesian",

    interpolation_order = 2 ,

    cell_length = [0.1],
    grid_length  = [48*0.1],

    number_of_patches = [4],

    timestep = 0.09,
    simulation_time = 0.1,

    EM_boundary_conditions = [ ["periodic","periodic"] ],

    random_seed = smilei_mpi_rank,

    print_every = 1,

    reference_angular_frequency_SI = wr,

)


Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.0,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "vay",
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    boundary_conditions = [
    	["periodic", "periodic"],
    ],
)

Species(
    name = "positron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = 0.0,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "vay",
    radiation_model = "Niel",
    boundary_conditions = [
    	["periodic", "periodic"],
    ],
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    c_part_max = 20.0,
    mass = 0,
    charge = 0.,
    number_density = 0.0,
    temperature = [0.],
    mean_velocity = [0.0 ,0.0, 0.0],
    pusher = "norm",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = [["periodic", "periodic"]],
)

RadiationReaction(

    # General radiation parameters
    minimum_chi_discontinuous = 1e-2,
    compute_table = True,
    output_format = "hdf5",
    table_path = "./",

    # Parameter to generate the table integfochi used by the Monte-Carlo model
    integfochi_chipa_min = 1e-4,
    integfochi_chipa_max = 1e1,
    integfochi_dim = 16,

    # Parameters to generate the table h used by Niel et al.
    h_chipa_min = 1E-3,
    h_chipa_max = 1E1,
    h_dim = 16,

    # Parameter to generate the table xip used by the Monte-Carlo model
    #xip_chipa_min = 1e-4,
    #xip_chipa_max = 1e1,
    #xip_power = 4,
    xip_threshold = 1e-3,
    xip_chipa_dim = 4,
    xip_chiph_dim = 4,

)

MultiphotonBreitWheeler(
    compute_table = True,
    table_path = "./",
    #table_path = "/home/mathieu/Documents/Codes/particle_merging/databases"

    T_dim = 16,

    xip_chipa_dim = 4,
    xip_chiph_dim = 4,

)
