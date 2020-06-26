/*! @file Params.h

 @brief Params.h is the class that hold the simulation parameters and can read from a file the namelist

 @date 2013-02-15
 */

#ifndef Params_H
#define Params_H

#undef _POSIX_C_SOURCE
#undef _XOPEN_SOURCE

#include "Timer.h"
#include "codeConstants.h"

#include <vector>
#include <string>
#include <cstdlib>

#include <fstream>
#include <sstream>
#include <map>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <iterator>
#include <random>

class SmileiMPI;
class Species;
class Profile;

namespace Rand
{
extern std::random_device device;
extern std::mt19937 gen;

extern std::uniform_real_distribution<double> uniform_distribution;
extern double uniform();

extern std::uniform_real_distribution<double> uniform_distribution1;
extern double uniform1();

extern std::uniform_real_distribution<double> uniform_distribution2;
extern double uniform2();

extern std::normal_distribution<double> normal_distribution;
extern double normal( double stddev );
}


// ---------------------------------------------------------------------------------------------------------------------
//! Params class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class Params
{

public:
    //! Creator for Params
    Params( SmileiMPI *, std::vector<std::string> );

    //! destructor
    ~Params();

    //! compute grid-related parameters & apply normalization
    void compute();

    //! check if input parameters & apply normalizationare coherent
    void check_consistency();

    //! print a summary of the values in txt
    void print_init();
    //! Printing out some data at a given timestep
    void print_timestep( unsigned int itime, double time_dual, Timer &timer );
    void print_timestep_headers();

    //! Print information about the parallel aspects
    void print_parallelism_params( SmileiMPI *smpi );

    //! Tells whether standard output is this timestep
    bool printNow( int current_timestep )
    {
        return ( current_timestep % print_every == 0 );
    }

    //! sets nDim_particle and nDim_field based on the geometry
    void setDimensions();
    
    //! Return the species related to this field, or "" if not a species field
    std::string speciesField( std::string field_name );
    
    //! defines the geometry of the simulation
    std::string geometry;

    //! defines the interpolation/projection order
    unsigned int interpolation_order;

    //! number of space dimensions for the particles
    unsigned int nDim_particle;

    //! number of space dimensions for the fields
    unsigned int nDim_field;

    /*! \brief Time resolution
     Number of timesteps in \f$ 2\pi/\omega_N \f$ where \f$ \omega_N \f$ is the normalization (plasma or laser) frequency
     */
    double res_time;

    //! simulation exit time in units of \f$ 2\pi/\omega_N \f$
    double simulation_time;

    /*! \brief Space resolution
     Number of cells in every direction in \f$ 2\pi/k_N \f$ where \f$ k_N=\omega_N/c \f$ is the normalization wavenumber
     */
    std::vector<double> res_space;

    //! local simulation box size in \f$2\pi/k_N \f$
    std::vector<double> grid_length;

    //! time during which the Maxwell's equations are not solved
    double time_fields_frozen;

    //! Boundary conditions for ElectroMagnetic Fields
    std::vector< std::vector<std::string> > EM_BCs;
    //! k parameters for some kinds of ElectroMagnetic boundary conditions
    std::vector< std::vector<double> > EM_BCs_k;
    //! Are open boundaries used ?
    bool open_boundaries;
    bool save_magnectic_fields_for_SM;

    //! Boundary conditions for Envelope Field
    std::vector< std::vector<std::string> > Env_BCs;

    //! Define if the ponderomotive force is computed (default = false)
    //bool ponderomotive_force;

    //! Define if laser envelope model is used (default = false)
    bool Laser_Envelope_model=false;

    //! Define if there is at least one species ionized by envelope
    bool envelope_ionization_is_active = false;
    double envelope_ellipticity = 0.; // (0: linear polarization, 1: circular polarization)
    double envelope_polarization_phi = 0.; // used only for envelope ionization; in radians, angle with the xy plane
    // define the solver for the envelope equation
    std::string envelope_solver;

    
    //Poisson solver
    //! Do we solve poisson
    bool solve_poisson;
    //! Maxium number of poisson iteration
    unsigned int poisson_max_iteration;
    //! Maxium poisson error tolerated
    double poisson_max_error;

    //"Relativistic" Poisson solver
    //! Do we solve "relativistic poisson problem" for relativistic species
    bool solve_relativistic_poisson;
    //! Maxium number of relativistic poisson iteration
    unsigned int relativistic_poisson_max_iteration;
    //! Maxium relativistic poisson error tolerated
    double relativistic_poisson_max_error;

    //! Do we need to exchange full B (default=0 <=> only 2 components are exchanged by dimension)
    bool full_B_exchange;
    //! Do we need to exchange full A,Phi,Chi (default=0 <=> only 2 components are exchanged by dimension)
    bool full_Envelope_exchange;
    
    //! Maxwell Solver (default='Yee')
    std::string maxwell_sol;

    //! Current spatial filter: number of binomial passes
    std::vector<unsigned int> currentFilter_passes;
    std::string currentFilter_model;
    std::vector<double> currentFilter_kernelFIR;

    //! is Friedman filter applied [Greenwood et al., J. Comp. Phys. 201, 665 (2004)]
    bool Friedman_filter;

    //! Fridman filtering parameter [real between 0 and 1]
    double Friedman_theta;

    //! Clusters width
    //unsigned int clrw;
    int clrw;
    //! Number of cells per cluster
    int n_cell_per_patch;

    //! initial number of particles
    unsigned int n_particles;

    //! number of total timesteps to perform in the simulation
    unsigned int n_time;

    //! dt for the simulation
    double timestep;

    //! Number of modes
    unsigned int nmodes;

    //! Number of modes for relativistic field initialization
    unsigned int nmodes_rel_field_init;

    //! max value for dt (due to usual FDTD CFL condition: should be moved to ElectroMagn solver (MG))
    double dtCFL;

    //! number of cells in every direction of the local domain
    std::vector<unsigned int> n_space;
    
    //! number of cells in every direction of the local domain (can be different from 1 MPI process to another)
    std::vector<unsigned int> n_space_region;
    
    std::vector<unsigned int> number_of_region;
    std::vector< std::vector<int> > offset_map;
    std::vector< std::vector< std::vector<int> > > map_rank;
    std::vector<int> coordinates;
    
    //! number of cells in every direction of the global domain
    std::vector<unsigned int> n_space_global;

    //! spatial step (cell dimension in every direction)
    std::vector<double> cell_length;

    //! Size of a patch in each direction
    std::vector<double> patch_dimensions;

    //! volume of cell (this will be removed by untructured mesh!)
    double cell_volume;

    //! wavelength (in SI units)
    double reference_angular_frequency_SI;

    //! Oversize domain to exchange less particles
    std::vector<unsigned int> oversize;
    unsigned int custom_oversize ;
    std::vector<unsigned int> region_oversize;
    unsigned int custom_region_oversize ;
    //! Number of damping cells
    std::vector<unsigned int> number_of_damping_cells;

    unsigned int pseudo_spectral_guardells;
    bool apply_rotational_cleaning;
    
    //! True if restart requested
    bool restart;

    //! frequency of exchange particles (default = 1, disabled for now, incompatible with sort)
    int exchange_particles_each;
    
    //! frequency to apply shrinkToFit on particles structure
    int every_clean_particles_overhead;

    //! Total number of patches
    unsigned int tot_number_of_patches;
    //! Number of patches per direction
    std::vector<unsigned int> number_of_patches;
    //! Domain decomposition
    std::string patch_arrangement;

    //! Time selection for adaptive vectorization
    TimeSelection *adaptive_vecto_time_selection;
    //! Flag for the adaptive vectorization
    bool has_adaptive_vectorization;

    //! Time selection for load balancing
    TimeSelection *load_balancing_time_selection;
    //! True if must balance at some point
    bool has_load_balancing;
    //! Load coefficient applied to a cell (default = 1)
    double cell_load;
    //! Load coefficient applied to a frozen particle (default = 0.1)
    double frozen_particle_load;
    //! Return if number of patch = number of MPI process, to tune IO //ism
    bool one_patch_per_MPI;
    //! Compute an initially balanced patch distribution right from the start
    bool initial_balance;

    //! String containing the vectorization mode: off, on, adaptive, adaptive_mixed_sort
    std::string vectorization_mode;
    //! Initial state of the patches in adaptive mode
    std::string adaptive_default_mode;

    //! Tells whether there is a moving window
    bool hasWindow;

    //! Tells whether there is a species with Monte-Carlo Compton radiation
    bool hasMCRadiation;
    //! Tells whether there is a species with Continuous radiation loss.
    bool hasLLRadiation;
    //! Tells whether there is a species with the stochastic radiation loss
    //! of Niel et al.
    bool hasNielRadiation;
    //! Tells whether there is w/out radiation reaction but for which a RadiationSpectrum diag is called
    bool hasDiagRadiationSpectrum;

    //! Tells whether there is a species with multiphoton Breit-Wheeler
    bool hasMultiphotonBreitWheeler;

    //! Log2 of the number of patch in the whole simulation box in every direction.
    //! The number of patch in a given direction MUST be a power of 2 and is 2^(mi[i]).
    std::vector<unsigned int> mi;

    //! string containing the whole clean namelist
    std::string namelist;

    //! call the python cleanup function and
    //! check if python can be closed (e.g. there is no laser python profile)
    //! by calling the _keep_python_running python function (part of pycontrol.pyh)
    void cleanup( SmileiMPI * );

    //! Method to find the numbers of requested species, sorted, and duplicates removed
    static std::vector<unsigned int> FindSpecies( std::vector<Species *> &, std::vector<std::string> );

    //! every for the standard pic timeloop output
    unsigned int print_every;
    
    // Double grids parameters (particles and fields)
    std::vector<unsigned int> global_factor;
    void uncoupled_decomposition();
    void uncoupled_decomposition_1D();
    void uncoupled_decomposition_2D();
    void uncoupled_decomposition_3D();
    void print_uncoupled_params();
    bool uncoupled_grids;

    // PXR parameters
    bool  is_spectral=false ;
    bool  is_pxr=false ;
    int   norderx = 2;
    int   nordery = 2;
    int   norderz = 2;
    std::vector<int> norder;

    //! Boolean for printing the expected disk usage or not
    bool print_expected_disk_usage;

    //! Random seed
    unsigned int random_seed;

    // ---------------------------------------------
    // Constants
    // ---------------------------------------------

    //! Fine structure constant
    const double fine_struct_cst = 7.2973525698e-3;

    //! Reduced Planck Constant (J.s)
    const double red_planck_cst = 1.054571628E-34;

    //! Electron mass
    const double electron_mass = 9.109382616e-31;

    //! Speed of light in vacuum (m/s)
    const double c_vacuum = 299792458;

    //! passing named command to python
    void runScript( std::string command, std::string name, PyObject * );

    //! Characters width for timestep output
    unsigned int timestep_width;

    bool cell_sorting;
};

#endif
