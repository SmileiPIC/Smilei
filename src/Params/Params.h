/*! @file Params.h

 @brief Params.h is the class that hold the simulation parameters and can read from a file the namelist

 @date 2013-02-15
 */

#ifndef Params_H
#define Params_H

#undef _POSIX_C_SOURCE
#undef _XOPEN_SOURCE

#include "Profile.h"
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
    extern double normal(double stddev);
}


// ---------------------------------------------------------------------------------------------------------------------
//! Params class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class Params {

public:
    //! Creator for Params
    Params(SmileiMPI*, std::vector<std::string>);

    //! destructor
    ~Params();

    //! compute grid-related parameters & apply normalization
    void compute();

    //! print a summary of the values in txt
    void print_init();
    //! Printing out some data at a given timestep
    void print_timestep(unsigned int itime, double time_dual, Timer & timer);
    void print_timestep_headers();

    //! Print information about the parallel aspects
    void print_parallelism_params(SmileiMPI* smpi);

    //! Tells whether standard output is this timestep
    bool printNow( int timestep ) {
        return (timestep % print_every == 0);
    }

    //! sets nDim_particle and nDim_field based on the geometry
    void setDimensions();

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
    double sim_time;

    /*! \brief Space resolution
     Number of cells in every direction in \f$ 2\pi/k_N \f$ where \f$ k_N=\omega_N/c \f$ is the normalization wavenumber
     */
    std::vector<double> res_space;

    //! local simulation box size in \f$2\pi/k_N \f$
    std::vector<double> sim_length;

    //! time during which the Maxwell's equations are not solved
    double time_fields_frozen;

    //! Boundary conditions for ElectroMagnetic Fields
    std::vector<std::string> bc_em_type_x;
    std::vector<std::string> bc_em_type_y;
    std::vector<std::string> bc_em_type_z;

    //Poisson solver
    //! Do we solve poisson
    bool solve_poisson;
    //! Maxium number of poisson iteration
    unsigned int poisson_iter_max;
    //! Maxium poisson error tolerated
    double poisson_error_max;

    //! Maxwell Solver (default='Yee')
    std::string maxwell_sol;

    //! Current spatial filter parameter: number of binomial pass
    unsigned int currentFilter_int;

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

    //! max value for dt (due to usual FDTD CFL condition: should be moved to ElectroMagn solver (MG))
    double dtCFL;

    //! number of cells in every direction of the local domain
    std::vector<unsigned int> n_space;

    //! number of cells in every direction of the global domain
    std::vector<unsigned int> n_space_global;

    //! spatial step (cell dimension in every direction)
    std::vector<double> cell_length;

    //! volume of cell (this will be removed by untructured mesh!)
    double cell_volume;

    //! wavelength (in SI units)
    double referenceAngularFrequency_SI;

    //! Oversize domain to exchange less particles
    std::vector<unsigned int> oversize;

    //! True if restart requested
    bool restart;

    //! frequency of exchange particles (default = 1, disabled for now, incompatible with sort)
    int exchange_particles_each;

    //! frequency to apply shrink_to_fit on particles structure
    int every_clean_particles_overhead;

    //! Total number of patches
    unsigned int tot_number_of_patches;
    //! Number of patches per direction
    std::vector<unsigned int> number_of_patches;
    //! Load balancing frequency
    int balancing_every;
    //! Load coefficient applied to a cell (default = 1)
    double coef_cell;
    //! Load coefficient applied to a frozen particle (default = 0.1)
    double coef_frozen;
    //! Return if number of patch = number of MPI process, to tune IO //ism
    bool one_patch_per_MPI;
    //! Compute an initially balanced patch distribution right from the start
    bool initial_balance;

    //! Tells whether there is a moving window
    bool hasWindow;

    //! Tells whether there is a species with Monte-Carlo Compton radiation
    bool hasMCRadiation;
    //! Tells whether there is a species with Continuous radiation loss.
    bool hasLLRadiation;
    //! Tells whether there is a species with the stochastic radiation loss
    //! of Niel et al.
    bool hasNielRadiation;

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
    void cleanup(SmileiMPI*);
    //! Method to find the numbers of requested species, sorted, and duplicates removed
    static std::vector<unsigned int> FindSpecies(std::vector<Species*>&, std::vector<std::string>);

    //! every for the standard pic timeloop output
    unsigned int print_every;

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

private:
    //! passing named command to python
    void runScript(std::string command, std::string name=std::string(""));

    //! Characters width for timestep output
    unsigned int timestep_width;
};

#endif
