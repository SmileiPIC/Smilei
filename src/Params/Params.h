/*! @file Params.h
 
 @brief Params.h is the class that hold the simulation parameters and can read from a file the namelist
 
 @date 2013-02-15
 */

#ifndef Params_H
#define Params_H

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

class SmileiMPI;
class Species;


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
    void print();
    
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
    
    //! normalization (used in the input files only)
    std::string sim_units;
    
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
    
    
    //! window simulation box size in number of cells
    int nspace_win_x;
    //! Time at which the moving window starts.
    double t_move_win;
    //! Velocity of the moving window along x in c.
    double vx_win;
    
    
    //! Clusters width
    //unsigned int clrw;
    int clrw;
    //! Number of cells per cluster
    int n_cell_per_cluster;
    
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
    double wavelength_SI;
    
    //! Oversize domain to exchange less particles
    std::vector<unsigned int> oversize;
    
    //! restart namelist
    bool restart;
    
    //! frequency of exchange particles (default = 1, disabled for now, incompatible with sort) 
    int exchange_particles_each;
    
    //! Number of MPI process per direction (default : as square as possible)
    std::vector<int> number_of_procs;
    
    //! global number of time exits (it will be used if not specified in various diags/fields)
    unsigned int global_every;
    
    //! string containing the whole clean namelist
    std::string namelist;
    
    //! call the python cleanup function and 
    //! check if python can be closed (e.g. there is no laser python profile)
    //! by calling the _keep_python_running python function (part of pycontrol.pyh)
    void cleanup();
    
    //! close Python
    static void closePython();
    
    
    //! Method to find the numbers of requested species, sorted, and duplicates removed
    static std::vector<unsigned int> FindSpecies(std::vector<Species*>&, std::vector<std::string>);
    

private:
    //! init python RTE
    void initPython(SmileiMPI*, std::vector<std::string>);
    
    //! passing named command to python
    void pyRunScript(std::string command, std::string name=std::string(""));
    
};

#endif
