/*! @file PicParams.h
 
 @brief PicParams.h is the class that hold the simulation parameters and can read from a file the namelist
 
 @date 2013-02-15
 */

#ifndef PICPARAMS_H
#define PICPARAMS_H

#include <vector>
#include <string>
#include "InputData.h"



// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each species
// ---------------------------------------------------------------------------------------------------------------------
struct SpeciesStructure {
    //! kind of species possible values: "ion" "eon" "test"
    std::string species_type;
    
    //! density profile
    std::string density_profile;
    
    //! initialization type. Possible values: "regular" "cold" "Maxwell-Juettner"
    std::string initialization_type;
    
    //! number of particles per cell
    unsigned int n_part_per_cell;
    
    //! coefficient on the maximum number of particles for the species
    double c_part_max;
    
    //! mass [electron mass]
    double mass;
    
    //! atomic number
    unsigned int atomic_number;
    
    //! charge [proton charge]
    short charge;
    
    //! density [\f$n_N=\epsilon_0\,m_e\,\omega_N^{2}/e^2\f$ ]
    double density;
    //! mean velocity in units of light velocity
    std::vector<double> mean_velocity; // must be params.nDim_field
    //! temperature [\f$m_e\,c^2\f$ ]
    std::vector<double> temperature;
    
    //! dynamics type. Possible values: "Norm" "Radiation Reaction"
    std::string dynamics_type;
    
    //! Time for which the species is frozen
    double time_frozen;
    
    //! logical true if particles radiate
    bool radiating;

    //! logical true if particles radiate
    bool isTest;
    
    //! Boundary conditions for particules
    std::string bc_part_type_long;
    std::string bc_part_type_trans;
    
    //! Ionization model per Specie (tunnel)
    std::string ionization_model;

    //! species geometry
    std::string species_geometry;

    //! vacuum lengths
    std::vector<double> vacuum_length;
    
    //! lengths related to the density profile definition
    std::vector<double> dens_length_x;
    
    //! lengths related to the density profile definition
    std::vector<double> dens_length_y;
    
    //! lengths related to the density profile definition
    std::vector<double> dens_length_z;
    
    //! doubles related to the density profile definition
    std::vector<double> dens_dbl_params;
    
    //! integer related to the density profile definition
    std::vector<short int> dens_int_params;
    
    //! slope lengths (symmetric for trapezoidal geometry, general for triangular geometry)
    std::vector<double> slope_length;
    
    //! left slope lengths(not symmetric for trapezoidal case)
    std::vector<double> left_slope_length;
    
    //! right slope lengths(not symmetric for trapezoidal case)
    std::vector<double> right_slope_length;
    
    //! cut parameter for a gaussian profile
    std::vector<double> cut;
    
    //! sigma parameter for a gaussian profile
    std::vector<double> sigma;
    
    //! plateau for a gaussian profile
    std::vector<double> plateau;
    
    //! polygonal density profile in x direction
    std::vector<double> x_density_coor;
    
    //! polygonal density profile relative values in x direction
    std::vector<double> density_rel_values_x;
    
    //! mode for 1D cos density profile
    double mode;
    
    //! fase  for 1D cos density profile
    double thetax;
    
    //! amplitude  for 1D cos density profile
    double ampl;
    
};



// ---------------------------------------------------------------------------------------------------------------------
//! PicParams class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class PicParams {
    
public:
    //! Creator for PicParams
    PicParams(InputData &);
    
    //! compute grid-related parameters & apply normalization
    void compute();
    
    //! compute species-related parameters & apply normalization
    void computeSpecies();
    
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
    
    //! conversion factor (=1 when normalized units, 2\pi when wavelength-related normalisations)
    double conv_fac;

    
    /*! \brief Time resolution.
     Number of timesteps in \f$ 2\pi/\omega_N \f$ where \f$ \omega_N \f$ is the normalization (plasma or laser) frequency
     */
    double res_time;
    
    //! simulation exit time in units of \f$ 2\pi/\omega_N \f$
    double sim_time;
    
    /*! \brief Space resolution.
     Number of cells in every direction in \f$ 2\pi/k_N \f$ where \f$ k_N=\omega_N/c \f$ is the normalization wavenumber
     */
    std::vector<double> res_space;
    
    //! local simulation box size in \f$2\pi/k_N \f$
    std::vector<double> sim_length;
    
    //! Boundary conditions for ElectroMagnetic Fields
    std::string bc_em_type_long;
    std::string bc_em_type_trans;
    
    
    //! window simulation box size in number of cells
    int nspace_win_x;
    //! Time at which the moving window starts.
    double t_move_win;
    //! Velocity of the moving window along x in c.
    double vx_win;
    
    //! Clusters width
    int clrw;
    
    //! initial number of species
    unsigned int n_species;
    
    //! parameters of the species
    std::vector<SpeciesStructure> species_param;
    
    //! initial number of particles
    unsigned int n_particles;
    
    //! number of total timesteps to perform in the simulation
    unsigned int n_time;
    
    //! dt for the simulation (CFL)
    double timestep;
    
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
	
	//! Timestep to dump everything
	unsigned int dump_step;
    
	//! Human minutes to dump everything
	double dump_minutes;
    
	//! exit once dump done
	bool exit_after_dump;
	
	//! check for file named "stop"
	bool check_stop_file;
	
	//! keep the last dump_file_sequence dump files
	unsigned int dump_file_sequence;
	
	//! restart namelist
	bool restart;
	
	//! frequency of exchange particles (default = 1, disabled for now, incompatible with sort) 
	int exchange_particles_each;
    
    //! Number of MPI process per direction (default : as square as possible)
    std::vector<int> number_of_procs;
    
    //! global number of time exits (it will be used if not specified in various diags/fields)
    unsigned int global_every;
};

#endif
