/*! @file PicParams.h

  @brief PicParams.h is the class that hold the simulation parameters and can read from a file the namelist

  @author tommaso vinci
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
struct LaserStructure {

    //! Laser field amplitude
    double a0;

    //! Laser angle
    double angle;

    //! Laser delta (ellipticity parameter)
    double delta;

    //! Laser profile
    std::string time_profile;

    //! int vector for laser parameters
    std::vector<int> int_params;

    //! double vector for laser parameters
    std::vector<double> double_params;
};



// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each species
// ---------------------------------------------------------------------------------------------------------------------
struct SpeciesStructure {
    //! kind of species possible values: "ion" "electron" "test"
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

    //! Boundary conditions for particules
    std::string bc_part_type;

    //! Ionization model per Specie (tunnel)
    std::string ionization_model;

};



// ---------------------------------------------------------------------------------------------------------------------
//! PicParams class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class PicParams {

public:
    //! Creator for PicParams
    PicParams(InputData &);
    void compute();
    void print();
    void setDimensions();

    //! defines the geometry of the simulation
    std::string geometry;

    //! defines the interpolation/projection order
    unsigned int interpolation_order;

    //! number of space dimensions for the particles
    unsigned int nDim_particle;

    //! number of space dimensions for the fields
    unsigned int nDim_field;

    /*! \brief Time resolution.
      Number of timesteps in \f$ 2\pi/\omega_N \f$ where \f$ \omega_N \f$ is the normalization (plasma or laser) frequency
    */
    unsigned int res_time;

    //! simulation exit time in units of \f$ 2\pi/\omega_N \f$
    double sim_time;

    /*! \brief Space resolution.
      Number of cells in every direction in \f$ 2\pi/k_N \f$ where \f$ k_N=\omega_N/c \f$ is the normalization wavenumber
    */
    std::vector<unsigned int> res_space;

    //! local simulation box size in \f$2\pi/k_N \f$
    std::vector<double> sim_length;

    //! plasma geometry
    std::string plasma_geometry;

    //!
    std::vector<double> density_double_params;

    //! plasma lengths
    std::vector<double> plasma_length;

    //! vacuum lengths
    std::vector<double> vacuum_length;

    //! slope lengths (symmetric for trapezoidal geometry, general for triangular geometry)
    std::vector<double> slope_length;
    
    //! left slope lengths(not symmetric for trapezoidal case)
    std::vector<double> left_slope_length;
    
    //! right slope lengths(not symmetric for trapezoidal case)
    std::vector<double> right_slope_length;

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

    //! physicist (cgs with temperatures in eV) theorist
    std::string sim_units;


    //! initial number of laser pulses
    unsigned int n_laser;

    //! laser parameters
    std::vector<LaserStructure> laser_param;

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
	
	//! enable sort particles (default = yes) 
	bool use_sort_particles;

	//! frequency of exchange particles (default = 1) 
	int exchange_particles_each;

};

#endif
