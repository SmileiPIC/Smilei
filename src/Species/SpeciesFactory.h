// -----------------------------------------------------------------------------
//
//! \file SpeciesFactory.h
//
//! \brief Contains the class SpeciesFactory that manages the
//!        initialization of the species. For the moment,
//!        only one kind of species exist.
//
// -----------------------------------------------------------------------------

#ifndef SPECIESFACTORY_H
#define SPECIESFACTORY_H

#include "Species.h"

#ifdef _VECTO
#include "SpeciesV.h"
#include "SpeciesVAdaptiveMixedSort.h"
#include "SpeciesVAdaptive.h"
#endif

#include "ParticlesFactory.h"
#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "PartBoundCond.h"

#include "Params.h"
#include "Patch.h"

#include "ParticleData.h"

#include "ParticleCreator.h"

#include "Tools.h"
#include "H5.h"
#ifdef SMILEI_USE_NUMPY
#include <numpy/arrayobject.h>
#endif

class SpeciesFactory
{
public:
    static Species *create( Params &params, int ispec, Patch *patch )
    {

        unsigned int tot_species_number = PyTools::nComponents( "Species" );

        std::string species_name;
        if( ! PyTools::extractOrNone( "name", species_name, "Species", ispec ) ) {
            std::ostringstream name( "" );
            name << "species" << std::setfill( '0' ) << std::setw( log10( tot_species_number )+1 ) << ispec;
            species_name = name.str();
        }

        if( patch->isMaster() ) {
            MESSAGE( 1,"");
            MESSAGE( 1, "Creating Species #" << ispec << ": " << species_name );
        }

        if( species_name.size() < 2 ) {
            ERROR_NAMELIST("For species #" << ispec << ", name cannot be only 1 character",
            LINK_NAMELIST + std::string("#name"));
        }

        if( species_name.substr(0,2) == "m_" ) {
            ERROR_NAMELIST("For species #" << ispec << ", name cannot start  with `m_`",
            LINK_NAMELIST + std::string("#name"));
        }

        // Extract type of species dynamics from namelist
        std::string pusher = "boris"; // default value
        PyTools::extract( "pusher", pusher, "Species", ispec );

        // Extract type of species radiation from namelist
        std::string radiation_model = "none"; // default value
        PyTools::extract( "radiation_model", radiation_model, "Species", ispec );
        // Cancelation of the letter case for `radiation_model`
        std::transform( radiation_model.begin(), radiation_model.end(), radiation_model.begin(), ::tolower );
        // Name simplification
        if( radiation_model=="monte-carlo" ) {
            radiation_model="mc";
        }
        if( radiation_model=="landau-lifshitz" ) {
            radiation_model="ll";
        }
        if( radiation_model=="corrected-landau-lifshitz" ) {
            radiation_model="cll";
        }


        // Extract mass from namelist
        double mass;
        PyTools::extract( "mass", mass, "Species", ispec );

        // Create species object
        Species *this_species = NULL;

        // Species type
        if ( params.vectorization_mode == "off" ) {
            this_species = new Species( params, patch );
        } 
        #ifdef _VECTO
        else if( params.vectorization_mode == "on" ) {
            this_species = new SpeciesV( params, patch );
        } else if( params.vectorization_mode == "adaptive_mixed_sort" ) {
            this_species = new SpeciesVAdaptiveMixedSort( params, patch );
        } else if( params.vectorization_mode == "adaptive" ) {
            this_species = new SpeciesVAdaptive( params, patch );
        }
        #endif

        // Particles
        if( mass > 0. ) {
            // Dynamics of the species
            if( pusher == "boris"
                    || pusher == "borisnr"
                    || pusher == "vay"
                    || pusher == "higueracary"
                    || pusher == "ponderomotive_boris" ) {
                // Species with relativistic Boris pusher if  =='boris'
                // Species with nonrelativistic Boris pusher == 'borisnr'
                // Species with J.L. Vay pusher if == "vay"
                // Species with Higuary Cary pusher if == "higueracary"
            } else {
                ERROR_NAMELIST( "For species `" << species_name << "`, pusher must be 'boris', 'borisnr', 'vay', 'higueracary', 'ponderomotive_boris'",
                LINK_NAMELIST + std::string("#pusher") );
            }
            this_species->pusher_name_ = pusher;
            MESSAGE( 2, "> Pusher: " << this_species->pusher_name_ );

            // Radiation model of the species
            // Species with a Monte-Carlo process for the radiation loss
            if( radiation_model=="mc" ) {
                this_species->particles->isQuantumParameter = true;
                this_species->particles->isMonteCarlo = true;
                this_species->radiating_ = true;
            }
            // Species with another radiation loss model
            else if( ( radiation_model=="ll" )
                     || ( radiation_model=="cll" )
                     || ( radiation_model=="niel" ) ) {
                this_species->particles->isQuantumParameter = true;
                this_species->radiating_ = true;
            } else if( radiation_model=="diagradiationspectrum" ) {
                    this_species->particles->isQuantumParameter = true;
                    this_species->radiating_ = true;
            } else if( radiation_model != "none" ) {
                ERROR_NAMELIST( "For species `" << species_name
                       << " radiation_model must be 'none',"
                       << " 'Landau-Lifshitz' ('ll'),"
                       << " 'corrected-Landau-Lifshitz' ('cll'),"
                       << " 'Niel' ('niel') or 'Monte-Carlo' ('mc')",
                       LINK_NAMELIST + std::string("#radiation_model") );
            }

            this_species->radiation_model_ = radiation_model;

            if( radiation_model == "ll" ) {
                MESSAGE( 2, "> Radiating species with the classical Landau-Lifshitz radiating model" );
            } else if( radiation_model == "cll" ) {
                MESSAGE( 2, "> Radiating species with the quantum corrected Landau-Lifshitz radiating model" );
            } else if( radiation_model == "niel" ) {
                MESSAGE( 2, "> Radiating species with the stochastic model of Niel et al." );
            } else if( radiation_model == "mc" ) {
                MESSAGE( 2, "> Radiating species with the stochastic Monte-Carlo model" );
            } else if( radiation_model == "diagradiationspectrum" ) {
                MESSAGE(2,"> Radiating species without backreaction: a `DiagRadiationSpectrum` can be applied to this species");
            } else if( radiation_model != "none" ) {
                MESSAGE( 2, "> Radiating species with model: `" << radiation_model << "`" );
            }

            // Non compatibility
            if( ( pusher=="borisnr" )
                    && ( radiation_model=="mc"
                         || radiation_model=="ll"
                         || radiation_model=="cll"
                         || radiation_model=="niel"
                         || radiation_model=="diagradiationspectrum") ) {
                ERROR_NAMELIST( "For species `" << species_name
                       << "` radiation_model `"
                       << radiation_model
                       << "` is not compatible with pusher "
                       << pusher,
                   LINK_NAMELIST + std::string("#radiation_model") );

            }

            // Non compatibility
            if( ( mass>1 )
                    && ( radiation_model=="mc"
                         || radiation_model=="ll"
                         || radiation_model=="cll"
                         || radiation_model=="niel"
                         || radiation_model=="diagradiationspectrum") ) {
                ERROR_NAMELIST( "For species `" << species_name
                       << "` radiation_model `"
                       << radiation_model
                       << "` is only compatible with electron and positron species (charge <= 1).",
                    LINK_NAMELIST + std::string("#radiation_model"));

            }

        }

        // Photon species
        else if( mass == 0 ) {
            // Photon can not radiate
            this_species->radiation_model_ = "none";
            this_species-> pusher_name_ = "norm";

            MESSAGE( 2, "> " <<species_name <<" is a photon species (mass==0)." );
            //MESSAGE( 2, "> Radiation model set to none." );
            MESSAGE( 2, "> Pusher set to norm." );
        }

        this_species->name_ = species_name;
        this_species->mass_ = mass;
        this_species->species_number_ = ispec;

        // Vectorized operators
        if( params.vectorization_mode == "off" ) {
            this_species->vectorized_operators = false;
        } else if( params.vectorization_mode == "on" || params.vectorization_mode == "adaptive_mixed_sort" || params.vectorization_mode == "adaptive" ) {
            this_species->vectorized_operators = true;
        }

        // Extract various parameters from the namelist

        // Monte-Carlo Photon emission properties
        if( mass > 0. ) {
            if( this_species->radiation_model_ == "mc" ) {
                if( PyTools::extractOrNone( "radiation_photon_species", this_species->radiation_photon_species, "Species", ispec ) ) {

                    MESSAGE( 3, "| Macro-photon emission activated" );

                    // Species that will receive the emitted photons
                    if( !this_species->radiation_photon_species.empty() ) {
                        MESSAGE( 3, "| Emitted photon species set to `" << this_species->radiation_photon_species << "`" );
                    } else {
                        ERROR_NAMELIST( " The radiation photon species (`radiation_photon_species`) is not specified.",
                        LINK_NAMELIST + std::string("#radiation_photon_species"))
                    }

                    // Number of photons emitted per Monte-Carlo event
                    PyTools::extract( "radiation_photon_sampling",
                        this_species->radiation_photon_sampling_, "Species", ispec );
                    if( this_species->radiation_photon_sampling_ < 1 ) {
                        ERROR_NAMELIST( "For species '" << species_name
                               << "' radiation_photon_sampling should be > 1",
                           LINK_NAMELIST + std::string("#radiation_photon_sampling") );
                    }
                    MESSAGE( 3, "| Number of macro-photons emitted per MC event: "
                             << this_species->radiation_photon_sampling_ );

                    // Photon energy threshold
                    PyTools::extract( "radiation_photon_gamma_threshold",
                        this_species->radiation_photon_gamma_threshold_, "Species", ispec );

                    MESSAGE( 3, "| Photon energy threshold for macro-photon emission: "
                             << this_species->radiation_photon_gamma_threshold_ );
                    
                    // Creation of the photon particles object to receive the emitted photons
                    // if( !this_species->radiation_photon_species.empty() ) {
                    //     this_species->radiated_photons_ = ParticlesFactory::create( params );
                    //     //this_species->radiated_photons_->initialize( 0, params.nDim_particle, params.keep_position_old);
                    // }
                             
                    
                // else, no emitted macro-photons
                } else {
                    MESSAGE( 3, "| Macro-photon emission not activated" );
                    this_species->radiated_photons_ = NULL;
                    this_species->photon_species_   = NULL;
                }

            }
        }

        // Multiphoton Breit-Wheeler
        if( mass == 0 ) {
            
            //Photons can not radiate
            this_species->radiated_photons_ = NULL;
            this_species->photon_species_   = NULL;
            
            // If this_species->multiphoton_Breit_Wheeler
            if( PyTools::extractV( "multiphoton_Breit_Wheeler", this_species->multiphoton_Breit_Wheeler_, "Species", ispec ) ) {
                // If one of the species is empty
                if( this_species->multiphoton_Breit_Wheeler_[1].empty() || this_species->multiphoton_Breit_Wheeler_[0].empty() ) {
                    ERROR_NAMELIST(  "For species '" << species_name
                            << "' multiphoton_Breit_Wheeler can not be empty,"
                            << " select electron and positron species.",
                            LINK_NAMELIST + std::string("#multiphoton_Breit_Wheeler"));
                } else {
                    // Activation of the additional variables
                    this_species->particles->isQuantumParameter = true;
                    this_species->particles->isMonteCarlo = true;

                    MESSAGE( 2, "> Decay into pair via the multiphoton Breit-Wheeler activated" );
                    MESSAGE( 3, "| Generated electrons and positrons go to species: "
                             << this_species->multiphoton_Breit_Wheeler_[0]
                             << " & " << this_species->multiphoton_Breit_Wheeler_[1] );

                    // Number of emitted particles per MC event
                    this_species->mBW_pair_creation_sampling_.resize( 2 );
                    if( !PyTools::extractV( "multiphoton_Breit_Wheeler_sampling",
                                           this_species->mBW_pair_creation_sampling_, "Species", ispec ) ) {
                        this_species->mBW_pair_creation_sampling_[0] = 1;
                        this_species->mBW_pair_creation_sampling_[1] = 1;
                    }
                    MESSAGE( 3, "| Number of emitted macro-particles per MC event: "
                             << this_species->mBW_pair_creation_sampling_[0]
                             << " & " << this_species->mBW_pair_creation_sampling_[1] );
                }
            }
        }

        // Particle Merging

        // Extract merging method
        this_species->merging_method_ = "none"; // default value
        this_species->has_merging_ = false; // default value
        PyTools::extract( "merging_method", this_species->merging_method_, "Species", ispec );

        // Cancelation of the letter case for `merging_method_`
        std::transform( this_species->merging_method_.begin(),
                        this_species->merging_method_.end(),
                        this_species->merging_method_.begin(), ::tolower );

        if( ( this_species->merging_method_ != "vranic_spherical" ) &&
            ( this_species->merging_method_ != "vranic_cartesian" ) &&
            ( this_species->merging_method_ != "none" ) ) {
            ERROR_NAMELIST( "In Species " << this_species->name_
                            << ": merging method not valid, must be `vranic_spherical`, `vranic_cartesian` or `none`",
                            LINK_NAMELIST + std::string("#particle-merging") );
        }

        // if( params.vectorization_mode == "off" && this_species->merging_method_ != "none" ) {
        //     ERROR_NAMELIST( "In Species " << this_species->name_
        //            << ": particle merging only available with `vectorization_mode` = `on` or `adaptive`" );
        // }

        if ( this_species->merging_method_ != "none" ) {

            // get parameter "every" which describes a timestep selection
            if( !this_species->merging_time_selection_ ) {
                this_species->merging_time_selection_ = new TimeSelection(
                    PyTools::extract_py( "merge_every", "Species", ispec ), "Particle merging"
                );
            }

            // get extra parameters
            // Minimum particle number per packet to merge
            PyTools::extract( "merge_min_packet_size", this_species->merge_min_packet_size_ , "Species", ispec );
            if (this_species->merge_min_packet_size_ < 4 && this_species->mass_ > 0)
            {
                ERROR_NAMELIST( "In Species " << this_species->name_
                       << ": minimum number of particles per merging packet "
                       << "(`merge_min_packet_size`)"
                       << "must be above or equal to 4.",
                       LINK_NAMELIST + std::string("#particle-merging"));
            }
            if (this_species->merge_min_packet_size_ < 4 && this_species->mass_ == 0)
            {
                ERROR_NAMELIST( "In Species " << this_species->name_
                       << " of type photon"
                       << ": minimum number of particles per merging packet "
                       << "(`merge_min_packet_size`)"
                       << "must be above or equal to 4.",
                       LINK_NAMELIST + std::string("#particle-merging"));
            }

            // Maximum particle number per packet to merge
            PyTools::extract( "merge_max_packet_size", this_species->merge_max_packet_size_ , "Species", ispec );
            if (this_species->merge_max_packet_size_ < 4 && this_species->mass_ > 0)
            {
                ERROR_NAMELIST( "In Species " << this_species->name_
                       << ": maximum number of particles per merging packet "
                       << "(`merge_max_packet_size`)"
                       << "must be above or equal to 4.",
                       LINK_NAMELIST + std::string("#particle-merging"));
            }
            if (this_species->merge_max_packet_size_ < 4 && this_species->mass_ == 0)
            {
                ERROR_NAMELIST( "In Species " << this_species->name_
                       << " of type photon"
                       << ": maximum number of particles per merging packet "
                       << "(`merge_max_packet_size`)"
                       << "must be above or equal to 4.",
                       LINK_NAMELIST + std::string("#particle-merging"));
            }
            if (this_species->merge_max_packet_size_ < this_species->merge_min_packet_size_) {
                ERROR_NAMELIST( "In Species " << this_species->name_
                       << ": maximum number of particles per merging packet "
                       << "(`merge_max_packet_size`)"
                       << "must be below or equal to the minimum particle number"
                       << " per merging packet (`merge_min_packet_size`)",
                       LINK_NAMELIST + std::string("#particle-merging"));
            }

            // Minimum momentum cell length for the momentum discretization
            if( PyTools::extractV( "merge_min_momentum_cell_length",
                                  this_species->merge_min_momentum_cell_length_ ,
                                  "Species", ispec ) ) {
                for (unsigned int i = 0 ; i < 3 ; i++) {
                    if (this_species->merge_min_momentum_cell_length_[i] <= 0) {
                        ERROR_NAMELIST( "In Species " << this_species->name_
                                 << ": The minimal momentum cell length "
                                 << "(`merge_min_particles_per_cell`)"
                                 << " must be above 0 ("
                                 << this_species->merge_min_momentum_cell_length_[i]
                                 << ")",
                             LINK_NAMELIST + std::string("#particle-merging"));
                    }
                }
            } else {
                ERROR_NAMELIST(
                    "In Species " << this_species->name_ 
                    << ": merge_min_momentum_cell_length should be a list of floats",
                    LINK_NAMELIST + std::string("#particle-merging")
                );
            }

            // Read and check the threshold on the number of particles per cell
            PyTools::extract( "merge_min_particles_per_cell", this_species->merge_min_particles_per_cell_ , "Species", ispec );
            if (this_species->merge_min_particles_per_cell_ < 4) {
                ERROR_NAMELIST(
                    "In Species " << this_species->name_
                    << ": The threshold on the number of particles per cell "
                    << "(`merge_min_particles_per_cell`)"
                    << "must be above or equal to 4",
                    LINK_NAMELIST + std::string("#particle-merging")
                );
            }

            // Read flag to activate the accumulation correction
            PyTools::extract( "merge_accumulation_correction", this_species->merge_accumulation_correction_ , "Species", ispec );

            // Momentum cell discretization
            if( PyTools::extractV( "merge_momentum_cell_size",
                                  this_species->merge_momentum_cell_size_ ,
                                  "Species", ispec ) ) {
                for (unsigned int i = 0 ; i < 3 ; i++) {
                    if (this_species->merge_momentum_cell_size_[i] <= 0) {
                        ERROR_NAMELIST(
                            "In Species " << this_species->name_
                            << ": The momentum cell discretization can not be equal or below 0 "
                            << "(`merge_momentum_cell_size_`).",
                            LINK_NAMELIST + std::string("#particle-merging")
                        );
                    }
                }
            } else {
                ERROR_NAMELIST( 
                    "In Species " << this_species->name_ 
                    << ": merge_momentum_cell_size should be a list of integers",
                    LINK_NAMELIST + std::string("#particle-merging")
                );
            }


            // Momentum cell discretization
            std::string discretization_scale;
            PyTools::extract( "merge_discretization_scale",
                              discretization_scale, "Species", ispec );
            if (discretization_scale == "linear") {
                this_species->merge_log_scale_ = false;
            } else if (discretization_scale == "log") {
                this_species->merge_log_scale_ = true;
                if (this_species->merge_accumulation_correction_ == true)
                {
                    this_species->merge_accumulation_correction_ = false;
                }
            } else {
                ERROR_NAMELIST(
                    "In Species " << this_species->name_
                    << ": discretization_scale should be `linear` or `log`",
                    LINK_NAMELIST + std::string("#particle-merging")
                );
            }

            // Minimum momentum in log scale
            PyTools::extract( "merge_min_momentum",
                              this_species->merge_min_momentum_log_scale_,
                              "Species", ispec );
            if (this_species->merge_min_momentum_log_scale_ <= 0) {
                ERROR_NAMELIST(
                    "In Species " << this_species->name_
                    << ": merge_min_momentum should be above 0.",
                    LINK_NAMELIST + std::string("#particle-merging")
                );
            }

            // We activate the merging
            this_species->has_merging_ = true;
        }

        // Information about the merging process
        if( this_species->merging_method_ != "none" ) {
            MESSAGE( 2, "> Particle merging with the method: "
                     << this_species->merging_method_ );
            MESSAGE( 3, "| Merging time selection: "
                     << this_species->merging_time_selection_->info() );
            if (this_species->merge_log_scale_) {
                MESSAGE( 3, "| Discretization scale: log");
                MESSAGE( 3, "| Minimum momentum: " << std::scientific << std::setprecision(5)
                << this_species->merge_min_momentum_log_scale_);
            } else {
                MESSAGE( 3, "| Discretization scale: linear");
                if (this_species->merge_accumulation_correction_) {
                    MESSAGE( 3, "| Accumulation correction activated");
                } else {
                    MESSAGE( 3, "| Accumulation correction disabled");
                }
            }
            MESSAGE( 3, "| Momentum cell discretization: "
                     << this_species->merge_momentum_cell_size_[0] << " "
                     << this_species->merge_momentum_cell_size_[1] << " "
                     << this_species->merge_momentum_cell_size_[2] << " ");
            MESSAGE( 3, "| Minimum momentum cell length: "
                    << std::scientific
                    << this_species->merge_min_momentum_cell_length_[0] << " "
                    << this_species->merge_min_momentum_cell_length_[1] << " "
                    << this_species->merge_min_momentum_cell_length_[2] << " ");
            MESSAGE( 3, "| Minimum particle number per cell: "
                     << std::fixed
                     << this_species->merge_min_particles_per_cell_ );
            MESSAGE( 3, "| Minimum particle packet size: "
                     << this_species->merge_min_packet_size_ );
            MESSAGE( 3, "| Maximum particle packet size: "
                     << this_species->merge_max_packet_size_ );
        }

        // Position initialization
        PyObject *py_pos_init = PyTools::extract_py( "position_initialization", "Species", ispec );
        if( PyTools::py2scalar( py_pos_init, this_species->position_initialization_ ) ) {
            if( this_species->position_initialization_.empty() ) {
                ERROR_NAMELIST(
                    "For species '" << species_name << "' empty position_initialization.",
                    LINK_NAMELIST + std::string("#species")
                );
            // Regular, random, centered
            } else if(    this_species->position_initialization_=="centered" and params.geometry == "AMcylindrical"){
                ERROR_NAMELIST(
                    "For species '" << species_name 
                    << "' centered position_initialization is not supported in AM geometry.",
                    LINK_NAMELIST + std::string("#species")
                );
            } else if(    this_species->position_initialization_=="regular"
                       || this_species->position_initialization_=="random"
                       || this_species->position_initialization_=="centered" ) {
            //    ;
            // Copy positions of other species
            } else if( PyTools::isSpecies( this_species->position_initialization_ ) ) {
                // Find the linked species
                bool ok = false;
                for( unsigned int ispec2 = 0; ispec2<patch->vecSpecies.size(); ispec2++ ) {
                    if( patch->vecSpecies[ispec2]->name_ == this_species->position_initialization_ ) {
                        ok = true;
                        this_species->position_initialization_on_species_index = ispec2;
                        break;
                    }
                }
                // std::cerr << this_species->position_initialization_on_species_index
                //           << std::endl;
                // The link species must already exist
                if( ok == false ) {
                    ERROR_NAMELIST( 
                        "For species '" << species_name 
                        << "' cannot initialize positions on a species ('"
                        <<this_species->position_initialization_
                        <<"') defined afterwards",
                        LINK_NAMELIST + std::string("#species")
                    );
                }
                this_species->position_initialization_on_species_ = true;
            // HDF5 file where arrays are stored
            } else {
                H5Read f( this_species->position_initialization_ );
                std::vector<std::string> ax = {"weight", "position/x", "position/y", "position/z"};
                for( unsigned int i=0; i<params.nDim_particle+1; i++ ) {
                    std::vector<hsize_t> shape = f.shape( ax[i] );
                    if( shape.size() == 0 ) {
                        ERROR_NAMELIST(
                            "For species '" << species_name 
                            << "': " << ax[i] 
                            << " not found in file " 
                            << this_species->position_initialization_,
                            LINK_NAMELIST + std::string("#species")
                        );
                    }
                    if( this_species->file_position_npart_ == 0 ) {
                        if( shape[0] == 0 ) {
                            ERROR_NAMELIST( 
                                "For species '" << species_name << "': " 
                                << ax[i] << " is empty in file " 
                                << this_species->position_initialization_ ,
                                LINK_NAMELIST + std::string("#species")
                            );
                        }
                        this_species->file_position_npart_ = shape[0];
                    }
                    if( shape.size() > 1 || shape[0] != this_species->file_position_npart_ ) {
                        ERROR_NAMELIST( 
                            "For species '" << species_name 
                            << "': wrong size for " << ax[i] 
                            << " in file " << this_species->position_initialization_,
                            LINK_NAMELIST + std::string("#species") 
                        );
                    }
                }
            }
        }
#ifdef SMILEI_USE_NUMPY
        else if( PyArray_Check( py_pos_init ) ) {
            //Initialize position from this array
            PyArrayObject *np_ret = reinterpret_cast<PyArrayObject *>( py_pos_init );

            //Check dimensions
            unsigned int ndim_local = PyArray_NDIM( np_ret );
            if( ndim_local != 2 ) {
                ERROR_NAMELIST(
                    "For species '" << species_name 
                    << "' Provide a 2-dimensional array in order to init particle position from a numpy array." ,
                    LINK_NAMELIST + std::string("#species") 
                );
            }

            //Check number of coordinates provided
            ndim_local = PyArray_SHAPE( np_ret )[0];
            if( ndim_local != params.nDim_particle + 1 ) {
                ERROR_NAMELIST(
                    "For species '" << species_name
                    << "' position_initializtion must provide a 2-dimensional array with " 
                    <<  params.nDim_particle + 1 << " columns.",
                    LINK_NAMELIST + std::string("#species")  
                );
            }

            //Get number of particles. Do not initialize any more if this is a restart.
            if( !params.restart ) {
                this_species->n_numpy_particles_ =  PyArray_SHAPE( np_ret )[1];
            }
            this_species->position_initialization_array_ = new double[ndim_local*this_species->n_numpy_particles_] ;
            for( unsigned int idim = 0; idim < ndim_local ; idim++ ) {
                for( unsigned int ipart = 0; ipart < ( unsigned int )this_species->n_numpy_particles_; ipart++ ) {
                    this_species->position_initialization_array_[idim*this_species->n_numpy_particles_+ipart] = *( ( double * )PyArray_GETPTR2( np_ret, idim, ipart ) );
                }
            }
        }
#endif
        else {
            ERROR_NAMELIST(
                "For species '" << species_name 
                << "' non valid position_initialization. It must be either a string or a numpy array.",
                LINK_NAMELIST + std::string("#species") 
            );
        }
        Py_DECREF( py_pos_init );

        if( PyTools::extractV( "regular_number", this_species->regular_number_array_, "Species", ispec )){
            if( this_species->position_initialization_ != "regular" ) {
                ERROR_NAMELIST(
                    "regular_number may not be provided if species position_initialization is not set to 'regular'.",
                    LINK_NAMELIST + std::string("#species") 
                );
            }
            if( this_species->regular_number_array_.size() != this_species->nDim_particle ) {
                ERROR_NAMELIST(
                    "Please provide as many regular numbers of particles as there are particle dimensions in the domain ("
                    << this_species->nDim_particle <<").",
                    LINK_NAMELIST + std::string("#species") 
                );
            }
        }

        // Momentum initialisation types
        PyObject *py_mom_init = PyTools::extract_py( "momentum_initialization", "Species", ispec );
        if( PyTools::py2scalar( py_mom_init, this_species->momentum_initialization_ ) ) {
            if(    this_species->momentum_initialization_=="mj"
                || this_species->momentum_initialization_=="maxj" ) {
                this_species->momentum_initialization_="maxwell-juettner";
            }
            // Cold or rectangular
            if(    this_species->momentum_initialization_=="cold"
                || this_species->momentum_initialization_=="rectangular" ) {
                ;
            // Maxwell-juettner
            } else if( this_species->momentum_initialization_=="maxwell-juettner" ) {
                if( this_species->mass_ == 0 ) {
                    ERROR_NAMELIST(
                        "For photon species '" << species_name
                        << "' momentum_initialization cannot be maxwell-juettner",
                        LINK_NAMELIST + std::string("#species")
                    );
                }
            // HDF5 file where arrays are stored
            } else if( this_species->file_position_npart_ > 0 ) {
                H5Read f( this_species->momentum_initialization_ );
                std::vector<std::string> ax = {"momentum/x", "momentum/y", "momentum/z"};
                for( unsigned int i=0; i<ax.size(); i++ ) {
                    std::vector<hsize_t> shape = f.shape( ax[i] );
                    if( shape.size() == 0 ) {
                        ERROR_NAMELIST(
                            "For species '" << species_name 
                            << "': " << ax[i] << " not found in file " 
                            << this_species->momentum_initialization_,
                            LINK_NAMELIST + std::string("#species")
                        );
                    }
                    if( shape.size() > 1 || shape[0] != this_species->file_position_npart_ ) {
                        ERROR_NAMELIST(
                            "For species '" << species_name 
                            << "': wrong shape for " << ax[i] 
                            << " in file " << this_species->momentum_initialization_,
                            LINK_NAMELIST + std::string("#species")
                        );
                    }
                    this_species->file_momentum_npart_ = shape[0];
                }
            // Otherwise, error
            } else {
                ERROR_NAMELIST( 
                    "For species '" << species_name
                    << "' unknown momentum_initialization: "
                    << this_species->momentum_initialization_,
                    LINK_NAMELIST + std::string("#species")
                );
            }
        }
#ifdef SMILEI_USE_NUMPY
        else if( PyArray_Check( py_mom_init ) ) {

            if( !this_species->position_initialization_array_ ) {
                ERROR_NAMELIST( 
                    "For species '" << species_name 
                    << "'. Momentum initialization by a numpy array is only possible if positions are initialized with a numpy array as well. ",
                    LINK_NAMELIST + std::string("#species")
                );
            }

            PyArrayObject *np_ret_mom = reinterpret_cast<PyArrayObject *>( py_mom_init );
            //Check dimensions
            unsigned int ndim_local = PyArray_NDIM( np_ret_mom ) ; //Ok
            if( ndim_local != 2 ) {
                ERROR_NAMELIST(
                    "For species '" << species_name 
                    << "' Provide a 2-dimensional array in order to init particle momentum from a numpy array.",
                    LINK_NAMELIST + std::string("#species")
                );
            }

            //Check number of coordinates provided
            ndim_local =  PyArray_SHAPE( np_ret_mom )[0]; // ok
            if( ndim_local != 3 ) {
                ERROR_NAMELIST( 
                    "For species '" << species_name 
                    << "' momentum_initialization must provide a 2-dimensional array with " 
                    <<  3 << " columns.",
                    LINK_NAMELIST + std::string("#species")
                );
            }

            //Get number of particles
            if( !params.restart && this_species->n_numpy_particles_ != PyArray_SHAPE( np_ret_mom )[1] ) {
                ERROR_NAMELIST(
                    "For species '" << species_name 
                    << "' momentum_initialization must provide as many particles as position_initialization.",
                    LINK_NAMELIST + std::string("#species") 
                );
            }
            this_species->momentum_initialization_array_ = new double[ndim_local*this_species->n_numpy_particles_] ;
            for( unsigned int idim = 0; idim < ndim_local ; idim++ ) {
                for( unsigned int ipart = 0; ipart < ( unsigned int )this_species->n_numpy_particles_; ipart++ ) {
                    this_species->momentum_initialization_array_[idim*this_species->n_numpy_particles_+ipart] = *( ( double * )PyArray_GETPTR2( np_ret_mom, idim, ipart ) );
                }
            }
        }
#endif
        else {
            ERROR_NAMELIST( 
                "For species '" << species_name 
                << "' invalid momentum_initialization. It must be either a string or a numpy array." ,
                LINK_NAMELIST + std::string("#species")
            );
        }
        Py_DECREF( py_mom_init );

        PyTools::extract( "c_part_max", this_species->c_part_max_, "Species", ispec );

        PyTools::extract( "time_frozen", this_species->time_frozen_, "Species", ispec );
        if( this_species->time_frozen_ > 0 && this_species->momentum_initialization_!="cold" ) {
            if( patch->isMaster() ) {
                WARNING( "For species '" << species_name << "' possible conflict between time-frozen & not cold initialization" );
            }
        }
        if (this_species->time_frozen_ > 0) {
            MESSAGE( 2, "> Species frozen until time: " << this_species->time_frozen_ );
        }

        // iteration when the relativistic field initialization is applied, if enabled
        this_species->iter_relativistic_initialization_ = ( int )( this_species->time_frozen_/params.timestep );

        if( !PyTools::extractVV( "boundary_conditions", this_species->boundary_conditions_, "Species", ispec ) ) {
            ERROR_NAMELIST(
                "For species '" << species_name 
                << "', boundary_conditions not defined",
                LINK_NAMELIST + std::string("#species")
            );
        }

        unsigned int number_of_boundaries = (params.geometry=="AMcylindrical") ? 2 : params.nDim_particle;

        if( this_species->boundary_conditions_.size() == 0 ) {
            ERROR_NAMELIST(
                "For species '" << species_name 
                << "', boundary_conditions cannot be empty" ,
                LINK_NAMELIST + std::string("#species")
            );
        } else if( this_species->boundary_conditions_.size() == 1 ) {
            while( this_species->boundary_conditions_.size() < number_of_boundaries ) {
                this_species->boundary_conditions_.push_back( this_species->boundary_conditions_[0] );
            }
        } else if( this_species->boundary_conditions_.size() != number_of_boundaries ) {
            ERROR_NAMELIST( 
                "For species '" << species_name 
                << "', boundary_conditions must be of size "
                << number_of_boundaries <<"." ,
                LINK_NAMELIST + std::string("#species")
            );
        }


        bool has_thermalize = false;
        std::ostringstream t;
        for( unsigned int iDim=0; iDim<number_of_boundaries; iDim++ ) {
            if( this_species->boundary_conditions_[iDim].size() == 1 ) {
                this_species->boundary_conditions_[iDim].push_back( this_species->boundary_conditions_[iDim][0] );
            }
            if( this_species->boundary_conditions_[iDim].size() != 2 ) {
                ERROR_NAMELIST( 
                    "For species '" << species_name 
                    << "', boundary_conditions["<<iDim
                    <<"] must have one or two arguments" ,
                    LINK_NAMELIST + std::string("#species")
                );
            }
            for( unsigned int ii=0; ii<2; ii++ ) {
                if( this_species->boundary_conditions_[iDim][ii] == "thermalize" ) {
                    has_thermalize = true;
                    if( this_species->mass_ == 0 ) {
                        ERROR_NAMELIST( 
                            "For photon species '" << species_name 
                            << "' Thermalizing BCs are not available." ,
                            LINK_NAMELIST + std::string("#species")
                        );
                    }
                } else if( this_species->boundary_conditions_[iDim][ii] == "stop" ) {
                    if( this_species->mass_ == 0 ) {
                        ERROR_NAMELIST( 
                            "For photon species '" << species_name
                             << "' stop BCs are not physical.",
                             LINK_NAMELIST + std::string("#species") 
                         );
                    }
                } else if( this_species->boundary_conditions_[iDim][ii] == "periodic" && params.EM_BCs[iDim][ii] != "periodic" ) {
                    ERROR_NAMELIST( 
                        "For species '" << species_name 
                        << "',  boundary_conditions["<<iDim
                        <<"] cannot be periodic as the EM boundary conditions are not periodic",
                        LINK_NAMELIST + std::string("#species")  
                    );
                }
                t << " " << this_species->boundary_conditions_[iDim][ii];
            }
        }
        if( (params.geometry=="AMcylindrical") && ( this_species->boundary_conditions_[1][1] != "remove" ) && ( this_species->boundary_conditions_[1][1] != "stop" ) ) {
            ERROR_NAMELIST(
                " In AM geometry particle boundary conditions supported in Rmax are 'remove' and 'stop' ",
                LINK_NAMELIST + std::string("#species") 
            );
        }
        if( (params.hasWindow) && (( this_species->boundary_conditions_[0][1] != "remove" ) || ( this_species->boundary_conditions_[0][0] != "remove" ) )) {
            ERROR_NAMELIST(
                " When MovingWindow is activated 'remove' boundary conditions along x is mandatory for all species. ",
                LINK_NAMELIST + std::string("#species")
            );
        }
        MESSAGE( 2, "> Boundary conditions:" << t.str() );

        // for thermalizing BCs on particles check if thermal_boundary_temperature is correctly defined
        bool has_temperature = PyTools::extractV( "thermal_boundary_temperature", this_species->thermal_boundary_temperature_, "Species", ispec ) > 0;
        bool has_velocity    = PyTools::extractV( "thermal_boundary_velocity", this_species->thermal_boundary_velocity_, "Species", ispec ) > 0;
        if( has_thermalize ) {
            if( !has_temperature ) {
                ERROR_NAMELIST( 
                    "For species '" << species_name 
                    << "' thermal_boundary_temperature (thermalizing BC) should be a list of floats",
                    LINK_NAMELIST + std::string("#species") 
                );
            }
            if( !has_velocity ) {
                ERROR_NAMELIST( 
                    "For species '" << species_name 
                    << "' thermal_boundary_velocity (thermalizing BC) should be a list of floats",
                    LINK_NAMELIST + std::string("#species") 
                );
            }
            if( this_species->thermal_boundary_velocity_.size()!=3 ) {
                ERROR_NAMELIST( 
                    "For species '" << species_name 
                    << "' thermal_boundary_velocity (thermalizing BC) should have 3 components",
                    LINK_NAMELIST + std::string("#species") 
                );
            }
            if( this_species->thermal_boundary_temperature_.size()==1 ) {
                WARNING( "For species '" << species_name << "' Using thermal_boundary_temperature[0] in all directions" );
                this_species->thermal_boundary_temperature_.resize( 3 );
                this_species->thermal_boundary_temperature_[1] = this_species->thermal_boundary_temperature_[0];
                this_species->thermal_boundary_temperature_[2] = this_species->thermal_boundary_temperature_[0];
            }

            // Compute the thermal_velocity_ & Momentum for thermalizing bcs
            this_species->thermal_velocity_.resize( 3 );
            this_species->thermal_momentum_.resize( 3 );
            for( unsigned int i=0; i<3; i++ ) {
                this_species->thermal_velocity_[i] = sqrt( 2.*this_species->thermal_boundary_temperature_[i]/this_species->mass_ );
                this_species->thermal_momentum_[i] = this_species->thermal_velocity_[i];
                // Caution: momentum in SMILEI actually correspond to p/m
                if( this_species->thermal_velocity_[i]>0.3 ) {
                    ERROR_NAMELIST(
                        "For species '" << species_name 
                        << "' Thermalizing BCs require non-relativistic thermal_boundary_temperature",
                        LINK_NAMELIST + std::string("#species") 
                    );
                }
            }
        }

        // Manage the ionization parameters
        if( this_species->mass_ > 0 ) {
            this_species->atomic_number_ = 0;
            PyTools::extractOrNone( "atomic_number", this_species->atomic_number_, "Species", ispec );

            this_species->maximum_charge_state_ = 0;
            PyTools::extract( "maximum_charge_state", this_species->maximum_charge_state_, "Species", ispec);

            std::string model;
            PyTools::extract( "ionization_model", model, "Species", ispec );
            if( model!="none" ) {

                this_species->ionization_model = model;

                if( this_species->particles->is_test ) {
                    ERROR_NAMELIST( 
                        "For species '" << species_name 
                        << ": cannot ionize test species",
                        LINK_NAMELIST + std::string("#species") );
                }

                if( ( this_species->atomic_number_==0 )&&( this_species->maximum_charge_state_==0 ) ) {
                    ERROR_NAMELIST( 
                        "For species '" << species_name 
                        << ": undefined atomic_number & maximum_charge_state (required for ionization)",
                        LINK_NAMELIST + std::string("#species") );
                }

                if( model == "tunnel" ){
                    if (params.Laser_Envelope_model){
                        ERROR_NAMELIST("An envelope is present, so tunnel_envelope or tunnel_envelope_averaged ionization model should be selected for species "<<species_name,
                        LINK_NAMELIST + std::string("#species"));
                    }
                }else if ( model == "tunnel_envelope_averaged" ){

                    if (!params.Laser_Envelope_model) {
                        ERROR_NAMELIST("An envelope ionization model has been selected but no envelope is present",
                        LINK_NAMELIST + std::string("#laser-envelope-model"));
                    }

                }else if( model == "from_rate" ) {

                    if( this_species->maximum_charge_state_ == 0 ) {
                        this_species->maximum_charge_state_ = this_species->atomic_number_;
                        WARNING( "For species '" << species_name << ": ionization 'from_rate' is used with maximum_charge_state = "<<this_species->maximum_charge_state_ << " taken from atomic_number" );
                    }
                    this_species->ionization_rate_ = PyTools::extract_py( "ionization_rate", "Species", ispec );
                    if( this_species->ionization_rate_ == Py_None ) {
                        ERROR_NAMELIST( "For species '" << species_name << " ionization 'from_rate' requires 'ionization_rate' ",
                        LINK_NAMELIST + std::string("#species") );
                    } else {
#ifdef SMILEI_USE_NUMPY
                        // Test the ionization_rate function with temporary, "fake" particles
                        std::ostringstream name( "" );
                        name << " ionization_rate:";
                        double *dummy = NULL;
                        ParticleData test( params.nDim_particle, this_species->ionization_rate_, name.str(), dummy );
#else
                        ERROR_NAMELIST( "For species '" << species_name << " ionization 'from_rate' requires Numpy",
                        LINK_NAMELIST + std::string("#species") );
#endif
                    }

                } else if( model != "none" ) {
                    ERROR_NAMELIST( "For species " << species_name << ": unknown ionization model `" << model,
                    LINK_NAMELIST + std::string("#species") );
                }

                if( params.vectorization_mode != "off" ) {
                    WARNING( "Performances of advanced physical processes which generates new particles could be degraded for the moment!" );
                    WARNING( "\t The improvement of their integration in vectorized algorithms is in progress." );
                }

                PyTools::extract( "ionization_electrons", this_species->ionization_electrons, "Species", ispec );

            }

        }

        // Extract if the species is relativistic and needs ad hoc fields initialization
        bool relativistic_field_initialization = false;
        PyTools::extract( "relativistic_field_initialization", relativistic_field_initialization, "Species", ispec );
        this_species->relativistic_field_initialization_ = relativistic_field_initialization;



        // Species geometry
        // ----------------

        // Density
        bool ok1, ok2;
        PyObject * profile1 = nullptr;
        if( this_species->position_initialization_array_ == NULL
         && this_species->file_position_npart_ == 0 ) {
            //These quantities are disregarded if positioning of the species is directly specified by the user
            // Matter particles
            ok1 = PyTools::extract_pyProfile( "number_density", profile1, "Species", ispec );
            ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "Species", ispec );
            if( this_species->mass_ > 0 ) {
                if( ok1 &&  ok2 ) {
                    ERROR_NAMELIST( "For species '" << species_name << "', cannot define both `number_density ` and `charge_density`.",
                    LINK_NAMELIST + std::string("#species") );
                }
                if( !ok1 && !ok2 ) {
                    ERROR_NAMELIST( "For species '" << species_name << "', must define `number_density ` or `charge_density`.",
                    LINK_NAMELIST + std::string("#species") );
                }
                if( ok1 ) {
                    this_species->density_profile_type_ = "nb";
                }
                if( ok2 ) {
                    this_species->density_profile_type_ = "charge";
                }
                //MESSAGE(this_species->density_profile_type_);
            }
            // Photons
            else if( this_species->mass_ == 0 ) {
                if( ok2 ) {
                    ERROR_NAMELIST( "For photon species '" << species_name << "', `charge_density` has no meaning."
                            << "You must use `number_density`.",
                            LINK_NAMELIST + std::string("#species") );
                }
                if( !ok1 ) {
                    ERROR_NAMELIST( "For photon species '" << species_name << "', must define `number_density`.",
                    LINK_NAMELIST + std::string("#species") );
                }
                this_species->density_profile_type_ = "nb";
            }

            this_species->density_profile_ = new Profile( profile1, params.nDim_field, Tools::merge( this_species->density_profile_type_, "_density ", species_name ), params, true, true );
            MESSAGE(2, "> Density profile: " << this_species->density_profile_->getInfo());

            // Number of particles per cell
            if( !PyTools::extract_pyProfile( "particles_per_cell", profile1, "Species", ispec ) ) {
                ERROR_NAMELIST( "For species '" << species_name << "', particles_per_cell not found or not understood",
                LINK_NAMELIST + std::string("#species") );
            }
            this_species->particles_per_cell_profile_ = new Profile( profile1, params.nDim_field, Tools::merge( "particles_per_cell ", species_name ), params, true, true );
        } else {
            if( PyTools::extract_pyProfile( "particles_per_cell", profile1, "Species", ispec ) ) {
                ERROR_NAMELIST( "For species '" << species_name << "', cannot define both `particles_per_cell` and  `position_initialization` array.",
                LINK_NAMELIST + std::string("#species") );
            }
            ok1 = PyTools::extract_pyProfile( "number_density", profile1, "Species", ispec );
            ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "Species", ispec );
            if( ok1 ||  ok2 ) {
                ERROR_NAMELIST( "For species '" << species_name << "', cannot define both `density` and `position_initialization` array.",
                LINK_NAMELIST + std::string("#species") );
            }
        }

        // Charge
        if( !PyTools::extract_pyProfile( "charge", profile1, "Species", ispec ) ) {
            ERROR_NAMELIST( "For species '" << species_name << "', charge not found or not understood",
            LINK_NAMELIST + std::string("#species") );
        }
        this_species->charge_profile_ = new Profile( profile1, params.nDim_field, Tools::merge( "charge ", species_name ), params, true, true );

        std::vector<PyObject *> prof;
        if( this_species->momentum_initialization_array_ == NULL
         && this_species->file_momentum_npart_ == 0 ) {
            // Mean velocity
            if( PyTools::extract_1or3Profiles( "mean_velocity", "Species", ispec, prof ) ) {
                this_species->velocity_profile_[0] = new Profile( prof[0], params.nDim_field, Tools::merge( "mean_velocity[0] ", species_name ), params, true, true );
                this_species->velocity_profile_[1] = new Profile( prof[1], params.nDim_field, Tools::merge( "mean_velocity[1] ", species_name ), params, true, true );
                this_species->velocity_profile_[2] = new Profile( prof[2], params.nDim_field, Tools::merge( "mean_velocity[2] ", species_name ), params, true, true );
            }
            // Temperature
            if( PyTools::extract_1or3Profiles( "temperature", "Species", ispec, prof ) ) {
                this_species->temperature_profile_[0] = new Profile( prof[0], params.nDim_field, Tools::merge( "temperature[0] ", species_name ), params, true, true );
                this_species->temperature_profile_[1] = new Profile( prof[1], params.nDim_field, Tools::merge( "temperature[1] ", species_name ), params, true, true );
                this_species->temperature_profile_[2] = new Profile( prof[2], params.nDim_field, Tools::merge( "temperature[2] ", species_name ), params, true, true );
            }
        } else {
            ok1 = PyTools::extract_1or3Profiles( "mean_velocity", "Species", ispec, prof ) ;
            ok2 = PyTools::extract_1or3Profiles( "temperature", "Species", ispec, prof ) ;
            if( ok1 ) {
                ERROR_NAMELIST( "For species '" << species_name << "', cannot define both `mean_velocity` and `momentum_initialization` array.",
                LINK_NAMELIST + std::string("#species") );
            }
            if( ok2 ) {
                ERROR_NAMELIST( "For species '" << species_name << "', cannot define both `temperature` and `momentum_initialization` array.",
                LINK_NAMELIST + std::string("#species") );
            }
        }


        // Get info about tracking
        unsigned int ntrack = PyTools::nComponents( "DiagTrackParticles" );
        this_species->particles->tracked = false;
        for( unsigned int itrack=0; itrack<ntrack; itrack++ ) {
            std::string track_species;
            PyTools::extract( "species", track_species, "DiagTrackParticles", itrack );
            if( track_species==species_name ) {
                if( this_species->particles->tracked ) {
                    ERROR_NAMELIST( "In this version, species '" << species_name << "' cannot be tracked by two DiagTrackParticles",
                    LINK_NAMELIST + std::string("#species") );
                }
                this_species->particles->tracked  = true;
            }
        }

        // Extract test Species flag
        PyTools::extract( "is_test", this_species->particles->is_test, "Species", ispec );

        // Verify they don't ionize
        if( this_species->ionization_model!="none" && this_species->particles->is_test ) {
            ERROR_NAMELIST( "For species '" << species_name << "' test & ionized is currently impossible",
            LINK_NAMELIST + std::string("#species") );
        }

        return this_species;
    } // End Species* create()


    // Method to clone a species from an existing one
    // Note that this must be only called from cloneVector, because additional init is needed
    static Species *clone( Species *species, Params &params, Patch *patch )
    {

        // Create new species object
        Species *new_species = NULL;


        // Type of species
        if ( params.vectorization_mode == "off" ) {
            new_species = new Species( params, patch );
        }
#ifdef _VECTO
        else if( params.vectorization_mode == "on" ) {
            new_species = new SpeciesV( params, patch );
        } else if( params.vectorization_mode == "adaptive" ) {
            new_species = new SpeciesVAdaptive( params, patch );
        } else if( params.vectorization_mode == "adaptive_mixed_sort" ) {
            new_species = new SpeciesVAdaptiveMixedSort( params, patch );
        }
#endif

        // Copy members
        new_species->name_                                     = species->name_;
        new_species->pusher_name_                              = species->pusher_name_;
        new_species->radiation_model_                          = species->radiation_model_;
        new_species->radiation_photon_species                  = species->radiation_photon_species;
        new_species->radiation_photon_sampling_                = species->radiation_photon_sampling_;
        new_species->radiation_photon_gamma_threshold_         = species->radiation_photon_gamma_threshold_;
        new_species->photon_species_                           = species->photon_species_;
        new_species->species_number_                           = species->species_number_;
        new_species->position_initialization_on_species_       = species->position_initialization_on_species_;
        new_species->position_initialization_on_species_index  = species->position_initialization_on_species_index;
        new_species->position_initialization_                  = species->position_initialization_;
        new_species->position_initialization_array_            = species->position_initialization_array_;
        new_species->file_position_npart_                      = species->file_position_npart_;
        new_species->file_momentum_npart_                      = species->file_momentum_npart_;
        new_species->regular_number_array_                     = species->regular_number_array_;
        new_species->n_numpy_particles_                        = species->n_numpy_particles_            ;
        new_species->momentum_initialization_                  = species->momentum_initialization_;
        new_species->momentum_initialization_array_            = species->momentum_initialization_array_;
        new_species->c_part_max_                               = species->c_part_max_;
        new_species->mass_                                     = species->mass_;
        new_species->time_frozen_                              = species->time_frozen_;
        new_species->radiating_                                = species->radiating_;
        new_species->relativistic_field_initialization_        = species->relativistic_field_initialization_;
        new_species->iter_relativistic_initialization_         = species->iter_relativistic_initialization_;
        new_species->boundary_conditions_                      = species->boundary_conditions_;
        new_species->thermal_boundary_temperature_             = species->thermal_boundary_temperature_;
        new_species->thermal_boundary_velocity_                = species->thermal_boundary_velocity_;
        new_species->thermal_velocity_                          = species->thermal_velocity_;
        new_species->thermal_momentum_                          = species->thermal_momentum_;
        new_species->atomic_number_                            = species->atomic_number_;
        new_species->maximum_charge_state_                     = species->maximum_charge_state_;
        new_species->ionization_rate_                          = species->ionization_rate_;
        if( new_species->ionization_rate_!=Py_None ) {
            Py_INCREF( new_species->ionization_rate_ );
        }
        new_species->ionization_model                         = species->ionization_model;
        new_species->density_profile_type_                       = species->density_profile_type_;
        new_species->vectorized_operators                     = species->vectorized_operators;
        new_species->merging_method_                          = species->merging_method_;
        new_species->has_merging_                             = species->has_merging_;
        new_species->merging_time_selection_                  = species->merging_time_selection_;
        new_species->merge_log_scale_                         = species->merge_log_scale_;
        new_species->merge_min_momentum_log_scale_            = species->merge_min_momentum_log_scale_;
        new_species->merge_min_particles_per_cell_            = species->merge_min_particles_per_cell_;
        new_species->merge_min_packet_size_                   = species->merge_min_packet_size_;
        new_species->merge_max_packet_size_                   = species->merge_max_packet_size_;
        new_species->merge_accumulation_correction_           = species->merge_accumulation_correction_;
        new_species->merge_momentum_cell_size_[0]             = species->merge_momentum_cell_size_[0];
        new_species->merge_momentum_cell_size_[1]             = species->merge_momentum_cell_size_[1];
        new_species->merge_momentum_cell_size_[2]             = species->merge_momentum_cell_size_[2];
        new_species->merge_min_momentum_cell_length_[0]       = species->merge_min_momentum_cell_length_[0];
        new_species->merge_min_momentum_cell_length_[1]       = species->merge_min_momentum_cell_length_[1];
        new_species->merge_min_momentum_cell_length_[2]       = species->merge_min_momentum_cell_length_[2];


        new_species->charge_profile_                            = new Profile( species->charge_profile_ );
        if( species->density_profile_ ) {
            new_species->density_profile_                     = new Profile( species->density_profile_ );
            new_species->particles_per_cell_profile_          = new Profile( species->particles_per_cell_profile_ );
        }
        new_species->velocity_profile_.resize( 3 );
        new_species->temperature_profile_.resize( 3 );

        if( species->velocity_profile_[0] ) {
            new_species->velocity_profile_[0]                 = new Profile( species->velocity_profile_[0] );
            new_species->velocity_profile_[1]                 = new Profile( species->velocity_profile_[1] );
            new_species->velocity_profile_[2]                 = new Profile( species->velocity_profile_[2] );
        }
        if( species->temperature_profile_[0] ) {
            new_species->temperature_profile_[0]              = new Profile( species->temperature_profile_[0] );
            new_species->temperature_profile_[1]              = new Profile( species->temperature_profile_[1] );
            new_species->temperature_profile_[2]              = new Profile( species->temperature_profile_[2] );
        }
        new_species->max_charge_                               = species->max_charge_;
        new_species->tracking_diagnostic                      = species->tracking_diagnostic;

        if( new_species->mass_==0 ) {
            new_species->multiphoton_Breit_Wheeler_[0]         = species->multiphoton_Breit_Wheeler_[0];
            new_species->multiphoton_Breit_Wheeler_[1]         = species->multiphoton_Breit_Wheeler_[1];
            new_species->mBW_pair_creation_sampling_[0]        = species->mBW_pair_creation_sampling_[0];
            new_species->mBW_pair_creation_sampling_[1]        = species->mBW_pair_creation_sampling_[1];
        }

        new_species->particles->is_test                       = species->particles->is_test;
        new_species->particles->tracked                       = species->particles->tracked;
        new_species->particles->isQuantumParameter            = species->particles->isQuantumParameter;
        new_species->particles->isMonteCarlo                  = species->particles->isMonteCarlo;

        return new_species;
    } // End Species* clone()


    static void createVector( Params &params, Patch *patch )
    {
        // read number of species from python namelist
        unsigned int tot_species_number = PyTools::nComponents( "Species" );
        if (tot_species_number > 0) {
            TITLE("Initializing species");
        }

        // Create all species
        for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {
            Species *this_species = SpeciesFactory::create( params, ispec, patch );
            // Verify the new species does not have the same name as a previous one
            for( unsigned int i = 0; i < ispec; i++ ) {
                if( this_species->name_ == patch->vecSpecies[i]->name_ ) {
                    ERROR_NAMELIST("Two species cannot have the same name `"<<this_species->name_<<"`",
                    LINK_NAMELIST + std::string("#species"));
                }
            }
            // Put the newly created species in the vector of species
            patch->vecSpecies.push_back( this_species );
        }

        // Initialize particles & operators
        for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {
            patch->vecSpecies[ispec]->initParticles( params, patch );
            patch->vecSpecies[ispec]->initOperators( params, patch );
        }

        // Loop species to find related species
        for( unsigned int ispec1 = 0; ispec1<patch->vecSpecies.size(); ispec1++ ) {

            // Ionizable species
            if( patch->vecSpecies[ispec1]->Ionize ) {
                // Loop all other species
                for( unsigned int ispec2 = 0; ispec2<patch->vecSpecies.size(); ispec2++ ) {
                    if( patch->vecSpecies[ispec1]->ionization_electrons == patch->vecSpecies[ispec2]->name_ ) {
                        if( ispec1==ispec2 ) {
                            ERROR_NAMELIST( "For species '"<<patch->vecSpecies[ispec1]->name_<<"' ionization_electrons must be a distinct species",
                            LINK_NAMELIST + std::string("#species") );
                        }
                        if( patch->vecSpecies[ispec2]->mass_!=1 ) {
                            ERROR_NAMELIST( "For species '"<<patch->vecSpecies[ispec1]->name_<<"' ionization_electrons must be a species with mass==1",
                            LINK_NAMELIST + std::string("#species") );
                        }
                        patch->vecSpecies[ispec1]->electron_species_index = ispec2;
                        patch->vecSpecies[ispec1]->electron_species = patch->vecSpecies[ispec2];

                        int max_eon_number =
                            patch->vecSpecies[ispec1]->getNbrOfParticles()
                            * ( patch->vecSpecies[ispec1]->atomic_number_ || patch->vecSpecies[ispec1]->maximum_charge_state_ );
                        patch->vecSpecies[ispec1]->Ionize->new_electrons.initializeReserve(
                            max_eon_number, *patch->vecSpecies[ispec1]->electron_species->particles
                        );
                        break;
                    }
                }
                if( patch->vecSpecies[ispec1]->electron_species_index==-1 ) {
                    ERROR_NAMELIST( "For species '"<<patch->vecSpecies[ispec1]->name_
                    <<"' ionization_electrons named " << patch->vecSpecies[ispec1]->ionization_electrons 
                    << " could not be found",
                    LINK_NAMELIST + std::string("#species") );
                }
            }

            // Radiating species
            if( patch->vecSpecies[ispec1]->Radiate ) {
                // No emission of discrete photon, only scalar diagnostics are updated
                if( patch->vecSpecies[ispec1]->radiation_photon_species.empty() ) {
                    patch->vecSpecies[ispec1]->photon_species_index = -1;
                    patch->vecSpecies[ispec1]->photon_species_ = nullptr;
                    patch->vecSpecies[ispec1]->radiated_photons_ = nullptr;
                }
                // Else, there will be emission of macro-photons.
                else {
                    unsigned int ispec2 = 0;
                    for( ispec2 = 0; ispec2<patch->vecSpecies.size(); ispec2++ ) {
                        if( patch->vecSpecies[ispec1]->radiation_photon_species == patch->vecSpecies[ispec2]->name_ ) {
                            if( ispec1==ispec2 ) {
                                ERROR_NAMELIST( "For species '"<<patch->vecSpecies[ispec1]->name_
                                <<"' radiation_photon_species must be a distinct photon species",
                                LINK_NAMELIST + std::string("#species") );
                            }
                            if( patch->vecSpecies[ispec2]->mass_!=0 ) {
                                ERROR_NAMELIST( "For species '"<<patch->vecSpecies[ispec1]->name_
                                <<"' radiation_photon_species must be a photon species with mass==0",
                                LINK_NAMELIST + std::string("#species") );
                            }
                            patch->vecSpecies[ispec1]->photon_species_index = ispec2;
                            patch->vecSpecies[ispec1]->photon_species_ = patch->vecSpecies[ispec2];
                            patch->vecSpecies[ispec1]->radiated_photons_ = ParticlesFactory::create( params );
                            patch->vecSpecies[ispec1]->radiated_photons_->initializeReserve(
                                patch->vecSpecies[ispec1]->getNbrOfParticles(),
                                *patch->vecSpecies[ispec1]->photon_species_->particles
                            );
                            // patch->vecSpecies[ispec1]->radiated_photons_->initialize(
                            //     0,params.keep_position_old);
                            //patch->vecSpecies[ispec1]->radiated_photons_->reserve( 10, params.nDim_particle );
                            break;
                        }
                    }
                    if( ispec2 == patch->vecSpecies.size() ) {
                        ERROR_NAMELIST( "Species '" << patch->vecSpecies[ispec1]->radiation_photon_species 
                        << "' does not exist.",
                        LINK_NAMELIST + std::string("#species") )
                    }
                }
            }

            // Breit-Wheeler species
            if( patch->vecSpecies[ispec1]->Multiphoton_Breit_Wheeler_process ) {
                unsigned int ispec2;
                for( int k=0; k<2; k++ ) {
                    ispec2 = 0;
                    while( ispec2<patch->vecSpecies.size()) {
                        // We look for the pair species multiphoton_Breit_Wheeler_[k]
                        if( patch->vecSpecies[ispec1]->multiphoton_Breit_Wheeler_[k] == patch->vecSpecies[ispec2]->name_ ) {
                            if( ispec1==ispec2 ) {
                                ERROR_NAMELIST( "For species '" << patch->vecSpecies[ispec1]->name_
                                       << "' pair species must be a distinct particle species",
                                       LINK_NAMELIST + std::string("#species") );
                            }
                            if( patch->vecSpecies[ispec2]->mass_ != 1 ) {
                                ERROR_NAMELIST( "For species '"<<patch->vecSpecies[ispec1]->name_
                                  <<"' pair species must be an electron and positron species (mass = 1)",
                                  LINK_NAMELIST + std::string("#species") );
                            }
                            patch->vecSpecies[ispec1]->mBW_pair_species_index[k] = ispec2;
                            patch->vecSpecies[ispec1]->mBW_pair_species[k] = patch->vecSpecies[ispec2];
                            patch->vecSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].initializeReserve(
                                patch->vecSpecies[ispec1]->getNbrOfParticles(),
                                *patch->vecSpecies[ispec1]->mBW_pair_species[k]->particles
                            );
                            ispec2 = patch->vecSpecies.size() + 1;
                        }
                        ispec2++ ;
                    }
                    // This means that one of the pair species has not been fould
                    if( ispec2 == patch->vecSpecies.size() ) {
                        ERROR_NAMELIST( "In Species `" << patch->vecSpecies[ispec1]->name_ << "`"
                           << " the pair species `" << patch->vecSpecies[ispec1]->multiphoton_Breit_Wheeler_[k]
                           << "` does not exist.",
                           LINK_NAMELIST + std::string("#species") )
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    //! Method to clone the whole vector of species
    // -----------------------------------------------------------------------------------------------------------------
    static void cloneVector( std::vector<Species *> vector_species, Params &params, Patch *patch, bool with_particles = true )
    {
        patch->vecSpecies.resize( 0 );

        // Clone all species
        for( unsigned int i = 0; i < vector_species.size(); i++ ) {
            Species *new_species = SpeciesFactory::clone( vector_species[i], params, patch );
            patch->vecSpecies.push_back( new_species );
        }

        // Initialize particles & operators
        for( unsigned int i = 0; i < vector_species.size(); i++ ) {
            patch->vecSpecies[i]->initParticles( params, patch, with_particles, vector_species[i]->particles );
            patch->vecSpecies[i]->initOperators( params, patch );
        }

        // Ionization
        for( unsigned int i=0; i<patch->vecSpecies.size(); i++ ) {
            if( patch->vecSpecies[i]->Ionize ) {
                patch->vecSpecies[i]->electron_species_index = vector_species[i]->electron_species_index;
                patch->vecSpecies[i]->electron_species = patch->vecSpecies[patch->vecSpecies[i]->electron_species_index];
                patch->vecSpecies[i]->Ionize->new_electrons.tracked = patch->vecSpecies[i]->electron_species->particles->tracked;
                patch->vecSpecies[i]->Ionize->new_electrons.isQuantumParameter = patch->vecSpecies[i]->electron_species->particles->isQuantumParameter;
                patch->vecSpecies[i]->Ionize->new_electrons.isMonteCarlo = patch->vecSpecies[i]->electron_species->particles->isMonteCarlo;
                patch->vecSpecies[i]->Ionize->new_electrons.initialize( 0, params.nDim_particle, params.keep_position_old );
            }
        }

        // Synchrotron-like radiation
        for( unsigned int i=0; i<patch->vecSpecies.size(); i++ ) {
            if( patch->vecSpecies[i]->Radiate ) {
                patch->vecSpecies[i]->radiation_photon_species = vector_species[i]->radiation_photon_species;
                patch->vecSpecies[i]->photon_species_index = vector_species[i]->photon_species_index;
                // Photon emission activated:
                if( vector_species[i]->photon_species_ ) {
                    patch->vecSpecies[i]->photon_species_ = patch->vecSpecies[patch->vecSpecies[i]->photon_species_index];
                    //patch->vecSpecies[i]->Radiate->new_photons_.initialize(patch->vecSpecies[i]->getNbrOfParticles(),
                    //                                               params.nDim_particle );
                    patch->vecSpecies[i]->radiated_photons_ = ParticlesFactory::create( params );
                    patch->vecSpecies[i]->radiated_photons_->tracked = patch->vecSpecies[i]->photon_species_->particles->tracked;
                    patch->vecSpecies[i]->radiated_photons_->isQuantumParameter = patch->vecSpecies[i]->photon_species_->particles->isQuantumParameter;
                    patch->vecSpecies[i]->radiated_photons_->isMonteCarlo = patch->vecSpecies[i]->photon_species_->particles->isMonteCarlo;
                    patch->vecSpecies[i]->radiated_photons_->initialize( 0, params.nDim_particle, params.keep_position_old );
                    //patch->vecSpecies[i]->radiated_photons_->reserve( 10, params.nDim_particle );
                    // patch->vecSpecies[i]->radiated_photons_->initializeReserve( 0, *patch->vecSpecies[i]->photon_species_->particles );
                } else {
                    patch->vecSpecies[i]->photon_species_ = NULL;
                    patch->vecSpecies[i]->radiated_photons_ = NULL;
                }
            }
        }

        // multiphoton Breit-Wheeler
        for( unsigned int i=0; i<patch->vecSpecies.size(); i++ ) {
            if( patch->vecSpecies[i]->Multiphoton_Breit_Wheeler_process ) {
                // Loop on pairs
                for( int k=0; k<2; k++ ) {
                    patch->vecSpecies[i]->multiphoton_Breit_Wheeler_[k] = vector_species[i]->multiphoton_Breit_Wheeler_[k];
                    patch->vecSpecies[i]->mBW_pair_species_index[k] = vector_species[i]->mBW_pair_species_index[k];
                    patch->vecSpecies[i]->mBW_pair_species[k] = patch->vecSpecies[patch->vecSpecies[i]->mBW_pair_species_index[k]];
                    patch->vecSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].tracked = patch->vecSpecies[i]->mBW_pair_species[k]->particles->tracked;
                    patch->vecSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].isQuantumParameter = patch->vecSpecies[i]->mBW_pair_species[k]->particles->isQuantumParameter;
                    patch->vecSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].isMonteCarlo = patch->vecSpecies[i]->mBW_pair_species[k]->particles->isMonteCarlo;
                    patch->vecSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].initialize(
                        0, params.nDim_particle, params.keep_position_old );
                }
            } else {
                patch->vecSpecies[i]->mBW_pair_species[0] = NULL;
                patch->vecSpecies[i]->mBW_pair_species[1] = NULL;
            }
        }
    }

};

#endif
