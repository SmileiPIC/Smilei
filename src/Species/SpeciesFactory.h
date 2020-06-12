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
#include "SpeciesNorm.h"

#ifdef _VECTO
#include "SpeciesNormV.h"
#include "SpeciesV.h"
#include "SpeciesVAdaptiveMixedSort.h"
#include "SpeciesVAdaptive.h"
#include "SpeciesNormV.h"
#endif

#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "PartBoundCond.h"

#include "Params.h"
#include "Patch.h"

#include "ParticleData.h"

#include "ParticleCreator.h"

#include "Tools.h"
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
            ERROR("For species #" << ispec << ", name cannot be only 1 character");
        }

        if( species_name.substr(0,2) == "m_" ) {
            ERROR("For species #" << ispec << ", name cannot start  with `m_`");
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
                if( ( params.vectorization_mode == "off" ) && !params.cell_sorting ) {
                    this_species = new SpeciesNorm( params, patch );
                }

#ifdef _VECTO
                else if( ( params.vectorization_mode == "on" ) || params.cell_sorting ) {
                    this_species = new SpeciesNormV( params, patch );
                } else if( params.vectorization_mode == "adaptive_mixed_sort" ) {
                    this_species = new SpeciesVAdaptiveMixedSort( params, patch );
                } else if( params.vectorization_mode == "adaptive" ) {
                    this_species = new SpeciesVAdaptive( params, patch );
                }
#endif
            } else {
                ERROR( "For species `" << species_name << "`, pusher must be 'boris', 'borisnr', 'vay', 'higueracary', 'ponderomotive_boris'" );
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
                ERROR( "For species `" << species_name
                       << " radiation_model must be 'none',"
                       << " 'Landau-Lifshitz' ('ll'),"
                       << " 'corrected-Landau-Lifshitz' ('cll'),"
                       << " 'Niel' ('niel') or 'Monte-Carlo' ('mc')" );
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
                ERROR( "For species `" << species_name
                       << "` radiation_model `"
                       << radiation_model
                       << "` is not compatible with pusher "
                       << pusher );

            }

        }

        // Photon species
        else if( mass == 0 ) {
            if( ( params.vectorization_mode == "off" ) && !params.cell_sorting ) {
                this_species = new SpeciesNorm( params, patch );
            }
#ifdef _VECTO
            else if( ( params.vectorization_mode == "on" ) || params.cell_sorting ) {
                this_species = new SpeciesNormV( params, patch );
            } else if( params.vectorization_mode == "adaptive_mixed_sort" ) {
                this_species = new SpeciesVAdaptiveMixedSort( params, patch );
            } else if( params.vectorization_mode == "adaptive" ) {
                this_species = new SpeciesVAdaptive( params, patch );
            }
#endif
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
                        ERROR( " The radiation photon species is not specified." )
                    }

                    // Number of photons emitted per Monte-Carlo event
                    PyTools::extract( "radiation_photon_sampling",
                        this_species->radiation_photon_sampling_, "Species", ispec );
                    if( this_species->radiation_photon_sampling_ < 1 ) {
                        ERROR( "For species '" << species_name
                               << "' radiation_photon_sampling should be > 1" );
                    }
                    MESSAGE( 3, "| Number of macro-photons emitted per MC event: "
                             << this_species->radiation_photon_sampling_ );

                    // Photon energy threshold
                    PyTools::extract( "radiation_photon_gamma_threshold",
                        this_species->radiation_photon_gamma_threshold_, "Species", ispec );
 
                    MESSAGE( 3, "| Photon energy threshold for macro-photon emission: "
                             << this_species->radiation_photon_gamma_threshold_ );
                } else {
                    MESSAGE( 3, "| Macro-photon emission not activated" );
                }

            }
        }

        // Multiphoton Breit-Wheeler
        if( mass == 0 ) {
            // If this_species->multiphoton_Breit_Wheeler
            if( PyTools::extractV( "multiphoton_Breit_Wheeler", this_species->multiphoton_Breit_Wheeler_, "Species", ispec ) ) {
                // If one of the species is empty
                if( this_species->multiphoton_Breit_Wheeler_[1].empty() || this_species->multiphoton_Breit_Wheeler_[0].empty() ) {
                    ERROR( "For species '" << species_name
                           << "' multiphoton_Breit_Wheeler can not be empty,"
                           << " select electron and positron species." );
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
            ERROR( "In Species " << this_species->name_
                   << ": merging method not valid, must be `vranic_spherical`, `vranic_cartesian` or `none`" );
        }

        // if( params.vectorization_mode == "off" && this_species->merging_method_ != "none" ) {
        //     ERROR( "In Species " << this_species->name_
        //            << ": particle merging only available with `vectorization_mode` = `on` or `adaptive`" );
        // }

        if ( this_species->merging_method_ != "none" ) {

            if (!params.cell_sorting && !this_species->vectorized_operators) {
                ERROR( "In Species " << this_species->name_
                       << ": merging required cell sorting to be "
                       << "activated (`cell_sorting = True` in the mains or vectorization on).");
            }

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
                ERROR( "In Species " << this_species->name_
                       << ": minimum number of particles per merging packet "
                       << "(`merge_min_packet_size`)"
                       << "must be above or equal to 4.");
            }
            if (this_species->merge_min_packet_size_ < 4 && this_species->mass_ == 0)
            {
                ERROR( "In Species " << this_species->name_
                       << " of type photon"
                       << ": minimum number of particles per merging packet "
                       << "(`merge_min_packet_size`)"
                       << "must be above or equal to 4.");
            }
            
            // Maximum particle number per packet to merge
            PyTools::extract( "merge_max_packet_size", this_species->merge_max_packet_size_ , "Species", ispec );
            if (this_species->merge_max_packet_size_ < 4 && this_species->mass_ > 0)
            {
                ERROR( "In Species " << this_species->name_
                       << ": maximum number of particles per merging packet "
                       << "(`merge_max_packet_size`)"
                       << "must be above or equal to 4.");
            }
            if (this_species->merge_max_packet_size_ < 4 && this_species->mass_ == 0)
            {
                ERROR( "In Species " << this_species->name_
                       << " of type photon"
                       << ": maximum number of particles per merging packet "
                       << "(`merge_max_packet_size`)"
                       << "must be above or equal to 4.");
            }
            if (this_species->merge_max_packet_size_ < this_species->merge_min_packet_size_) {
                ERROR( "In Species " << this_species->name_
                       << ": maximum number of particles per merging packet "
                       << "(`merge_max_packet_size`)"
                       << "must be below or equal to the minimum particle number"
                       << " per merging packet (`merge_min_packet_size`)");
            }
            
            // Minimum momentum cell length for the momentum discretization
            if( PyTools::extractV( "merge_min_momentum_cell_length",
                                  this_species->merge_min_momentum_cell_length_ ,
                                  "Species", ispec ) ) {
                for (unsigned int i = 0 ; i < 3 ; i++) {
                    if (this_species->merge_min_momentum_cell_length_[i] <= 0) {
                        ERROR( "In Species " << this_species->name_
                                 << ": The minimal momentum cell length "
                                 << "(`merge_min_particles_per_cell`)"
                                 << " must be above 0 ("
                                 << this_species->merge_min_momentum_cell_length_[i]
                                 << ")");
                    }
                }
            } else {
                ERROR( "In Species " << this_species->name_ << ": merge_min_momentum_cell_length should be a list of floats" );
            }

            // Read and check the threshold on the number of particles per cell
            PyTools::extract( "merge_min_particles_per_cell", this_species->merge_min_particles_per_cell_ , "Species", ispec );
            if (this_species->merge_min_particles_per_cell_ < 4) {
                ERROR( "In Species " << this_species->name_
                       << ": The threshold on the number of particles per cell "
                       << "(`merge_min_particles_per_cell`)"
                       << "must be above or equal to 4");
            }
            
            // Read flag to activate the accumulation correction
            PyTools::extract( "merge_accumulation_correction", this_species->merge_accumulation_correction_ , "Species", ispec );

            // Momentum cell discretization
            if( PyTools::extractV( "merge_momentum_cell_size",
                                  this_species->merge_momentum_cell_size_ ,
                                  "Species", ispec ) ) {
                for (unsigned int i = 0 ; i < 3 ; i++) {
                    if (this_species->merge_momentum_cell_size_[i] <= 0) {
                        ERROR( "In Species " << this_species->name_
                               << ": The momentum cell discretization can not be equal or below 0 "
                               << "(`merge_momentum_cell_size_`).");
                    }
                }
            } else {
                ERROR( "In Species " << this_species->name_ << ": merge_momentum_cell_size should be a list of integers" );
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
                ERROR( "In Species " << this_species->name_
                       << ": discretization_scale should be `linear` or `log`");
            }
            
            // Minimum momentum in log scale
            PyTools::extract( "merge_min_momentum",
                              this_species->merge_min_momentum_log_scale_,
                              "Species", ispec );
            if (this_species->merge_min_momentum_log_scale_ <= 0) {
                ERROR( "In Species " << this_species->name_
                       << ": merge_min_momentum should be above 0.");
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
                ERROR( "For species '" << species_name << "' empty position_initialization" );
            } else if( ( this_species->position_initialization_!="regular" )
                       &&( this_species->position_initialization_!="random" )
                       &&( this_species->position_initialization_!="centered" ) ) {
                this_species->position_initialization_on_species_=true;
            }
        }
        

#ifdef SMILEI_USE_NUMPY
        else if( PyArray_Check( py_pos_init ) ) {
            //Initialize position from this array
            
            PyArrayObject *np_ret = reinterpret_cast<PyArrayObject *>( py_pos_init );
            //Check dimensions
            unsigned int ndim_local = PyArray_NDIM( np_ret ); //Ok
            if( ndim_local != 2 ) ERROR( "For species '" << species_name << "' Provide a 2-dimensional array in order to init particle position from a numpy array." )
                
                //Check number of coordinates provided
                ndim_local = PyArray_SHAPE( np_ret )[0]; // ok
            if( ndim_local != params.nDim_particle + 1 )
                ERROR( "For species '" << species_name << "' position_initializtion must provide a 2-dimensional array with " <<  params.nDim_particle + 1 << " columns." )
                
                // OLD //Get number of particles
                // OLD this_species->n_numpy_particles_ =  PyArray_SHAPE(np_ret)[1];//  ok
                
                //Get number of particles. Do not initialize any more if this is a restart.
                if( !params.restart ) {
                    this_species->n_numpy_particles_ =  PyArray_SHAPE( np_ret )[1];    //  ok
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
            ERROR( "For species '" << species_name << "' non valid position_initialization. It must be either a string or a numpy array." );
        }
        Py_DECREF( py_pos_init );

        if( PyTools::extractV( "regular_number", this_species->regular_number_array_, "Species", ispec )){
             if (this_species->position_initialization_ != "regular") {
                 ERROR("regular_number may not be provided if species position_initialization is not set to 'regular'.");
             }
             if (this_species->regular_number_array_.size() != this_species->nDim_particle) {
                 ERROR("Please provide as many regular numbers of particles as there are particle dimensions in the domain ("<< this_species->nDim_particle <<").");
             }
        }

        PyTools::extract( "ponderomotive_dynamics", this_species->ponderomotive_dynamics, "Species", ispec );
      
        if( this_species->ponderomotive_dynamics && ! params.Laser_Envelope_model ) {
            MESSAGE( "No Laser Envelope is specified - Standard PIC dynamics will be used for all species" );
            this_species->ponderomotive_dynamics = false;
        }


        PyObject *py_mom_init = PyTools::extract_py( "momentum_initialization", "Species", ispec );
        if( PyTools::py2scalar( py_mom_init, this_species->momentum_initialization_ ) ) {
            if( ( this_species->momentum_initialization_=="mj" ) || ( this_species->momentum_initialization_=="maxj" ) ) {
                this_species->momentum_initialization_="maxwell-juettner";
            }
            // Matter particles
            if( this_species->mass_ > 0 ) {
                if( ( this_species->momentum_initialization_!="cold" )
                        && ( this_species->momentum_initialization_!="maxwell-juettner" )
                        && ( this_species->momentum_initialization_!="rectangular" ) ) {
                    ERROR( "For particle species '" << species_name
                           << "' unknown momentum_initialization: "
                           <<this_species->momentum_initialization_ );
                }
            }
            // Photons
            else if( this_species->mass_ == 0 ) {
                if( ( this_species->momentum_initialization_!="cold" )
                        && ( this_species->momentum_initialization_!="rectangular" ) ) {
                    ERROR( "For photon species '" << species_name
                           << "' unknown momentum_initialization: "
                           <<this_species->momentum_initialization_ );
                }
            }
        }
#ifdef SMILEI_USE_NUMPY
        else if( PyArray_Check( py_mom_init ) ) {

            if( !this_species->position_initialization_array_ ) {
                ERROR( "For species '" << species_name << "'. Momentum initialization by a numpy array is only possible if positions are initialized with a numpy array as well. " );
            }

            PyArrayObject *np_ret_mom = reinterpret_cast<PyArrayObject *>( py_mom_init );
            //Check dimensions
            unsigned int ndim_local = PyArray_NDIM( np_ret_mom ) ; //Ok
            if( ndim_local != 2 ) ERROR( "For species '" << species_name << "' Provide a 2-dimensional array in order to init particle momentum from a numpy array." )

                //Check number of coordinates provided
                ndim_local =  PyArray_SHAPE( np_ret_mom )[0]; // ok
            if( ndim_local != 3 )
                ERROR( "For species '" << species_name << "' momentum_initialization must provide a 2-dimensional array with " <<  3 << " columns." )

                //Get number of particles
                if( !params.restart && this_species->n_numpy_particles_ != PyArray_SHAPE( np_ret_mom )[1] )
                    ERROR( "For species '" << species_name << "' momentum_initialization must provide as many particles as position_initialization." )

                    this_species->momentum_initialization_array_ = new double[ndim_local*this_species->n_numpy_particles_] ;
            for( unsigned int idim = 0; idim < ndim_local ; idim++ ) {
                for( unsigned int ipart = 0; ipart < ( unsigned int )this_species->n_numpy_particles_; ipart++ ) {
                    this_species->momentum_initialization_array_[idim*this_species->n_numpy_particles_+ipart] = *( ( double * )PyArray_GETPTR2( np_ret_mom, idim, ipart ) );
                }
            }
        }
#endif
        else {
            ERROR( "For species '" << species_name << "' non valid momentum_initialization. It must be either a string or a numpy array." );
        }
        Py_DECREF( py_mom_init );

        PyTools::extract( "c_part_max", this_species->c_part_max_, "Species", ispec );

        PyTools::extract( "time_frozen", this_species->time_frozen_, "Species", ispec );
        if( this_species->time_frozen_ > 0 && this_species->momentum_initialization_!="cold" ) {
            if( patch->isMaster() ) {
                WARNING( "For species '" << species_name << "' possible conflict between time-frozen & not cold initialization" );
            }
        }
        // time when the relativistic field initialization is applied, if enabled
        int n_timesteps_relativistic_initialization   = ( int )( this_species->time_frozen_/params.timestep );
        this_species->time_relativistic_initialization_ = ( double )( n_timesteps_relativistic_initialization ) * params.timestep;

        if( !PyTools::extractVV( "boundary_conditions", this_species->boundary_conditions, "Species", ispec ) ) {
            ERROR( "For species '" << species_name << "', boundary_conditions not defined" );
        }

        unsigned int number_of_boundaries = (params.geometry=="AMcylindrical") ? 2 : params.nDim_particle;

        if( this_species->boundary_conditions.size() == 0 ) {
            ERROR( "For species '" << species_name << "', boundary_conditions cannot be empty" );
        } else if( this_species->boundary_conditions.size() == 1 ) {
            while( this_species->boundary_conditions.size() < number_of_boundaries ) {
                this_species->boundary_conditions.push_back( this_species->boundary_conditions[0] );
            }
        } else if( this_species->boundary_conditions.size() != number_of_boundaries ) {
            ERROR( "For species '" << species_name << "', boundary_conditions must be of size "<< number_of_boundaries <<"." );
        }


        bool has_thermalize = false;

        for( unsigned int iDim=0; iDim<number_of_boundaries; iDim++ ) {
            if( this_species->boundary_conditions[iDim].size() == 1 ) {
                this_species->boundary_conditions[iDim].push_back( this_species->boundary_conditions[iDim][0] );
            }
            if( this_species->boundary_conditions[iDim].size() != 2 )
                ERROR( "For species '" << species_name << "', boundary_conditions["<<iDim<<"] must have one or two arguments" )
            if( this_species->boundary_conditions[iDim][0] == "thermalize"
                    || this_species->boundary_conditions[iDim][1] == "thermalize" ) {
                has_thermalize = true;
                if( this_species->mass_ == 0 ) {
                    ERROR( "For photon species '" << species_name << "' Thermalizing BCs are not available." );
                }
            }
            if( this_species->boundary_conditions[iDim][0] == "stop"
                    || this_species->boundary_conditions[iDim][1] == "stop" ) {
                if( this_species->mass_ == 0 ) {
                    ERROR( "For photon species '" << species_name << "' stop BCs are not physical." );
                }
            }
        }
        if( (params.geometry=="AMcylindrical") && ( this_species->boundary_conditions[1][1] != "remove" ) && ( this_species->boundary_conditions[1][1] != "stop" ) ) {
            ERROR( " In AM geometry particle boundary conditions supported in Rmax are 'remove' and 'stop' " );
        }
        if( (params.hasWindow) && (( this_species->boundary_conditions[0][1] != "remove" ) || ( this_species->boundary_conditions[0][0] != "remove" ) )) {
            ERROR( " When MovingWindow is activated 'remove' boundary conditions along x is mandatory for all species. " );
        }

        // for thermalizing BCs on particles check if thermal_boundary_temperature is correctly defined
        bool has_temperature = PyTools::extractV( "thermal_boundary_temperature", this_species->thermal_boundary_temperature_, "Species", ispec ) > 0;
        bool has_velocity    = PyTools::extractV( "thermal_boundary_velocity", this_species->thermal_boundary_velocity_, "Species", ispec ) > 0;
        if( has_thermalize ) {
            if( !has_temperature ) {
                ERROR( "For species '" << species_name << "' thermal_boundary_temperature (thermalizing BC) should be a list of floats" );
            }
            if( !has_velocity ) {
                ERROR( "For species '" << species_name << "' thermal_boundary_velocity (thermalizing BC) should be a list of floats" );
            }
            if( this_species->thermal_boundary_velocity_.size()!=3 ) {
                ERROR( "For species '" << species_name << "' thermal_boundary_velocity (thermalizing BC) should have 3 components" );
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
                    ERROR( "For species '" << species_name << "' Thermalizing BCs require non-relativistic thermal_boundary_temperature" );
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
                    ERROR( "For species '" << species_name << ": cannot ionize test species" );
                }
                
                if( ( this_species->atomic_number_==0 )&&( this_species->maximum_charge_state_==0 ) ) {
                    ERROR( "For species '" << species_name << ": undefined atomic_number & maximum_charge_state (required for ionization)" );
                }

                if( model == "tunnel" ){
                    if (params.Laser_Envelope_model){
                        ERROR("An envelope is present, so tunnel_envelope or tunnel_envelope_averaged ionization model should be selected for species "<<species_name);
                    }
                }else if ( (model == "tunnel_envelope") or (model == "tunnel_envelope_averaged") ){

                    if (!params.Laser_Envelope_model) ERROR("An envelope ionization model has been selected but no envelope is present");

                    if (!this_species->ponderomotive_dynamics) ERROR("An envelope ionization model has been selected, but species" <<species_name<<" has ponderomotive_dynamics = False ");
                }else if( model == "from_rate" ) {
                    
                    if( this_species->maximum_charge_state_ == 0 ) {
                        this_species->maximum_charge_state_ = this_species->atomic_number_;
                        WARNING( "For species '" << species_name << ": ionization 'from_rate' is used with maximum_charge_state = "<<this_species->maximum_charge_state_ << " taken from atomic_number" );
                    }
                    this_species->ionization_rate_ = PyTools::extract_py( "ionization_rate", "Species", ispec );
                    if( this_species->ionization_rate_ == Py_None ) {
                        ERROR( "For species '" << species_name << " ionization 'from_rate' requires 'ionization_rate' " );
                    } else {
#ifdef SMILEI_USE_NUMPY
                        PyTools::setIteration( 0 );
                        // Test the ionization_rate function with temporary, "fake" particles
                        std::ostringstream name( "" );
                        name << " ionization_rate:";
                        double *dummy = NULL;
                        ParticleData test( params.nDim_particle, this_species->ionization_rate_, name.str(), dummy );
#else
                        ERROR( "For species '" << species_name << " ionization 'from_rate' requires Numpy" );
#endif
                    }

                } else if( model != "none" ) {
                    ERROR( "For species " << species_name << ": unknown ionization model `" << model );
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
        if( this_species->position_initialization_array_ == NULL ) {
            //These quantities are disregarded if positioning of the species is directly specified by the user
            // Matter particles
            ok1 = PyTools::extract_pyProfile( "number_density", profile1, "Species", ispec );
            ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "Species", ispec );
            if( this_species->mass_ > 0 ) {
                if( ok1 &&  ok2 ) {
                    ERROR( "For species '" << species_name << "', cannot define both `number_density ` and `charge_density`." );
                }
                if( !ok1 && !ok2 ) {
                    ERROR( "For species '" << species_name << "', must define `number_density ` or `charge_density`." );
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
                    ERROR( "For photon species '" << species_name << "', `charge_density` has no meaning."
                            << "You must use `number_density`." );
                }
                if( !ok1 ) {
                    ERROR( "For photon species '" << species_name << "', must define `number_density`." );
                }
                this_species->density_profile_type_ = "nb";
            }

            this_species->density_profile_ = new Profile( profile1, params.nDim_field, Tools::merge( this_species->density_profile_type_, "_density ", species_name ), true );
            MESSAGE(2, "> Density profile: " << this_species->density_profile_->getInfo());
            
            // Number of particles per cell
            if( !PyTools::extract_pyProfile( "particles_per_cell", profile1, "Species", ispec ) ) {
                ERROR( "For species '" << species_name << "', particles_per_cell not found or not understood" );
            }
            this_species->particles_per_cell_profile_ = new Profile( profile1, params.nDim_field, Tools::merge( "particles_per_cell ", species_name ), true );
        } else {
            if( PyTools::extract_pyProfile( "particles_per_cell", profile1, "Species", ispec ) ) {
                ERROR( "For species '" << species_name << "', cannot define both `particles_per_cell` and  `position_initialization` array." );
            }
            ok1 = PyTools::extract_pyProfile( "number_density", profile1, "Species", ispec );
            ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "Species", ispec );
            if( ok1 ||  ok2 ) {
                ERROR( "For species '" << species_name << "', cannot define both `density` and `position_initialization` array." );
            }
        }

        // Charge
        if( !PyTools::extract_pyProfile( "charge", profile1, "Species", ispec ) ) {
            ERROR( "For species '" << species_name << "', charge not found or not understood" );
        }
        this_species->charge_profile_ = new Profile( profile1, params.nDim_field, Tools::merge( "charge ", species_name ), true );
        
        std::vector<PyObject *> prof;
        if( this_species->momentum_initialization_array_ == NULL ) {
            // Mean velocity
            if( PyTools::extract_1or3Profiles( "mean_velocity", "Species", ispec, prof ) ) {
                this_species->velocity_profile_[0] = new Profile( prof[0], params.nDim_field, Tools::merge( "mean_velocity[0] ", species_name ), true );
                this_species->velocity_profile_[1] = new Profile( prof[1], params.nDim_field, Tools::merge( "mean_velocity[1] ", species_name ), true );
                this_species->velocity_profile_[2] = new Profile( prof[2], params.nDim_field, Tools::merge( "mean_velocity[2] ", species_name ), true );
            }
            // Temperature
            if( PyTools::extract_1or3Profiles( "temperature", "Species", ispec, prof ) ) {
                this_species->temperature_profile_[0] = new Profile( prof[0], params.nDim_field, Tools::merge( "temperature[0] ", species_name ), true );
                this_species->temperature_profile_[1] = new Profile( prof[1], params.nDim_field, Tools::merge( "temperature[1] ", species_name ), true );
                this_species->temperature_profile_[2] = new Profile( prof[2], params.nDim_field, Tools::merge( "temperature[2] ", species_name ), true );
            }
        } else {
            ok1 = PyTools::extract_1or3Profiles( "mean_velocity", "Species", ispec, prof ) ;
            ok2 = PyTools::extract_1or3Profiles( "temperature", "Species", ispec, prof ) ;
            if( ok1 ) {
                ERROR( "For species '" << species_name << "', cannot define both `mean_velocity` and `momentum_initialization` array." );
            }
            if( ok2 ) {
                ERROR( "For species '" << species_name << "', cannot define both `temperature` and `momentum_initialization` array." );
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
                    ERROR( "In this version, species '" << species_name << "' cannot be tracked by two DiagTrackParticles" );
                }
                this_species->particles->tracked  = true;
            }
        }

        // Extract test Species flag
        PyTools::extract( "is_test", this_species->particles->is_test, "Species", ispec );

        // Verify they don't ionize
        if( this_species->ionization_model!="none" && this_species->particles->is_test ) {
            ERROR( "For species '" << species_name << "' test & ionized is currently impossible" );
        }

        // Create the particles
        if( !params.restart ) {
            // does a loop over all cells in the simulation
            // considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]
            // Particle creator object
            ParticleCreator particle_creator;
            particle_creator.associate(this_species);
            particle_creator.create( params.n_space, params, patch, 0, 0 );
            
            // this_species->ParticleCreator(params.n_space, params, patch, 0 );

            //MESSAGE(" PARTICLES");
        } else {
            this_species->particles->initialize( 0, params.nDim_particle );
        }

        this_species->initOperators( params, patch );
        //MESSAGE("init operators");
        return this_species;
    } // End Species* create()


    // Method to clone a species from an existing one
    // Note that this must be only called from cloneVector, because additional init is needed
    static Species *clone( Species *species, Params &params, Patch *patch, bool with_particles = true )
    {

        // Create new species object
        Species *new_species = NULL;

        // Boris, Vay or Higuera-Cary
        if ( ( params.vectorization_mode == "off" ) && !params.cell_sorting ) {
            new_species = new SpeciesNorm( params, patch );
        }
#ifdef _VECTO
        else if( ( params.vectorization_mode == "on" ) || params.cell_sorting  ) {
            new_species = new SpeciesNormV( params, patch );
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
        new_species->photon_species_                            = species->photon_species_;
        new_species->species_number_                           = species->species_number_;
        new_species->position_initialization_on_species_       = species->position_initialization_on_species_;
        new_species->position_initialization_on_species_type_  = species->position_initialization_on_species_type_;
        new_species->position_initialization_on_species_index  = species->position_initialization_on_species_index;
        new_species->position_initialization_                  = species->position_initialization_;
        new_species->position_initialization_array_            = species->position_initialization_array_;
        new_species->regular_number_array_                     = species->regular_number_array_;
        new_species->n_numpy_particles_                        = species->n_numpy_particles_            ;
        new_species->momentum_initialization_                  = species->momentum_initialization_;
        new_species->momentum_initialization_array_            = species->momentum_initialization_array_;
        new_species->c_part_max_                               = species->c_part_max_;
        new_species->mass_                                     = species->mass_;
        new_species->time_frozen_                              = species->time_frozen_;
        new_species->radiating_                                = species->radiating_;
        new_species->relativistic_field_initialization_        = species->relativistic_field_initialization_;
        new_species->time_relativistic_initialization_         = species->time_relativistic_initialization_;
        new_species->boundary_conditions                      = species->boundary_conditions;
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
        new_species->ponderomotive_dynamics                   = species->ponderomotive_dynamics;

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


        // \todo : NOT SURE HOW THIS BEHAVES WITH RESTART
        if( ( !params.restart ) && ( with_particles ) ) {
            ParticleCreator particle_creator;
            particle_creator.associate(new_species);
            particle_creator.create( params.n_space, params, patch, 0, 0 );
            
            // new_species->ParticleCreator( params.n_space, params, patch, 0 );

        } else {
            new_species->particles->initialize( 0, ( *species->particles ) );
        }

        new_species->initOperators( params, patch );

        return new_species;
    } // End Species* clone()


    static std::vector<Species *> createVector( Params &params, Patch *patch )
    {
        // this will be returned
        std::vector<Species *> returned_species;
        returned_species.resize( 0 );

        // read from python namelist
        unsigned int tot_species_number = PyTools::nComponents( "Species" );
        if (tot_species_number > 0) {
            TITLE("Initializing species");
        }
        for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {
            Species *this_species = SpeciesFactory::create( params, ispec, patch );
            // Verify the new species does not have the same name as a previous one
            for( unsigned int i = 0; i < ispec; i++ ) {
                if( this_species->name_ == returned_species[i]->name_ ) {
                    ERROR("Two species cannot have the same name `"<<this_species->name_<<"`");
                }
            }
            // Put the newly created species in the vector of species
            returned_species.push_back( this_species );
        }

        // Loop species to find species which their particles positions is on another species
        for( unsigned int ispec1 = 0; ispec1<returned_species.size(); ispec1++ ) {
            if( returned_species[ispec1]->position_initialization_on_species_==true ) {
                // If true then position_initialization of spec1 is not 'centered', 'regular' or 'random'
                // So we have to check if :
                // - 'position_initialization' of spec1 is another already created specie name;
                // - 'position_initialization' of spec1 is not the spec1 name;
                // - 'position_initialization' of spec2 is centered,regular,random;
                // - The number of particle of spec1 is equal to spec2

                // Loop all other species
                for( unsigned int ispec2 = 0; ispec2<returned_species.size(); ispec2++ ) {
                    if( returned_species[ispec1]->position_initialization_ == returned_species[ispec2]->name_ ) {
                        if( returned_species[ispec1]->position_initialization_==returned_species[ispec1]->name_ ) {
                            ERROR( "For species '"<<returned_species[ispec1]->name_<<"' position_initialization must be different from '"<<returned_species[ispec1]->name_<<"'." );
                        }
                        if( returned_species[ispec2]->position_initialization_on_species_==true ) {
                            ERROR( "For species '"<<returned_species[ispec2]->name_<<"' position_initialization must be 'centered', 'regular' or 'random' (pre-defined position) in order to attach '"<<returned_species[ispec1]->name_<<"' to its initial position." );
                        }
                        if( returned_species[ispec1]->getNbrOfParticles() != returned_species[ispec2]->getNbrOfParticles() ) {
                            ERROR( "Number of particles in species '"<<returned_species[ispec1]->name_<<"' is not equal to the number of particles in species '"<<returned_species[ispec2]->name_<<"'." );
                        }
                        // We copy ispec2 which is the index of the species, already created, on which initialize particle of the new created species
                        returned_species[ispec1]->position_initialization_on_species_index=ispec2;
                        returned_species[ispec1]->position_initialization_on_species_type_ = returned_species[ispec2]->position_initialization_;
                        // We copy position of species 2 (index ispec2), for position on species 1 (index ispec1)
                        returned_species[ispec1]->particles->Position=returned_species[ispec2]->particles->Position;
                    }
                }
                if( returned_species[ispec1]->position_initialization_on_species_index==-1 ) {
                    ERROR( "Specie '"<<returned_species[ispec1]->position_initialization_<<"' doesn't exist. We can't initialize position on this species. Choose an already created specie or 'centered', 'regular', 'random'." )
                }
            } else {
                returned_species[ispec1]->position_initialization_on_species_type_ = returned_species[ispec1]->position_initialization_;
            }
        }

        // Update particles weight in specific case
        if (params.geometry=="AMcylindrical") {
            for( unsigned int ispec1 = 0; ispec1<returned_species.size(); ispec1++ ) {
                ParticleCreator::regulateWeightwithPositionAM( returned_species[ispec1]->particles, returned_species[ispec1]->position_initialization_on_species_type_, returned_species[ispec1]->cell_length[1]  );
            }
        }

        // Loop species to find related species
        for( unsigned int ispec1 = 0; ispec1<returned_species.size(); ispec1++ ) {
            
            // Ionizable species
            if( returned_species[ispec1]->Ionize ) {
                // Loop all other species
                for( unsigned int ispec2 = 0; ispec2<returned_species.size(); ispec2++ ) {
                    if( returned_species[ispec1]->ionization_electrons == returned_species[ispec2]->name_ ) {
                        if( ispec1==ispec2 ) {
                            ERROR( "For species '"<<returned_species[ispec1]->name_<<"' ionization_electrons must be a distinct species" );
                        }
                        if( returned_species[ispec2]->mass_!=1 ) {
                            ERROR( "For species '"<<returned_species[ispec1]->name_<<"' ionization_electrons must be a species with mass==1" );
                        }
                        returned_species[ispec1]->electron_species_index = ispec2;
                        returned_species[ispec1]->electron_species = returned_species[ispec2];
                        
                        int max_eon_number =
                            returned_species[ispec1]->getNbrOfParticles()
                            * ( returned_species[ispec1]->atomic_number_ || returned_species[ispec1]->maximum_charge_state_ );
                        returned_species[ispec1]->Ionize->new_electrons.initializeReserve(
                            max_eon_number, *returned_species[ispec1]->electron_species->particles
                        );
                        break;
                    }
                }
                if( returned_species[ispec1]->electron_species_index==-1 ) {
                    ERROR( "For species '"<<returned_species[ispec1]->name_<<"' ionization_electrons named " << returned_species[ispec1]->ionization_electrons << " could not be found" );
                }
            }
            
            // Radiating species
            if( returned_species[ispec1]->Radiate ) {
                // No emission of discrete photon, only scalar diagnostics are updated
                if( returned_species[ispec1]->radiation_photon_species.empty() ) {
                    returned_species[ispec1]->photon_species_index = -1;
                    returned_species[ispec1]->photon_species_ = NULL;
                }
                // Else, there will be emission of macro-photons.
                else {
                    unsigned int ispec2 = 0;
                    for( ispec2 = 0; ispec2<returned_species.size(); ispec2++ ) {
                        if( returned_species[ispec1]->radiation_photon_species == returned_species[ispec2]->name_ ) {
                            if( ispec1==ispec2 ) {
                                ERROR( "For species '"<<returned_species[ispec1]->name_<<"' radiation_photon_species must be a distinct photon species" );
                            }
                            if( returned_species[ispec2]->mass_!=0 ) {
                                ERROR( "For species '"<<returned_species[ispec1]->name_<<"' radiation_photon_species must be a photon species with mass==0" );
                            }
                            returned_species[ispec1]->photon_species_index = ispec2;
                            returned_species[ispec1]->photon_species_ = returned_species[ispec2];
                            returned_species[ispec1]->Radiate->new_photons_.initializeReserve(
                                returned_species[ispec1]->getNbrOfParticles(),
                                *returned_species[ispec1]->photon_species_->particles
                            );
                            break;
                        }
                    }
                    if( ispec2 == returned_species.size() ) {
                        ERROR( "Species '" << returned_species[ispec1]->radiation_photon_species << "' does not exist." )
                    }
                }
            }
            
            // Breit-Wheeler species
            if( returned_species[ispec1]->Multiphoton_Breit_Wheeler_process ) {
                unsigned int ispec2;
                for( int k=0; k<2; k++ ) {
                    ispec2 = 0;
                    while( ispec2<returned_species.size()) {
                        // We look for the pair species multiphoton_Breit_Wheeler_[k]
                        if( returned_species[ispec1]->multiphoton_Breit_Wheeler_[k] == returned_species[ispec2]->name_ ) {
                            if( ispec1==ispec2 ) {
                                ERROR( "For species '" << returned_species[ispec1]->name_
                                       << "' pair species must be a distinct particle species" );
                            }
                            if( returned_species[ispec2]->mass_ != 1 ) {
                                ERROR( "For species '"<<returned_species[ispec1]->name_
                                  <<"' pair species must be an electron and positron species (mass = 1)" );
                            }
                            returned_species[ispec1]->mBW_pair_species_index[k] = ispec2;
                            returned_species[ispec1]->mBW_pair_species[k] = returned_species[ispec2];
                            returned_species[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].initializeReserve(
                                returned_species[ispec1]->getNbrOfParticles(),
                                *returned_species[ispec1]->mBW_pair_species[k]->particles
                            );
                            ispec2 = returned_species.size() + 1;
                        }
                        ispec2++ ;
                    }
                    // This means that one of the pair species has not been fould
                    if( ispec2 == returned_species.size() ) {
                        ERROR( "In Species `" << returned_species[ispec1]->name_ << "`"
                           << " the pair species `" << returned_species[ispec1]->multiphoton_Breit_Wheeler_[k]
                           << "` does not exist." )
                    }
                }
            }
        }

        return returned_species;
    }

    // -----------------------------------------------------------------------------------------------------------------
    //! Method to clone the whole vector of species
    // -----------------------------------------------------------------------------------------------------------------
    static std::vector<Species *> cloneVector( std::vector<Species *> vector_species, Params &params, Patch *patch, bool with_particles = true )
    {
        std::vector<Species *> returned_species;
        returned_species.resize( 0 );

        for( unsigned int ispec = 0; ispec < vector_species.size(); ispec++ ) {
            Species *new_species = SpeciesFactory::clone( vector_species[ispec], params, patch, with_particles );
            returned_species.push_back( new_species );
        }
        patch->copyPositions(returned_species);
        
        // Update particles weight in specific case
        if (params.geometry=="AMcylindrical") {
            for( unsigned int ispec1 = 0; ispec1<returned_species.size(); ispec1++ ) {
                ParticleCreator::regulateWeightwithPositionAM( returned_species[ispec1]->particles, returned_species[ispec1]->position_initialization_on_species_type_, returned_species[ispec1]->cell_length[1]);
            }
        }

        // Ionization
        for( unsigned int i=0; i<returned_species.size(); i++ ) {
            if( returned_species[i]->Ionize ) {
                returned_species[i]->electron_species_index = vector_species[i]->electron_species_index;
                returned_species[i]->electron_species = returned_species[returned_species[i]->electron_species_index];
                returned_species[i]->Ionize->new_electrons.tracked = returned_species[i]->electron_species->particles->tracked;
                returned_species[i]->Ionize->new_electrons.isQuantumParameter = returned_species[i]->electron_species->particles->isQuantumParameter;
                returned_species[i]->Ionize->new_electrons.isMonteCarlo = returned_species[i]->electron_species->particles->isMonteCarlo;
                returned_species[i]->Ionize->new_electrons.initialize( 0, params.nDim_particle );
            }
        }

        // Synchrotron-like radiation
        for( unsigned int i=0; i<returned_species.size(); i++ ) {
            if( returned_species[i]->Radiate ) {
                returned_species[i]->radiation_photon_species = vector_species[i]->radiation_photon_species;
                returned_species[i]->photon_species_index = vector_species[i]->photon_species_index;
                if( vector_species[i]->photon_species_ ) {
                    returned_species[i]->photon_species_ = returned_species[returned_species[i]->photon_species_index];
                    returned_species[i]->Radiate->new_photons_.tracked = returned_species[i]->photon_species_->particles->tracked;
                    returned_species[i]->Radiate->new_photons_.isQuantumParameter = returned_species[i]->photon_species_->particles->isQuantumParameter;
                    returned_species[i]->Radiate->new_photons_.isMonteCarlo = returned_species[i]->photon_species_->particles->isMonteCarlo;
                    //returned_species[i]->Radiate->new_photons_.initialize(returned_species[i]->getNbrOfParticles(),
                    //                                               params.nDim_particle );
                    returned_species[i]->Radiate->new_photons_.initialize( 0, params.nDim_particle );
                } else {
                    returned_species[i]->photon_species_ = NULL;
                }
            }
        }

        // multiphoton Breit-Wheeler
        for( unsigned int i=0; i<returned_species.size(); i++ ) {
            if( returned_species[i]->Multiphoton_Breit_Wheeler_process ) {
                // Loop on pairs
                for( int k=0; k<2; k++ ) {
                    returned_species[i]->multiphoton_Breit_Wheeler_[k] = vector_species[i]->multiphoton_Breit_Wheeler_[k];
                    returned_species[i]->mBW_pair_species_index[k] = vector_species[i]->mBW_pair_species_index[k];
                    returned_species[i]->mBW_pair_species[k] = returned_species[returned_species[i]->mBW_pair_species_index[k]];
                    returned_species[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].tracked = returned_species[i]->mBW_pair_species[k]->particles->tracked;
                    returned_species[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].isQuantumParameter = returned_species[i]->mBW_pair_species[k]->particles->isQuantumParameter;
                    returned_species[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].isMonteCarlo = returned_species[i]->mBW_pair_species[k]->particles->isMonteCarlo;
                    returned_species[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].initialize(
                        0, params.nDim_particle );
                }
            } else {
                returned_species[i]->mBW_pair_species[0] = NULL;
                returned_species[i]->mBW_pair_species[1] = NULL;
            }
        }

        return returned_species;
    }

};

#endif
