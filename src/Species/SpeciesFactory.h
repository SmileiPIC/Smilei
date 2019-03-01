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
#include "SpeciesVAdaptive.h"
#include "SpeciesVAdaptiveMixedSort.h"
#include "SpeciesNormV.h"
#endif

#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "PartBoundCond.h"

#include "Params.h"
#include "Patch.h"

#include "ParticleData.h"

#include "Tools.h"
#ifdef SMILEI_USE_NUMPY
#include <numpy/arrayobject.h>
#endif

class SpeciesFactory
{
public:
    static Species *create( Params &params, int ispec, Patch *patch )
    {
    
        std::string species_name( "" );
        
        PyTools::extract( "name", species_name, "Species", ispec );
        if( patch->isMaster() ) {
            MESSAGE( 1, "Creating Species : " << species_name );
        }
        
        unsigned int tot_species_number = PyTools::nComponents( "Species" );
        if( species_name.empty() ) {
            std::ostringstream name( "" );
            name << "species" << std::setfill( '0' ) << std::setw( log10( tot_species_number )+1 ) << ispec;
            species_name=name.str();
            if( patch->isMaster() ) {
                MESSAGE( "For species #" << ispec << ", name will be " << species_name );
            }
        }
        
        // Extract type of species dynamics from namelist
        std::string pusher = "boris"; // default value
        PyTools::extract( "pusher", pusher, "Species", ispec );
        
        // Extract type of species radiation from namelist
        std::string radiation_model = "none"; // default value
        if( !PyTools::extract( "radiation_model", radiation_model, "Species", ispec ) ) {
            if( patch->isMaster() ) {
                WARNING( "For species '" << species_name << "', radiation_model not defined: assumed = 'none'." );
            }
        } else {
            // Cancelation of the letter case for `radiation_model`
            std::transform( radiation_model.begin(), radiation_model.end(), radiation_model.begin(), tolower );
            
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
        }
        
        // Extract mass from namelist
        double mass;
        if( !PyTools::extract( "mass", mass, "Species", ispec ) ) {
            if( patch->isMaster() ) {
                ERROR( "For species '" << species_name << "`, mass is not defined." );
            }
        }
        
        // Create species object
        Species *thisSpecies = NULL;
        
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
                if( params.vectorization_mode == "off" ) {
                    thisSpecies = new SpeciesNorm( params, patch );
                }
                
#ifdef _VECTO
                else if( params.vectorization_mode == "on" ) {
                    thisSpecies = new SpeciesNormV( params, patch );
                } else if( params.vectorization_mode == "adaptive_mixed_sort" ) {
                    thisSpecies = new SpeciesVAdaptive( params, patch );
                } else if( params.vectorization_mode == "adaptive" ) {
                    thisSpecies = new SpeciesVAdaptiveMixedSort( params, patch );
                }
#endif
            } else {
                ERROR( "For species `" << species_name << "`, pusher must be 'boris', 'borisnr', 'vay', 'higueracary', 'ponderomotive_boris'" );
            }
            thisSpecies->pusher = pusher;
            
            // Radiation model of the species
            // Species with a Monte-Carlo process for the radiation loss
            if( radiation_model=="mc" ) {
                thisSpecies->particles->isQuantumParameter = true;
                thisSpecies->particles->isMonteCarlo = true;
                thisSpecies->radiating = true;
            }
            // Species with another radiation loss model
            else if( ( radiation_model=="ll" )
                     || ( radiation_model=="cll" )
                     || ( radiation_model=="niel" ) ) {
                thisSpecies->particles->isQuantumParameter = true;
                thisSpecies->radiating = true;
            } else if( radiation_model != "none" ) {
                ERROR( "For species `" << species_name
                       << " radiation_model must be 'none',"
                       << " 'Landau-Lifshitz' ('ll'),"
                       << " 'corrected-Landau-Lifshitz' ('cll'),"
                       << " 'Niel' ('niel') or 'Monte-Carlo' ('mc')" );
            }
            
            thisSpecies->radiation_model = radiation_model;
            
            if( radiation_model == "ll" ) {
                MESSAGE( 2, "> Radiating species with the classical Landau-Lifshitz radiating model" );
            } else if( radiation_model == "cll" ) {
                MESSAGE( 2, "> Radiating species with the quantum corrected Landau-Lifshitz radiating model" );
            } else if( radiation_model == "niel" ) {
                MESSAGE( 2, "> Radiating species with the stochastic model of Niel et al." );
            } else if( radiation_model == "mc" ) {
                MESSAGE( 2, "> Radiating species with the stochastic Monte-Carlo model" );
            } else if( radiation_model != "none" ) {
                MESSAGE( 2, "> Radiating species with model: `" << radiation_model << "`" );
            }
            
            // Non compatibility
            if( ( pusher=="borisnr" )
                    && ( radiation_model=="mc"
                         || radiation_model=="ll"
                         || radiation_model=="cll"
                         || radiation_model=="niel" ) ) {
                ERROR( "For species `" << species_name
                       << "` radiation_model `"
                       << radiation_model
                       << "` is not compatible with pusher "
                       << pusher );
                       
            }
            
        }
        
        // Photon species
        else if( mass == 0 ) {
            if( params.vectorization_mode == "off" ) {
                thisSpecies = new SpeciesNorm( params, patch );
            }
#ifdef _VECTO
            else if( params.vectorization_mode == "on" ) {
                thisSpecies = new SpeciesNormV( params, patch );
            } else if( params.vectorization_mode == "adaptive_mixed_sort" ) {
                thisSpecies = new SpeciesVAdaptive( params, patch );
            } else if( params.vectorization_mode == "adaptive" ) {
                thisSpecies = new SpeciesVAdaptiveMixedSort( params, patch );
            }
#endif
            // Photon can not radiate
            thisSpecies->radiation_model = "none";
            thisSpecies-> pusher = "norm";
            
            MESSAGE( 2, "> " <<species_name <<" is a photon species (mass==0)." );
            MESSAGE( 2, "> Radiation model set to none." );
            MESSAGE( 2, "> Pusher set to norm." );
        }
        
        thisSpecies->name = species_name;
        thisSpecies->mass = mass;
        thisSpecies->speciesNumber = ispec;
        
        // Vectorized operators
        if( params.vectorization_mode == "off" ) {
            thisSpecies->vectorized_operators = false;
        } else if( params.vectorization_mode == "on" || params.vectorization_mode == "adaptive_mixed_sort" || params.vectorization_mode == "adaptive" ) {
            thisSpecies->vectorized_operators = true;
        }
        
        // Extract various parameters from the namelist
        
        // Monte-Carlo Photon emission properties
        if( mass > 0. ) {
            if( thisSpecies->radiation_model == "mc" ) {
                if( PyTools::extract( "radiation_photon_species", thisSpecies->radiation_photon_species, "Species", ispec ) ) {
                
                    MESSAGE( 2, "> Macro-photon emission activated" );
                    
                    // Species that will receive the emitted photons
                    if( !thisSpecies->radiation_photon_species.empty() ) {
                        MESSAGE( 2, "> Emitted photon species set to `" << thisSpecies->radiation_photon_species << "`" );
                    } else {
                        ERROR( " The radiation photon species is not specified." )
                    }
                    
                    // Number of photons emitted per Monte-Carlo event
                    if( PyTools::extract( "radiation_photon_sampling",
                                          thisSpecies->radiation_photon_sampling_, "Species", ispec ) ) {
                        if( thisSpecies->radiation_photon_sampling_ < 1 ) {
                            ERROR( "For species '" << species_name
                                   << "' radiation_photon_sampling should be > 1" );
                        }
                    } else {
                        thisSpecies->radiation_photon_sampling_ = 1;
                    }
                    MESSAGE( 2, "> Number of macro-photons emitted per MC event: "
                             << thisSpecies->radiation_photon_sampling_ );
                             
                    // Photon energy threshold
                    if( !PyTools::extract( "radiation_photon_gamma_threshold",
                                           thisSpecies->radiation_photon_gamma_threshold_, "Species", ispec ) ) {
                        thisSpecies->radiation_photon_gamma_threshold_ = 2.;
                    }
                    MESSAGE( 2, "> Photon energy threshold for macro-photon emission: "
                             << thisSpecies->radiation_photon_gamma_threshold_ );
                } else {
                    MESSAGE( 2, "> Macro-photon emission not activated" );
                }
                
            }
        }
        
        // Multiphoton Breit-Wheeler
        if( mass == 0 ) {
            // If thisSpecies->multiphoton_Breit_Wheeler
            if( PyTools::extract( "multiphoton_Breit_Wheeler", thisSpecies->multiphoton_Breit_Wheeler, "Species", ispec ) ) {
                // If one of the species is empty
                if( thisSpecies->multiphoton_Breit_Wheeler[1].empty() || thisSpecies->multiphoton_Breit_Wheeler[0].empty() ) {
                    ERROR( "For species '" << species_name
                           << "' multiphoton_Breit_Wheeler can not be empty,"
                           << " select electron and positron species." );
                } else {
                    // Activation of the additional variables
                    thisSpecies->particles->isQuantumParameter = true;
                    thisSpecies->particles->isMonteCarlo = true;
                    
                    MESSAGE( 2, "> Decay into pair via the multiphoton Breit-Wheeler activated" );
                    MESSAGE( 2, "> Generated electrons and positrons go to species: "
                             << thisSpecies->multiphoton_Breit_Wheeler[0]
                             << " & " << thisSpecies->multiphoton_Breit_Wheeler[1] );
                             
                    // Number of emitted particles per MC event
                    thisSpecies->mBW_pair_creation_sampling.resize( 2 );
                    if( !PyTools::extract( "multiphoton_Breit_Wheeler_sampling",
                                           thisSpecies->mBW_pair_creation_sampling, "Species", ispec ) ) {
                        thisSpecies->mBW_pair_creation_sampling[0] = 1;
                        thisSpecies->mBW_pair_creation_sampling[1] = 1;
                    }
                    MESSAGE( 2, "> Number of emitted macro-particles per MC event: "
                             << thisSpecies->mBW_pair_creation_sampling[0]
                             << " & " << thisSpecies->mBW_pair_creation_sampling[1] );
                }
            }
        }
        
        // Particle Merging
        
        // Extract merging method
        thisSpecies->merging_method_ = "none"; // default value
        if( PyTools::extract( "merging_method", thisSpecies->merging_method_, "Species", ispec ) ) {
            // Cancelation of the letter case for `merging_method_`
            std::transform( thisSpecies->merging_method_.begin(),
                            thisSpecies->merging_method_.end(),
                            thisSpecies->merging_method_.begin(), tolower );
                            
            if( ( thisSpecies->merging_method_ != "vranic" ) &&
                    ( thisSpecies->merging_method_ != "none" ) ) {
                ERROR( "In Species " << thisSpecies->name
                       << ": merging method not valid, must be `vranic` or `none`" );
            }
            
            if( params.vectorization_mode == "off" && thisSpecies->merging_method_ != "none" ) {
                ERROR( "In Species " << thisSpecies->name
                       << ": particle merging only available with `vectorization_mode` = `on` or `adpative`" );
            }
            
            // get parameter "every" which describes a timestep selection
            if( !thisSpecies->merging_time_selection_ ) {
                thisSpecies->merging_time_selection_ = new TimeSelection(
                    PyTools::extract_py( "merge_every", "Species", ispec ), "Particle merging"
                );
            }
            
            if( thisSpecies->merging_method_ != "none" ) {
                MESSAGE( 2, "> Particle merging with the method: "
                         << thisSpecies->merging_method_ );
                MESSAGE( 2, "> Merging time selection: "
                         << thisSpecies->merging_time_selection_->info() );
            }
        }
        
        PyObject *py_pos_init = PyTools::extract_py( "position_initialization", "Species", ispec );
        if( PyTools::convert( py_pos_init, thisSpecies->position_initialization ) ) {
            if( thisSpecies->position_initialization.empty() ) {
                ERROR( "For species '" << species_name << "' empty position_initialization" );
            } else if( ( thisSpecies->position_initialization!="regular" )
                       &&( thisSpecies->position_initialization!="random" )
                       &&( thisSpecies->position_initialization!="centered" ) ) {
                thisSpecies->position_initialization_on_species=true;
                
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
                // OLD thisSpecies->n_numpy_particles =  PyArray_SHAPE(np_ret)[1];//  ok
                
                //Get number of particles. Do not initialize any more if this is a restart.
                if( !params.restart ) {
                    thisSpecies->n_numpy_particles =  PyArray_SHAPE( np_ret )[1];    //  ok
                }
            thisSpecies->position_initialization_array = new double[ndim_local*thisSpecies->n_numpy_particles] ;
            for( unsigned int idim = 0; idim < ndim_local ; idim++ ) {
                for( unsigned int ipart = 0; ipart < ( unsigned int )thisSpecies->n_numpy_particles; ipart++ ) {
                    thisSpecies->position_initialization_array[idim*thisSpecies->n_numpy_particles+ipart] = *( ( double * )PyArray_GETPTR2( np_ret, idim, ipart ) );
                }
            }
        }
#endif
        else {
            ERROR( "For species '" << species_name << "' non valid position_initialization. It must be either a string or a numpy array." );
        }
        Py_DECREF( py_pos_init );
        
        
        PyTools::extract( "ponderomotive_dynamics", thisSpecies->ponderomotive_dynamics, "Species", ispec );
        if( thisSpecies->ponderomotive_dynamics && ( params.geometry != "1Dcartesian" ) && ( params.geometry != "2Dcartesian" ) && ( params.geometry != "3Dcartesian" ) ) {
            ERROR( "Ponderomotive/Envelope model only available in 1D, 2D, 3D cartesian geometry" );
        }
        int n_envlaser = PyTools::nComponents( "LaserEnvelope" );
        if( thisSpecies->ponderomotive_dynamics && ( n_envlaser < 1 ) ) {
            MESSAGE( "No Laser Envelope is specified - Standard PIC dynamics will be used for all species" );
            thisSpecies->ponderomotive_dynamics = false;
        }
        
        
        PyObject *py_mom_init = PyTools::extract_py( "momentum_initialization", "Species", ispec );
        if( PyTools::convert( py_mom_init, thisSpecies->momentum_initialization ) ) {
            if( ( thisSpecies->momentum_initialization=="mj" ) || ( thisSpecies->momentum_initialization=="maxj" ) ) {
                thisSpecies->momentum_initialization="maxwell-juettner";
            }
            // Matter particles
            if( thisSpecies->mass > 0 ) {
                if( ( thisSpecies->momentum_initialization!="cold" )
                        && ( thisSpecies->momentum_initialization!="maxwell-juettner" )
                        && ( thisSpecies->momentum_initialization!="rectangular" ) ) {
                    ERROR( "For particle species '" << species_name
                           << "' unknown momentum_initialization: "
                           <<thisSpecies->momentum_initialization );
                }
            }
            // Photons
            else if( thisSpecies->mass == 0 ) {
                if( ( thisSpecies->momentum_initialization!="cold" )
                        && ( thisSpecies->momentum_initialization!="rectangular" ) ) {
                    ERROR( "For photon species '" << species_name
                           << "' unknown momentum_initialization: "
                           <<thisSpecies->momentum_initialization );
                }
            }
        }
#ifdef SMILEI_USE_NUMPY
        else if( PyArray_Check( py_mom_init ) ) {
        
            if( !thisSpecies->position_initialization_array ) {
                ERROR( "For species '" << species_name << "'. Momentum initialization by a numpy array is only possible if positions are initialized with a numpy array as well. " );
            }
            
            PyArrayObject *np_ret_mom = reinterpret_cast<PyArrayObject *>( py_mom_init );
            //Check dimensions
            unsigned int ndim_local = PyArray_NDIM( np_ret_mom ) ; //Ok
            if( ndim_local != 2 ) ERROR( "For species '" << species_name << "' Provide a 2-dimensional array in order to init particle momentum from a numpy array." )
            
                //Check number of coordinates provided
                ndim_local =  PyArray_SHAPE( np_ret_mom )[0]; // ok
            if( ndim_local != 3 )
                ERROR( "For species '" << species_name << "' momentum_initializtion must provide a 2-dimensional array with " <<  3 << " columns." )
                
                //Get number of particles
                if( !params.restart && thisSpecies->n_numpy_particles != PyArray_SHAPE( np_ret_mom )[1] )
                    ERROR( "For species '" << species_name << "' momentum_initialization must provide as many particles as position_initialization." )
                    
                    thisSpecies->momentum_initialization_array = new double[ndim_local*thisSpecies->n_numpy_particles] ;
            for( unsigned int idim = 0; idim < ndim_local ; idim++ ) {
                for( unsigned int ipart = 0; ipart < ( unsigned int )thisSpecies->n_numpy_particles; ipart++ ) {
                    thisSpecies->momentum_initialization_array[idim*thisSpecies->n_numpy_particles+ipart] = *( ( double * )PyArray_GETPTR2( np_ret_mom, idim, ipart ) );
                }
            }
        }
#endif
        else {
            ERROR( "For species '" << species_name << "' non valid momentum_initialization. It must be either a string or a numpy array." );
        }
        Py_DECREF( py_mom_init );
        
        PyTools::extract( "c_part_max", thisSpecies->c_part_max, "Species", ispec );
        
        PyTools::extract( "time_frozen", thisSpecies->time_frozen, "Species", ispec );
        if( thisSpecies->time_frozen > 0 && thisSpecies->momentum_initialization!="cold" ) {
            if( patch->isMaster() ) {
                WARNING( "For species '" << species_name << "' possible conflict between time-frozen & not cold initialization" );
            }
        }
        // time when the relativistic field initialization is applied, if enabled
        int n_timesteps_relativistic_initialization   = ( int )( thisSpecies->time_frozen/params.timestep );
        thisSpecies->time_relativistic_initialization = ( double )( n_timesteps_relativistic_initialization ) * params.timestep;
        
        if( !PyTools::extract( "boundary_conditions", thisSpecies->boundary_conditions, "Species", ispec ) ) {
            ERROR( "For species '" << species_name << "', boundary_conditions not defined" );
        }
        if( params.geometry != "AMcylindrical" ) {
            if( thisSpecies->boundary_conditions.size() == 0 ) {
                ERROR( "For species '" << species_name << "', boundary_conditions cannot be empty" );
            } else if( thisSpecies->boundary_conditions.size() == 1 ) {
                while( thisSpecies->boundary_conditions.size() < params.nDim_particle ) {
                    thisSpecies->boundary_conditions.push_back( thisSpecies->boundary_conditions[0] );
                }
            } else if( thisSpecies->boundary_conditions.size() != params.nDim_particle ) {
                ERROR( "For species '" << species_name << "', boundary_conditions must be the same size as the number of dimensions" );
            }
        } else if( params.geometry == "AMcylindrical" ) {
            if( thisSpecies->boundary_conditions.size() == 0 ) {
                ERROR( "For species '" << species_name << "', boundary_conditions cannot be empty" );
            } else if( thisSpecies->boundary_conditions.size() == 1 ) {
                while( thisSpecies->boundary_conditions.size() < params.nDim_particle ) {
                    thisSpecies->boundary_conditions.push_back( thisSpecies->boundary_conditions[0] );
                }
            } else if( thisSpecies->boundary_conditions.size() != 2 ) {
                ERROR( "For AM geometry boundary_conditions must not be the same size as the number of dimensions it is applied only for Rmax Xmin and Xmax" );
            }
            if( ( thisSpecies->boundary_conditions[1][1] != "remove" ) && ( thisSpecies->boundary_conditions[1][1] != "stop" ) ) {
                ERROR( " In AM geometry particle boundary conditions supported in Rmax are 'remove' and 'stop' " );
            }
        }
        bool has_thermalize = false;
        if( params.geometry != "AMcylindrical" ) {
            for( unsigned int iDim=0; iDim<params.nDim_particle; iDim++ ) {
                if( thisSpecies->boundary_conditions[iDim].size() == 1 ) {
                    thisSpecies->boundary_conditions[iDim].push_back( thisSpecies->boundary_conditions[iDim][0] );
                }
                if( thisSpecies->boundary_conditions[iDim].size() != 2 )
                    ERROR( "For species '" << species_name << "', boundary_conditions["<<iDim<<"] must have one or two arguments" )
                    if( thisSpecies->boundary_conditions[iDim][0] == "thermalize"
                            || thisSpecies->boundary_conditions[iDim][1] == "thermalize" ) {
                        has_thermalize = true;
                        if( thisSpecies->mass == 0 ) {
                            ERROR( "For photon species '" << species_name << "' Thermalizing BCs are not available." );
                        }
                    }
                if( thisSpecies->boundary_conditions[iDim][0] == "stop"
                        || thisSpecies->boundary_conditions[iDim][1] == "stop" ) {
                    if( thisSpecies->mass == 0 ) {
                        ERROR( "For photon species '" << species_name << "' stop BCs are not physical." );
                    }
                }
            }
        }
        // for thermalizing BCs on particles check if thermal_boundary_temperature is correctly defined
        bool has_temperature = PyTools::extract( "thermal_boundary_temperature", thisSpecies->thermal_boundary_temperature, "Species", ispec );
        bool has_velocity    = PyTools::extract( "thermal_boundary_velocity", thisSpecies->thermal_boundary_velocity, "Species", ispec );
        if( has_thermalize ) {
            if( !has_temperature ) {
                ERROR( "For species '" << species_name << "' thermal_boundary_temperature (thermalizing BC) should be a list of floats" );
            }
            if( !has_velocity ) {
                ERROR( "For species '" << species_name << "' thermal_boundary_velocity (thermalizing BC) should be a list of floats" );
            }
            if( thisSpecies->thermal_boundary_velocity.size()!=3 ) {
                ERROR( "For species '" << species_name << "' thermal_boundary_velocity (thermalizing BC) should have 3 components" );
            }
            if( thisSpecies->thermal_boundary_temperature.size()==1 ) {
                WARNING( "For species '" << species_name << "' Using thermal_boundary_temperature[0] in all directions" );
                thisSpecies->thermal_boundary_temperature.resize( 3 );
                thisSpecies->thermal_boundary_temperature[1] = thisSpecies->thermal_boundary_temperature[0];
                thisSpecies->thermal_boundary_temperature[2] = thisSpecies->thermal_boundary_temperature[0];
            }
            
            // Compute the thermalVelocity & Momentum for thermalizing bcs
            thisSpecies->thermalVelocity.resize( 3 );
            thisSpecies->thermalMomentum.resize( 3 );
            for( unsigned int i=0; i<3; i++ ) {
                thisSpecies->thermalVelocity[i] = sqrt( 2.*thisSpecies->thermal_boundary_temperature[i]/thisSpecies->mass );
                thisSpecies->thermalMomentum[i] = thisSpecies->thermalVelocity[i];
                // Caution: momentum in SMILEI actually correspond to p/m
                if( thisSpecies->thermalVelocity[i]>0.3 ) {
                    ERROR( "For species '" << species_name << "' Thermalizing BCs require non-relativistic thermal_boundary_temperature" );
                }
            }
        }
        
        // Manage the ionization parameters
        if( thisSpecies->mass > 0 ) {
            thisSpecies->atomic_number = 0;
            PyTools::extract( "atomic_number", thisSpecies->atomic_number, "Species", ispec );
            
            thisSpecies->maximum_charge_state = 0;
            PyTools::extract( "maximum_charge_state", thisSpecies->maximum_charge_state, "Species", ispec );
            
            std::string model;
            if( PyTools::extract( "ionization_model", model, "Species", ispec ) && model!="none" ) {
            
                thisSpecies->ionization_model = model;
                
                if( ! PyTools::extract( "ionization_electrons", thisSpecies->ionization_electrons, "Species", ispec ) ) {
                    ERROR( "For species '" << species_name << " undefined ionization_electrons (required for ionization)" );
                }
                
                if( ( thisSpecies->atomic_number==0 )&&( thisSpecies->maximum_charge_state==0 ) ) {
                    ERROR( "For species '" << species_name << " undefined atomic_number & maximum_charge_state (required for ionization)" );
                }
                
                if( ( thisSpecies->ionization_model == "from_rate" ) && ( thisSpecies->maximum_charge_state == 0 ) ) {
                    thisSpecies->maximum_charge_state = thisSpecies->atomic_number;
                    WARNING( "For species '" << species_name << " ionization 'from_rate' is used with maximum_charge_state = "<<thisSpecies->maximum_charge_state << " taken from atomic_number" );
                }
                
                if( thisSpecies->ionization_model == "from_rate" ) {
                    thisSpecies->ionization_rate = PyTools::extract_py( "ionization_rate", "Species", ispec );
                    if( thisSpecies->ionization_rate==Py_None ) {
                        ERROR( "For species '" << species_name << " ionization 'from_rate' requires 'ionization_rate' " );
                    } else {
#ifdef SMILEI_USE_NUMPY
                        PyTools::setIteration( 0 );
                        // Test the ionization_rate function with temporary, "fake" particles
                        std::ostringstream name( "" );
                        name << " ionization_rate:";
                        double *dummy = NULL;
                        ParticleData test( params.nDim_particle, thisSpecies->ionization_rate, name.str(), dummy );
#else
                        ERROR( "For species '" << species_name << " ionization 'from_rate' requires Numpy package " );
#endif
                    }
                }
            }
            
            
        }
        
        // Extract if the species is relativistic and needs ad hoc fields initialization
        bool relativistic_field_initialization = false;
        if( !PyTools::extract( "relativistic_field_initialization", relativistic_field_initialization, "Species", ispec ) ) {
            if( patch->isMaster() ) {
                WARNING( "For species '" << species_name << "', relativistic_field_initialization not defined: assumed = 'false'." );
            }
        }
        thisSpecies->relativistic_field_initialization = relativistic_field_initialization;
        
        
        
        // Species geometry
        // ----------------
        
        // Density
        bool ok1, ok2;
        PyObject *profile1( nullptr ), *profile2( nullptr ), *profile3( nullptr );
        
        
        if( thisSpecies->position_initialization_array == NULL ) {
            //These quantities are disregarded if positioning of the species is directly specified by the user
            // Matter particles
            if( thisSpecies->mass > 0 ) {
                ok1 = PyTools::extract_pyProfile( "number_density", profile1, "Species", ispec );
                ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "Species", ispec );
                if( ok1 &&  ok2 ) {
                    ERROR( "For species '" << species_name << "', cannot define both `number_density ` and `charge_density`." );
                }
                if( !ok1 && !ok2 ) {
                    ERROR( "For species '" << species_name << "', must define `number_density ` or `charge_density`." );
                }
                if( ok1 ) {
                    thisSpecies->densityProfileType = "nb";
                }
                if( ok2 ) {
                    thisSpecies->densityProfileType = "charge";
                }
                //MESSAGE(thisSpecies->densityProfileType);
            }
            // Photons
            else if( thisSpecies->mass == 0 ) {
                ok1 = PyTools::extract_pyProfile( "number_density", profile1, "Species", ispec );
                ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "Species", ispec );
                if( ok2 ) {
                    ERROR( "For photon species '" << species_name << "', charge_density has no meaning." );
                }
                if( !ok1 ) {
                    ERROR( "For photon species '" << species_name << "', must define `number_density`." );
                }
                thisSpecies->densityProfileType = "nb";
            }
            
            thisSpecies->densityProfile = new Profile( profile1, params.nDim_field, Tools::merge( thisSpecies->densityProfileType, "_density ", species_name ), true );
            //MESSAGE("creating density profile");
            // Number of particles per cell
            if( !PyTools::extract_pyProfile( "particles_per_cell", profile1, "Species", ispec ) ) {
                ERROR( "For species '" << species_name << "', particles_per_cell not found or not understood" );
            }
            thisSpecies->ppcProfile = new Profile( profile1, params.nDim_field, Tools::merge( "particles_per_cell ", species_name ), true );
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
        thisSpecies->chargeProfile = new Profile( profile1, params.nDim_field, Tools::merge( "charge ", species_name ), true );
        //MESSAGE("creating charge profile");
        if( thisSpecies->momentum_initialization_array == NULL ) {
            // Mean velocity
            if( PyTools::extract3Profiles( "mean_velocity", ispec, profile1, profile2, profile3 ) ) {
                thisSpecies->velocityProfile[0] = new Profile( profile1, params.nDim_field, Tools::merge( "mean_velocity[0] ", species_name ), true );
                thisSpecies->velocityProfile[1] = new Profile( profile2, params.nDim_field, Tools::merge( "mean_velocity[1] ", species_name ), true );
                thisSpecies->velocityProfile[2] = new Profile( profile3, params.nDim_field, Tools::merge( "mean_velocity[2] ", species_name ), true );
            }
            //MESSAGE("velocity profile");
            // Temperature
            if( PyTools::extract3Profiles( "temperature", ispec, profile1, profile2, profile3 ) ) {
                thisSpecies->temperatureProfile[0] = new Profile( profile1, params.nDim_field, Tools::merge( "temperature[0] ", species_name ), true );
                thisSpecies->temperatureProfile[1] = new Profile( profile2, params.nDim_field, Tools::merge( "temperature[1] ", species_name ), true );
                thisSpecies->temperatureProfile[2] = new Profile( profile3, params.nDim_field, Tools::merge( "temperature[2] ", species_name ), true );
            } //MESSAGE("TEMPERATURE");
        } else {
            ok1 = PyTools::extract3Profiles( "mean_velocity", ispec, profile1, profile2, profile3 ) ;
            ok2 = PyTools::extract3Profiles( "temperature", ispec, profile1, profile2, profile3 ) ;
            if( ok1 ) {
                ERROR( "For species '" << species_name << "', cannot define both `mean_velocity` and `momentum_initialization` array." );
            }
            if( ok2 ) {
                ERROR( "For species '" << species_name << "', cannot define both `temperature` and `momentum_initialization` array." );
            }
        }
        
        
        // Get info about tracking
        unsigned int ntrack = PyTools::nComponents( "DiagTrackParticles" );
        thisSpecies->particles->tracked = false;
        for( unsigned int itrack=0; itrack<ntrack; itrack++ ) {
            std::string track_species;
            if( PyTools::extract( "species", track_species, "DiagTrackParticles", itrack ) && track_species==species_name ) {
                if( thisSpecies->particles->tracked ) {
                    ERROR( "In this version, species '" << species_name << "' cannot be tracked by two DiagTrackParticles" );
                }
                thisSpecies->particles->tracked  = true;
            }
        }
        
        // Extract test Species flag
        PyTools::extract( "is_test", thisSpecies->particles->is_test, "Species", ispec );
        
        // Verify they don't ionize
        if( thisSpecies->ionization_model!="none" && thisSpecies->particles->is_test ) {
            ERROR( "For species '" << species_name << "' test & ionized is currently impossible" );
        }
        
        // Create the particles
        if( !params.restart ) {
            // does a loop over all cells in the simulation
            // considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]
            thisSpecies->createParticles( params.n_space, params, patch, 0 );
            //MESSAGE(" PARTICLES");
        } else {
            thisSpecies->particles->initialize( 0, params.nDim_particle );
        }
        
        thisSpecies->initOperators( params, patch );
        //MESSAGE("init operators");
        return thisSpecies;
    } // End Species* create()
    
    
    // Method to clone a species from an existing one
    // Note that this must be only called from cloneVector, because additional init is needed
    static Species *clone( Species *species, Params &params, Patch *patch, bool with_particles = true )
    {
    
        // Create new species object
        Species *newSpecies = NULL;
        
        // Boris, Vay or Higuera-Cary
        if( params.vectorization_mode == "off" ) {
            newSpecies = new SpeciesNorm( params, patch );
        }
#ifdef _VECTO
        else if( params.vectorization_mode == "on" ) {
            newSpecies = new SpeciesNormV( params, patch );
        } else if( params.vectorization_mode == "adaptive_mixed_sort" ) {
            newSpecies = new SpeciesVAdaptive( params, patch );
        } else if( params.vectorization_mode == "adaptive" ) {
            newSpecies = new SpeciesVAdaptiveMixedSort( params, patch );
        }
#endif
        
        // Copy members
        newSpecies->name                                     = species->name;
        newSpecies->pusher                                   = species->pusher;
        newSpecies->radiation_model                          = species->radiation_model;
        newSpecies->radiation_photon_species                 = species->radiation_photon_species;
        newSpecies->radiation_photon_sampling_                = species->radiation_photon_sampling_;
        newSpecies->radiation_photon_gamma_threshold_         = species->radiation_photon_gamma_threshold_;
        newSpecies->photon_species                           = species->photon_species;
        newSpecies->speciesNumber                            = species->speciesNumber;
        newSpecies->position_initialization_on_species       = species->position_initialization_on_species;
        newSpecies->position_initialization_on_species_index = species->position_initialization_on_species_index;
        newSpecies->position_initialization                  = species->position_initialization;
        newSpecies->position_initialization_array            = species->position_initialization_array;
        newSpecies->n_numpy_particles                        = species->n_numpy_particles            ;
        newSpecies->momentum_initialization                  = species->momentum_initialization;
        newSpecies->momentum_initialization_array            = species->momentum_initialization_array;
        newSpecies->c_part_max                               = species->c_part_max;
        newSpecies->mass                                     = species->mass;
        newSpecies->time_frozen                              = species->time_frozen;
        newSpecies->radiating                                = species->radiating;
        newSpecies->relativistic_field_initialization        = species->relativistic_field_initialization;
        newSpecies->time_relativistic_initialization         = species->time_relativistic_initialization;
        newSpecies->boundary_conditions                      = species->boundary_conditions;
        newSpecies->thermal_boundary_temperature             = species->thermal_boundary_temperature;
        newSpecies->thermal_boundary_velocity                = species->thermal_boundary_velocity;
        newSpecies->thermalVelocity                          = species->thermalVelocity;
        newSpecies->thermalMomentum                          = species->thermalMomentum;
        newSpecies->atomic_number                            = species->atomic_number;
        newSpecies->maximum_charge_state                     = species->maximum_charge_state;
        newSpecies->ionization_rate                          = species->ionization_rate;
        if( newSpecies->ionization_rate!=Py_None ) {
            Py_INCREF( newSpecies->ionization_rate );
        }
        newSpecies->ionization_model                         = species->ionization_model;
        newSpecies->densityProfileType                       = species->densityProfileType;
        newSpecies->vectorized_operators                     = species->vectorized_operators;
        newSpecies->merging_method_                          = species->merging_method_;
        newSpecies->merging_time_selection_                  = species->merging_time_selection_;
        newSpecies->chargeProfile                            = new Profile( species->chargeProfile );
        if( species->densityProfile ) {
            newSpecies->densityProfile                       = new Profile( species->densityProfile );
            newSpecies->ppcProfile                           = new Profile( species->ppcProfile );
        }
        newSpecies->velocityProfile.resize( 3 );
        newSpecies->temperatureProfile.resize( 3 );
        
        if( species->velocityProfile[0] ) {
            newSpecies->velocityProfile[0]                   = new Profile( species->velocityProfile[0] );
            newSpecies->velocityProfile[1]                   = new Profile( species->velocityProfile[1] );
            newSpecies->velocityProfile[2]                   = new Profile( species->velocityProfile[2] );
        }
        if( species->temperatureProfile[0] ) {
            newSpecies->temperatureProfile[0]                = new Profile( species->temperatureProfile[0] );
            newSpecies->temperatureProfile[1]                = new Profile( species->temperatureProfile[1] );
            newSpecies->temperatureProfile[2]                = new Profile( species->temperatureProfile[2] );
        }
        newSpecies->max_charge                               = species->max_charge;
        newSpecies->tracking_diagnostic                      = species->tracking_diagnostic;
        newSpecies->ponderomotive_dynamics                   = species->ponderomotive_dynamics;
        
        if( newSpecies->mass==0 ) {
            newSpecies->multiphoton_Breit_Wheeler[0]         = species->multiphoton_Breit_Wheeler[0];
            newSpecies->multiphoton_Breit_Wheeler[1]         = species->multiphoton_Breit_Wheeler[1];
            newSpecies->mBW_pair_creation_sampling[0]        = species->mBW_pair_creation_sampling[0];
            newSpecies->mBW_pair_creation_sampling[1]        = species->mBW_pair_creation_sampling[1];
        }
        
        newSpecies->particles->is_test                       = species->particles->is_test;
        newSpecies->particles->tracked                       = species->particles->tracked;
        newSpecies->particles->isQuantumParameter            = species->particles->isQuantumParameter;
        newSpecies->particles->isMonteCarlo                  = species->particles->isMonteCarlo;
        
        
        // \todo : NOT SURE HOW THIS BEHAVES WITH RESTART
        if( ( !params.restart ) && ( with_particles ) ) {
            newSpecies->createParticles( params.n_space, params, patch, 0 );
        } else {
            newSpecies->particles->initialize( 0, ( *species->particles ) );
        }
        
        newSpecies->initOperators( params, patch );
        
        return newSpecies;
    } // End Species* clone()
    
    
    static std::vector<Species *> createVector( Params &params, Patch *patch )
    {
        // this will be returned
        std::vector<Species *> retSpecies;
        retSpecies.resize( 0 );
        
        // read from python namelist
        unsigned int tot_species_number = PyTools::nComponents( "Species" );
        for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {
            Species *thisSpecies = SpeciesFactory::create( params, ispec, patch );
            
            // Put the newly created species in the vector of species
            retSpecies.push_back( thisSpecies );
            
        }
        
        // Loop species to find species which their particles positions is on another species
        for( unsigned int ispec1 = 0; ispec1<retSpecies.size(); ispec1++ ) {
            if( retSpecies[ispec1]->position_initialization_on_species==true ) {
                // If true then position_initialization of spec1 is not 'centered', 'regular' or 'random'
                // So we have to check if :
                // - 'position_initialization' of spec1 is another already created specie name;
                // - 'position_initialization' of spec1 is not the spec1 name;
                // - 'position_initialization' of spec2 is centered,regular,random;
                // - The number of particle of spec1 is equal to spec2
                
                // Loop all other species
                for( unsigned int ispec2 = 0; ispec2<retSpecies.size(); ispec2++ ) {
                    if( retSpecies[ispec1]->position_initialization == retSpecies[ispec2]->name ) {
                        if( retSpecies[ispec1]->position_initialization==retSpecies[ispec1]->name ) {
                            ERROR( "For species '"<<retSpecies[ispec1]->name<<"' position_initialization must be different from '"<<retSpecies[ispec1]->name<<"'." );
                        }
                        if( retSpecies[ispec2]->position_initialization_on_species==true ) {
                            ERROR( "For species '"<<retSpecies[ispec2]->name<<"' position_initialization must be 'centered', 'regular' or 'random' (pre-defined position) in order to attach '"<<retSpecies[ispec1]->name<<"' to its initial position." );
                        }
                        if( retSpecies[ispec1]->getNbrOfParticles() != retSpecies[ispec2]->getNbrOfParticles() ) {
                            ERROR( "Number of particles in species '"<<retSpecies[ispec1]->name<<"' is not equal to the number of particles in species '"<<retSpecies[ispec2]->name<<"'." );
                        }
                        // We copy ispec2 which is the index of the species, already created, on which initialize particle of the new created species
                        retSpecies[ispec1]->position_initialization_on_species_index=ispec2;
                        // We copy position of species 2 (index ispec2), for position on species 1 (index ispec1)
                        retSpecies[ispec1]->particles->Position=retSpecies[ispec2]->particles->Position;
                    }
                }
                if( retSpecies[ispec1]->position_initialization_on_species_index==-1 ) {
                    ERROR( "Specie '"<<retSpecies[ispec1]->position_initialization<<"' doesn't exist. We can't initialize position on this species. Choose an already created specie or 'centered', 'regular', 'random'." )
                }
            }
        }
        
        // Loop species to find the electron species for ionizable species
        for( unsigned int ispec1 = 0; ispec1<retSpecies.size(); ispec1++ ) {
            if( ! retSpecies[ispec1]->Ionize ) {
                continue;
            }
            
            // Loop all other species
            for( unsigned int ispec2 = 0; ispec2<retSpecies.size(); ispec2++ ) {
                if( retSpecies[ispec1]->ionization_electrons == retSpecies[ispec2]->name ) {
                    if( ispec1==ispec2 ) {
                        ERROR( "For species '"<<retSpecies[ispec1]->name<<"' ionization_electrons must be a distinct species" );
                    }
                    if( retSpecies[ispec2]->mass!=1 ) {
                        ERROR( "For species '"<<retSpecies[ispec1]->name<<"' ionization_electrons must be a species with mass==1" );
                    }
                    retSpecies[ispec1]->electron_species_index = ispec2;
                    retSpecies[ispec1]->electron_species = retSpecies[ispec2];
                    retSpecies[ispec1]->Ionize->new_electrons.tracked = retSpecies[ispec1]->electron_species->particles->tracked;
                    retSpecies[ispec1]->Ionize->new_electrons.initialize( 0, params.nDim_particle );
                    if( ( !retSpecies[ispec1]->getNbrOfParticles() ) && ( !retSpecies[ispec2]->getNbrOfParticles() ) ) {
                        if( retSpecies[ispec1]->atomic_number!=0 ) {
                            int max_eon_number = retSpecies[ispec1]->getNbrOfParticles() * retSpecies[ispec1]->atomic_number;
                            retSpecies[ispec2]->particles->reserve( max_eon_number, retSpecies[ispec2]->particles->dimension() );
                        } else if( retSpecies[ispec1]->maximum_charge_state!=0 ) {
                            int max_eon_number = retSpecies[ispec1]->getNbrOfParticles() * retSpecies[ispec1]->maximum_charge_state;
                            retSpecies[ispec2]->particles->reserve( max_eon_number, retSpecies[ispec2]->particles->dimension() );
                        }
                    }
                    break;
                }
            }
            if( retSpecies[ispec1]->electron_species_index==-1 ) {
                ERROR( "For species '"<<retSpecies[ispec1]->name<<"' ionization_electrons named " << retSpecies[ispec1]->ionization_electrons << " could not be found" );
            }
        }
        
        // Loop species to find the photon species for radiation species
        for( unsigned int ispec1 = 0; ispec1<retSpecies.size(); ispec1++ ) {
            if( ! retSpecies[ispec1]->Radiate ) {
                continue;
            }
            
            // No emission of discrete photon, only scalar diagnostics are updated
            if( retSpecies[ispec1]->radiation_photon_species.empty() ) {
                retSpecies[ispec1]->photon_species_index = -1;
                retSpecies[ispec1]->photon_species = NULL;
            }
            // Else, there will be emission of macro-photons.
            else {
                unsigned int ispec2 = 0;
                for( ispec2 = 0; ispec2<retSpecies.size(); ispec2++ ) {
                    if( retSpecies[ispec1]->radiation_photon_species == retSpecies[ispec2]->name ) {
                        if( ispec1==ispec2 ) {
                            ERROR( "For species '"<<retSpecies[ispec1]->name<<"' radiation_photon_species must be a distinct photon species" );
                        }
                        if( retSpecies[ispec2]->mass!=0 ) {
                            ERROR( "For species '"<<retSpecies[ispec1]->name<<"' radiation_photon_species must be a photon species with mass==0" );
                        }
                        retSpecies[ispec1]->photon_species_index = ispec2;
                        retSpecies[ispec1]->photon_species = retSpecies[ispec2];
                        retSpecies[ispec1]->Radiate->new_photons_.tracked = retSpecies[ispec1]->photon_species->particles->tracked;
                        retSpecies[ispec1]->Radiate->new_photons_.isQuantumParameter = retSpecies[ispec1]->photon_species->particles->isQuantumParameter;
                        retSpecies[ispec1]->Radiate->new_photons_.isMonteCarlo = retSpecies[ispec1]->photon_species->particles->isMonteCarlo;
                        retSpecies[ispec1]->Radiate->new_photons_.initialize( 0,
                                params.nDim_particle );
                        //retSpecies[ispec1]->Radiate->new_photons_.initialize(retSpecies[ispec1]->getNbrOfParticles(),
                        //                                                    params.nDim_particle );
                        retSpecies[ispec2]->particles->reserve( retSpecies[ispec1]->getNbrOfParticles(),
                                                                retSpecies[ispec2]->particles->dimension() );
                        break;
                    }
                }
                if( ispec2 == retSpecies.size() ) {
                    ERROR( "Species '" << retSpecies[ispec1]->radiation_photon_species << "' does not exist." )
                }
            }
        }
        
        // Loop species to find the electron and positron
        // species for multiphoton Breit-wheeler
        for( unsigned int ispec1 = 0; ispec1<retSpecies.size(); ispec1++ ) {
            if( !retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process ) {
                continue;
            } else {
                for( unsigned int ispec2 = 0; ispec2<retSpecies.size(); ispec2++ ) {
                    for( int k=0; k<2; k++ ) {
                        if( retSpecies[ispec1]->multiphoton_Breit_Wheeler[k] == retSpecies[ispec2]->name ) {
                            if( ispec1==ispec2 ) {
                                ERROR( "For species '" << retSpecies[ispec1]->name
                                       << "' pair species must be a distinct particle species" );
                            }
                            if( retSpecies[ispec2]->mass != 1 ) {
                                ERROR( "For species '"<<retSpecies[ispec1]->name<<"' pair species must be an electron and positron species" );
                            }
                            retSpecies[ispec1]->mBW_pair_species_index[k] = ispec2;
                            retSpecies[ispec1]->mBW_pair_species[k] = retSpecies[ispec2];
                            retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].tracked = retSpecies[ispec1]->mBW_pair_species[k]->particles->tracked;
                            retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].isQuantumParameter = retSpecies[ispec1]->mBW_pair_species[k]->particles->isQuantumParameter;
                            retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].isMonteCarlo = retSpecies[ispec1]->mBW_pair_species[k]->particles->isMonteCarlo;
                            retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].initialize( 0,
                                    params.nDim_particle );
                            retSpecies[ispec2]->particles->reserve( retSpecies[ispec1]->getNbrOfParticles(),
                                                                    retSpecies[ispec2]->particles->dimension() );
                        }
                    }
                }
            }
        }
        
        return retSpecies;
    }
    
    // Method to clone the whole vector of species
    static std::vector<Species *> cloneVector( std::vector<Species *> vecSpecies, Params &params, Patch *patch, bool with_particles = true )
    {
        std::vector<Species *> retSpecies;
        retSpecies.resize( 0 );
        
        for( unsigned int ispec = 0; ispec < vecSpecies.size(); ispec++ ) {
            Species *newSpecies = SpeciesFactory::clone( vecSpecies[ispec], params, patch, with_particles );
            retSpecies.push_back( newSpecies );
        }
        
        // Init position on another specie
        for( unsigned int i=0; i<retSpecies.size(); i++ ) {
            if( retSpecies[i]->position_initialization_on_species==true ) {
                unsigned int pos_init_index = retSpecies[i]->position_initialization_on_species_index;
                if( retSpecies[i]->getNbrOfParticles() != retSpecies[pos_init_index]->getNbrOfParticles() ) {
                    ERROR( "Number of particles in species '"<<retSpecies[i]->name<<"' is not equal to the number of particles in species '"<<retSpecies[pos_init_index]->name<<"'." );
                }
                // We copy ispec2 which is the index of the species, already created, on which initialize particles of the new created species
                retSpecies[i]->particles->Position=retSpecies[pos_init_index]->particles->Position;
            }
        }
        
        // Ionization
        for( unsigned int i=0; i<retSpecies.size(); i++ ) {
            if( retSpecies[i]->Ionize ) {
                retSpecies[i]->electron_species_index = vecSpecies[i]->electron_species_index;
                retSpecies[i]->electron_species = retSpecies[retSpecies[i]->electron_species_index];
                retSpecies[i]->Ionize->new_electrons.tracked = retSpecies[i]->electron_species->particles->tracked;
                retSpecies[i]->Ionize->new_electrons.initialize( 0, params.nDim_particle );
            }
        }
        
        // Synchrotron-like radiation
        for( unsigned int i=0; i<retSpecies.size(); i++ ) {
            if( retSpecies[i]->Radiate ) {
                retSpecies[i]->radiation_photon_species = vecSpecies[i]->radiation_photon_species;
                retSpecies[i]->photon_species_index = vecSpecies[i]->photon_species_index;
                if( vecSpecies[i]->photon_species ) {
                    retSpecies[i]->photon_species = retSpecies[retSpecies[i]->photon_species_index];
                    retSpecies[i]->Radiate->new_photons_.tracked = retSpecies[i]->photon_species->particles->tracked;
                    retSpecies[i]->Radiate->new_photons_.isQuantumParameter = retSpecies[i]->photon_species->particles->isQuantumParameter;
                    retSpecies[i]->Radiate->new_photons_.isMonteCarlo = retSpecies[i]->photon_species->particles->isMonteCarlo;
                    //retSpecies[i]->Radiate->new_photons_.initialize(retSpecies[i]->getNbrOfParticles(),
                    //                                               params.nDim_particle );
                    retSpecies[i]->Radiate->new_photons_.initialize( 0, params.nDim_particle );
                } else {
                    retSpecies[i]->photon_species = NULL;
                }
            }
        }
        
        // multiphoton Breit-Wheeler
        for( unsigned int i=0; i<retSpecies.size(); i++ ) {
            if( retSpecies[i]->Multiphoton_Breit_Wheeler_process ) {
                // Loop on pairs
                for( int k=0; k<2; k++ ) {
                    retSpecies[i]->multiphoton_Breit_Wheeler[k] = vecSpecies[i]->multiphoton_Breit_Wheeler[k];
                    retSpecies[i]->mBW_pair_species_index[k] = vecSpecies[i]->mBW_pair_species_index[k];
                    retSpecies[i]->mBW_pair_species[k] = retSpecies[retSpecies[i]->mBW_pair_species_index[k]];
                    retSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].tracked = retSpecies[i]->mBW_pair_species[k]->particles->tracked;
                    retSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].isQuantumParameter = retSpecies[i]->mBW_pair_species[k]->particles->isQuantumParameter;
                    retSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].isMonteCarlo = retSpecies[i]->mBW_pair_species[k]->particles->isMonteCarlo;
                    retSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].initialize(
                        0, params.nDim_particle );
                }
            } else {
                retSpecies[i]->mBW_pair_species[0] = NULL;
                retSpecies[i]->mBW_pair_species[1] = NULL;
            }
        }
        
        return retSpecies;
    }
    
};

#endif
