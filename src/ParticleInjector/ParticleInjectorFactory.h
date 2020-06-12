// -----------------------------------------------------------------------------
//
//! \file ParticleInjectorFactory.h
//
//! \brief Contains the class ParticleInjectorFactory that manages the
//!        initialization of the particle injectors.
//
// -----------------------------------------------------------------------------

#ifndef PARTICLEINJECTORFACTORY_H
#define PARTICLEINJECTORFACTORY_H

#include "ParticleInjector.h"

#include "Params.h"
#include "Patch.h"
#include "Species.h"
#include "Tools.h"

class ParticleInjectorFactory
{
public:
    
    //! Method to create an injector from the namelist parameters
    static ParticleInjector *create( Params &params, Patch *patch, std::vector<Species *> species_vector, int injector_index )
    {
        // Create injector object
        ParticleInjector *this_particle_injector = NULL;
        
        if( patch->isMaster() ) {
            MESSAGE( 1, "" );
        }

        // Read the name of the interpolator
        std::string injector_name;
        if( ! PyTools::extractOrNone( "name", injector_name, "ParticleInjector", injector_index ) ) {
            // If no name, name it from the number
            unsigned int tot_injector_number = PyTools::nComponents( "ParticleInjector" );
            std::ostringstream name( "" );
            name << "injector" << std::setfill( '0' ) << std::setw( log10( tot_injector_number )+1 ) << injector_index;
            injector_name=name.str();
        }
        if( patch->isMaster() ) {
            MESSAGE( 1, "Creating Injector: " << injector_name );
        }
        
        // Read the species name associated to this interpolator
        std::string species_name( "" );
        PyTools::extract( "species", species_name, "ParticleInjector", injector_index );
        
        // Checks that the species exists and assigns the species index
        bool species_defined = false;
        unsigned int species_number = 0;
        while( (species_number < species_vector.size()) and (!species_defined)) {
            if (species_name == species_vector[species_number]->name_) {
                species_defined = true;
            } else {
                species_number++;
            }
        }
        if (!species_defined) {
            ERROR( "For particle injector "<< injector_name
            << " (# " << injector_index << "), the specified species does not exist (" << species_name << ")" );
        } else {
            MESSAGE( 2, "> Associated species: " << species_name << " (of index "<< species_number << ")");
        }
        
        // Read the box side associated to this interpolator
        std::string box_side( "" );
        PyTools::extract( "box_side", box_side, "ParticleInjector", injector_index );
        
        // check box side according to the dimension
        if ( box_side !="xmin" && box_side !="xmax") {
            ERROR( "For particle injector "<< injector_name << " (# " << injector_index
                   << "), box side must be 'xmin' or 'xmax'.");
        }

        if( patch->isMaster() ) {
            MESSAGE( 2, "> Injection from from the side: " << box_side);
        }
            
        // Creation of the injector object
        this_particle_injector = new ParticleInjector(params, patch);
        
        this_particle_injector->name_ = injector_name;
        this_particle_injector->species_name_ = species_name;
        this_particle_injector->species_number_ = species_number;
        this_particle_injector->box_side_ = box_side;
        this_particle_injector->position_initialization_on_injector_= false;
        
        // Shortcut pointer to the associated species
        Species * species = species_vector[this_particle_injector->species_number_];

        // Read the position initialization
        PyTools::extract( "position_initialization", this_particle_injector->position_initialization_, "ParticleInjector", injector_index );
        if ( this_particle_injector->position_initialization_=="species" || this_particle_injector->position_initialization_=="") {
            MESSAGE( 2, "> Position initialization defined as the species.");
            this_particle_injector->position_initialization_ = species->position_initialization_;
        } else if( ( this_particle_injector->position_initialization_!="regular" )
                   &&( this_particle_injector->position_initialization_!="random" )
                   &&( this_particle_injector->position_initialization_!="centered" ) ) {
            this_particle_injector->position_initialization_on_injector_=true;
            //ERROR("For particle injector " << injector_name << ", position initialization not or badly specified.");
        }

        if( patch->isMaster() ) {
            MESSAGE( 2, "> Position initialization: " << this_particle_injector->position_initialization_);
        }

        // Read the momentum initialization
        PyTools::extract( "momentum_initialization", this_particle_injector->momentum_initialization_, "ParticleInjector", injector_index );
        if( ( this_particle_injector->momentum_initialization_=="mj" ) || ( this_particle_injector->momentum_initialization_=="maxj" ) ) {
                this_particle_injector->momentum_initialization_="maxwell-juettner";
        }
        if ( this_particle_injector->momentum_initialization_=="species" || this_particle_injector->momentum_initialization_=="") {
            MESSAGE( 2, "> Momentum initialization defined as the species.");
            this_particle_injector->momentum_initialization_ = species->momentum_initialization_;
        }
        // Matter particles
        if( species_vector[this_particle_injector->species_number_]->mass_ > 0 ) {
            if( ( this_particle_injector->momentum_initialization_!="cold" )
                    && ( this_particle_injector->momentum_initialization_!="maxwell-juettner" )
                    && ( this_particle_injector->momentum_initialization_!="rectangular" ) ) {
                ERROR( "For particle injector '" << injector_name
                       << "' unknown momentum_initialization: "
                       <<this_particle_injector->momentum_initialization_ );
            }
        }
        // Photons
        else if( species_vector[this_particle_injector->species_number_]->mass_ == 0 ) {
            if ( this_particle_injector->momentum_initialization_ == "maxwell-juettner" ) {
                ERROR( "For photon injector '" << injector_name
                       << "' Maxwell-Juettner is not valid.");
            }
            if (( this_particle_injector->momentum_initialization_!="cold" )
                    && ( this_particle_injector->momentum_initialization_!="rectangular" ) ) {
                ERROR( "For photon injector '" << injector_name
                       << "' unknown momentum_initialization: "
                       <<this_particle_injector->momentum_initialization_ );
            }
        }

        if( patch->isMaster() ) {
            MESSAGE( 2, "> Momentum initialization: " << this_particle_injector->momentum_initialization_);
        }
        
        // Mean velocity
        // std::vector<double> mean_velocity_input;
        std::vector<PyObject *> prof;
        if( PyTools::extract_1or3Profiles( "mean_velocity", "ParticleInjector" , injector_index, prof ) ) {
            this_particle_injector->velocity_profile_[0] = new Profile( prof[0], params.nDim_field, Tools::merge( "mean_velocity[0] ", this_particle_injector->name_ ), true );
            this_particle_injector->velocity_profile_[1] = new Profile( prof[1], params.nDim_field, Tools::merge( "mean_velocity[1] ", this_particle_injector->name_ ), true );
            this_particle_injector->velocity_profile_[2] = new Profile( prof[2], params.nDim_field, Tools::merge( "mean_velocity[2] ", this_particle_injector->name_ ), true );
            MESSAGE( 2, "> Mean velocity redefined: " << this_particle_injector->velocity_profile_[0]->getInfo());
            // string message =  "> Mean velocity: ";
            // for (unsigned int i = 0 ; i < mean_velocity_input.size()-1 ; i++) {
            //     message += to_string(mean_velocity_input[i]) + ", ";
            // }
            // message += to_string(mean_velocity_input[mean_velocity_input.size()-1]);
            // MESSAGE(2, message);
        } else {
            MESSAGE( 2, "> Mean velocity defined as the species.");
            this_particle_injector->velocity_profile_[0] = new Profile(species->velocity_profile_[0]);
            this_particle_injector->velocity_profile_[1] = new Profile(species->velocity_profile_[1]);
            this_particle_injector->velocity_profile_[2] = new Profile(species->velocity_profile_[2]);
        }
        
        // Temperature
        // std::vector<double> temperature_input;
        if( PyTools::extract_1or3Profiles( "temperature", "ParticleInjector", injector_index, prof ) ) {
            this_particle_injector->temperature_profile_[0] = new Profile( prof[0], params.nDim_field, Tools::merge( "temperature[0] ", this_particle_injector->name_ ), true );
            this_particle_injector->temperature_profile_[1] = new Profile( prof[1], params.nDim_field, Tools::merge( "temperature[1] ", this_particle_injector->name_ ), true );
            this_particle_injector->temperature_profile_[2] = new Profile( prof[2], params.nDim_field, Tools::merge( "temperature[2] ", this_particle_injector->name_ ), true );
            
            // string message =  "> Temperature: ";
            // for (unsigned int i = 0 ; i < temperature_input.size()-1 ; i++) {
            //     message += to_string(temperature_input[i]) + ", ";
            // }
            // message += to_string(temperature_input[temperature_input.size()-1]);
            // MESSAGE(2, message);
        } else {
            MESSAGE( 2, "> Temperature defined as the species.");
            this_particle_injector->temperature_profile_[0] = new Profile(species->temperature_profile_[0]);
            this_particle_injector->temperature_profile_[1] = new Profile(species->temperature_profile_[1]);
            this_particle_injector->temperature_profile_[2] = new Profile(species->temperature_profile_[2]);
        }

        // We read the density profile specific for the injector
        // If nothing specified, similar to the species
        // First for mass particles:
        bool ok1, ok2;
        PyObject * profile1 = nullptr;
        if( species->mass_ > 0 ) {
            ok1 = PyTools::extract_pyProfile( "number_density", profile1, "ParticleInjector", injector_index );
            ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "ParticleInjector", injector_index );
            if (ok1 && ok2) {
                ERROR( "For injector '" << this_particle_injector->name_ << "', cannot define both `number_density ` and `charge_density`." );
            } else if( !ok1 && !ok2 ) {
                this_particle_injector->density_profile_type_ = species->density_profile_type_;
                if (species->density_profile_type_ == "nb") {
                    MESSAGE( 2, "> Number density profile defined as the species.");
                } else if (species->density_profile_type_ == "charge") {
                    MESSAGE( 2, "> Charge density profile defined as the species.");
                }
                this_particle_injector->density_profile_ = species->density_profile_;
            } else {
                if( ok1 ) {
                    this_particle_injector->density_profile_type_ = "nb";
                }
                if( ok2 ) {
                    this_particle_injector->density_profile_type_ = "charge";
                }
                this_particle_injector->density_profile_ =
                new Profile( profile1, params.nDim_field, Tools::merge( this_particle_injector->density_profile_type_, "_density ", injector_name ), true );
            }
        }
        // Photons
        else if( species->mass_ == 0 ) {
            ok1 = PyTools::extract_pyProfile( "number_density", profile1, "ParticleInjector", injector_index );
            ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "ParticleInjector", injector_index );
            if( ok2 ) {
                ERROR( "For photon injector '" << injector_name << "', `charge_density` has no meaning."
                        << "You must use `number_density`." );
            }
            if( !ok1 ) {
                species->density_profile_type_ = "nb";
                this_particle_injector->density_profile_ = species->density_profile_;
                MESSAGE( 2, "> Number density profile defined as the species.");
            } else {
                species->density_profile_type_ = "nb";
                this_particle_injector->density_profile_ =
                new Profile( profile1, params.nDim_field, Tools::merge( this_particle_injector->density_profile_type_, "_density ", injector_name ), true );
            }
        }
        
        PyTools::extract_pyProfile( "time_envelope", profile1, "ParticleInjector", injector_index );
        this_particle_injector->time_profile_ = new Profile( profile1, 1, Tools::merge( "time_profile_", injector_name ) );
        MESSAGE( 2, "> Time profile: " << this_particle_injector->time_profile_->getInfo());

        // Number of particles per cell
        if( !PyTools::extract_pyProfile( "particles_per_cell", profile1, "ParticleInjector", injector_index ) ) {
            this_particle_injector->particles_per_cell_profile_ = species->particles_per_cell_profile_;
            MESSAGE( 2, "> Particles per cell defined as the associated species: "
            << this_particle_injector->particles_per_cell_profile_->getInfo() << ".");
        } else {
            this_particle_injector->particles_per_cell_profile_ = new Profile( profile1, params.nDim_field,
                                          Tools::merge( "particles_per_cell ", injector_name ), true );
            MESSAGE( 2, "> Particles per cell profile: "
            << this_particle_injector->particles_per_cell_profile_->getInfo() << ".");
        }

        return this_particle_injector;
        
    }
    
    //! Create the vector of particle injectors
    static std::vector<ParticleInjector *> createVector( Params &params, Patch *patch, std::vector<Species *> species_vector )
    {
        // this will be returned
        std::vector<ParticleInjector *> particle_injector_vector;
        particle_injector_vector.resize( 0 );
        
        // read from python namelist
        unsigned int tot_injector_number = PyTools::nComponents( "ParticleInjector" );
        if (tot_injector_number > 0) {
            TITLE("Initializing particle injectors")
        }
        for( unsigned int i_inj = 0; i_inj < tot_injector_number; i_inj++ ) {
            ParticleInjector *this_particle_injector = ParticleInjectorFactory::create( params, patch, species_vector, i_inj );
            // Verify the new injector does not have the same name as a previous one
            for( unsigned int i = 0; i < i_inj; i++ ) {
                if( this_particle_injector->name_ == particle_injector_vector[i]->name_ ) {
                    ERROR("Two particle injectors cannot have the same name `"<<this_particle_injector->name_<<"`");
                }
            }
            // Put the newly created injector in the vector of injectors
            particle_injector_vector.push_back( this_particle_injector );
        }
        
        // Prepare injector with position initialization from another injector
        for( unsigned int i_inj = 0; i_inj < particle_injector_vector.size(); i_inj++ ) {
            // if we need another injector to initialize the positions...
            if (particle_injector_vector[i_inj]->position_initialization_on_injector_) {
                if( particle_injector_vector[i_inj]->position_initialization_==particle_injector_vector[i_inj]->name_ ) {
                    ERROR( "For injector '"<<particle_injector_vector[i_inj]->name_<<"' `position_initialization` can not be the same injector." );
                }
                // We look for this injector in the list
                for( unsigned int i = 0; i < particle_injector_vector.size(); i++ ) {
                    if (particle_injector_vector[i]->name_ == particle_injector_vector[i_inj]->position_initialization_) {
                        if( particle_injector_vector[i]->position_initialization_on_injector_ ) {
                            ERROR( "For injector '"<< particle_injector_vector[i]->name_
                                                   << "' position_initialization must be 'centered', 'regular' or 'random' (pre-defined position) in order to attach '"
                                                   << particle_injector_vector[i]->name_<<"' to its initial position." );
                        }
                        // We copy ispec2 which is the index of the species, already created, on which initialize particle of the new created species
                        particle_injector_vector[i_inj]->position_initialization_on_injector_index_=i;
                    }
                }
            }
        }
        
        return particle_injector_vector;
    }
    
    //! Method to clone a particle injector from an existing one
    //! Note that this must be only called from cloneVector, because additional init is needed
    static ParticleInjector *clone( ParticleInjector * particle_injector, Params &params, Patch *patch)
    {
        ParticleInjector * new_particle_injector = new ParticleInjector(params, patch);
        
        new_particle_injector->name_            = particle_injector->name_;
        new_particle_injector->injector_number_ = particle_injector->injector_number_;
        new_particle_injector->species_name_    = particle_injector->species_name_;
        new_particle_injector->species_number_  = particle_injector->species_number_;
        new_particle_injector->box_side_        = particle_injector->box_side_;
        new_particle_injector->position_initialization_                   = particle_injector->position_initialization_;
        new_particle_injector->momentum_initialization_                   = particle_injector->momentum_initialization_;
        new_particle_injector->position_initialization_on_injector_       = particle_injector->position_initialization_on_injector_;
        new_particle_injector->position_initialization_on_injector_index_ = particle_injector->position_initialization_on_injector_index_;
        
        new_particle_injector->velocity_profile_.resize( 3 );

        if( particle_injector->velocity_profile_[0] ) {
            new_particle_injector->velocity_profile_[0]                   = new Profile( particle_injector->velocity_profile_[0] );
            new_particle_injector->velocity_profile_[1]                   = new Profile( particle_injector->velocity_profile_[1] );
            new_particle_injector->velocity_profile_[2]                   = new Profile( particle_injector->velocity_profile_[2] );
        }
        
        new_particle_injector->temperature_profile_.resize( 3 );
        
        if( particle_injector->temperature_profile_[0] ) {
            new_particle_injector->temperature_profile_[0]                = new Profile( particle_injector->temperature_profile_[0] );
            new_particle_injector->temperature_profile_[1]                = new Profile( particle_injector->temperature_profile_[1] );
            new_particle_injector->temperature_profile_[2]                = new Profile( particle_injector->temperature_profile_[2] );
        }
        
        new_particle_injector->density_profile_type_ = particle_injector->density_profile_type_;
        new_particle_injector->density_profile_     = new Profile( particle_injector->density_profile_ );
        
        new_particle_injector->time_profile_     = new Profile( particle_injector->time_profile_ );
        
        new_particle_injector->particles_per_cell_profile_ = new Profile( particle_injector->particles_per_cell_profile_ );
        
        return new_particle_injector;
        
    }

    //! Method to clone the whole vector
    static std::vector<ParticleInjector *> cloneVector(std::vector<ParticleInjector *> particle_injector_vector, Params &params, Patch *patch )
    {
        
        std::vector<ParticleInjector *> new_vector_particle_injector;
        new_vector_particle_injector.resize( 0 );
        
        for( unsigned int i_inj = 0; i_inj < particle_injector_vector.size(); i_inj++ ) {
            ParticleInjector *new_particle_injector = ParticleInjectorFactory::clone( particle_injector_vector[i_inj], params, patch );
            new_vector_particle_injector.push_back( new_particle_injector );
        }
        
        return new_vector_particle_injector;
    }
};

#endif
