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


class ParticleInjectorFactory
{
public:
    
    //! Method to create an injector from the namelist parameters
    static ParticleInjector *create( Params &params, Patch *patch, std::vector<Species *> vecSpecies, int injector_index )
    {
        // Create injector object
        ParticleInjector *thisParticleInjector = NULL;
        
        if( patch->isMaster() ) {
            MESSAGE( 1, "" );
        }

        // Read the name of the interpolator
        std::string injector_name( "" );
        PyTools::extract( "name", injector_name, "ParticleInjector", injector_index );
        if( patch->isMaster() ) {
            MESSAGE( 1, "Creating Injector: " << injector_name );
        }
        
        // If no name, name it from the number
        unsigned int tot_injector_number = PyTools::nComponents( "ParticleInjector" );
        if( injector_name.empty() ) {
            std::ostringstream name( "" );
            name << "injector" << std::setfill( '0' ) << std::setw( log10( tot_injector_number )+1 ) << injector_index;
            injector_name=name.str();
            if( patch->isMaster() ) {
                MESSAGE( "For particle injector #" << injector_index << ", name will be " << injector_name );
            }
        }
        
        // Read the species name associated to this interpolator
        std::string species_name( "" );
        PyTools::extract( "species", species_name, "ParticleInjector", injector_index );
        
        // Checks that the species exists and assigns the species index
        bool species_defined = false;
        unsigned int species_number = 0;
        while( (species_number < vecSpecies.size()) and (!species_defined)) {
            if (species_name == vecSpecies[species_number]->name) {
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
        if (params.nDim_field == 1) {
            if ( box_side !="xmin" && box_side !="xmax") {
                ERROR( "For particle injector "<< injector_name << " (# " << injector_index
                << "), box side must be 'xmin' or 'xmax'.");
            }
        }
        if (params.nDim_field == 2) {
            if ( box_side !="xmin" && box_side !="xmax" && box_side !="ymin" && box_side !="ymax") {
                ERROR( "For particle injector "<< injector_name << " (# " << injector_index
                << "), box side must be 'xmin', 'xmax', 'ymin' or 'ymax'.");
            }
        }
        if( patch->isMaster() ) {
            MESSAGE( 2, "> Injection from from the side: " << box_side);
        }
            
        thisParticleInjector = new ParticleInjector(params, patch);
        
        thisParticleInjector->name_ = injector_name;
        thisParticleInjector->species_name_ = species_name;
        thisParticleInjector->species_number_ = species_number;
        thisParticleInjector->box_side_ = box_side;
        thisParticleInjector->position_initialization_on_injector_= false;

        // Read the position initialization
        PyTools::extract( "position_initialization", thisParticleInjector->position_initialization_, "ParticleInjector", injector_index );
        if ( thisParticleInjector->position_initialization_=="species" || thisParticleInjector->position_initialization_=="") {
            MESSAGE( 2, "> Position initialization defined as the species.");
            thisParticleInjector->position_initialization_ = vecSpecies[thisParticleInjector->species_number_]->position_initialization;
        } else if( ( thisParticleInjector->position_initialization_!="regular" )
                   &&( thisParticleInjector->position_initialization_!="random" )
                   &&( thisParticleInjector->position_initialization_!="centered" ) ) {
            thisParticleInjector->position_initialization_on_injector_=true;
            //ERROR("For particle injector " << injector_name << ", position initialization not or badly specified.");
        }

        if( patch->isMaster() ) {
            MESSAGE( 2, "> Position initialization: " << thisParticleInjector->position_initialization_);
        }

        return thisParticleInjector;
        
    }
    
    //! Create the vector of particle injectors
    static std::vector<ParticleInjector *> createVector( Params &params, Patch *patch, std::vector<Species *> vecSpecies )
    {
        // this will be returned
        std::vector<ParticleInjector *> vecParticleInjector;
        vecParticleInjector.resize( 0 );
        
        // read from python namelist
        unsigned int tot_injector_number = PyTools::nComponents( "ParticleInjector" );
        for( unsigned int i_inj = 0; i_inj < tot_injector_number; i_inj++ ) {
            ParticleInjector *thisParticleInjector = ParticleInjectorFactory::create( params, patch, vecSpecies, i_inj );
            // Verify the new injector does not have the same name as a previous one
            for( unsigned int i = 0; i < i_inj; i++ ) {
                if( thisParticleInjector->name_ == vecParticleInjector[i]->name_ ) {
                    ERROR("Two particle injectors cannot have the same name `"<<thisParticleInjector->name_<<"`");
                }
            }
            // Put the newly created injector in the vector of injectors
            vecParticleInjector.push_back( thisParticleInjector );
        }
        
        // Prepare injector with position initialization from another injector
        for( unsigned int i_inj = 0; i_inj < vecParticleInjector.size(); i_inj++ ) {
            // if we need another injector to initialize the positions...
            if (vecParticleInjector[i_inj]->position_initialization_on_injector_) {
                if( vecParticleInjector[i_inj]->position_initialization_==vecParticleInjector[i_inj]->name_ ) {
                    ERROR( "For injector '"<<vecParticleInjector[i_inj]->name_<<"' `position_initialization` can not be the same injector." );
                }
                // We look for this injector in the list
                for( unsigned int i = 0; i < vecParticleInjector.size(); i++ ) {
                    if (vecParticleInjector[i]->name_ == vecParticleInjector[i_inj]->position_initialization_) {
                        if( vecParticleInjector[i]->position_initialization_on_injector_ ) {
                            ERROR( "For injector '"<< vecParticleInjector[i]->name_
                                                   << "' position_initialization must be 'centered', 'regular' or 'random' (pre-defined position) in order to attach '"
                                                   << vecParticleInjector[i]->name_<<"' to its initial position." );
                        }
                        // We copy ispec2 which is the index of the species, already created, on which initialize particle of the new created species
                        vecParticleInjector[i_inj]->position_initialization_on_injector_index_=i;
                    }
                }
            }
        }
        
        return vecParticleInjector;
    }
    
    //! Method to clone a particle injector from an existing one
    //! Note that this must be only called from cloneVector, because additional init is needed
    static ParticleInjector *clone( ParticleInjector * particleInjector, Params &params, Patch *patch)
    {
        ParticleInjector * newParticleInjector = new ParticleInjector(params, patch);
        
        newParticleInjector->name_            = particleInjector->name_;
        newParticleInjector->injector_number_ = particleInjector->injector_number_;
        newParticleInjector->species_name_    = particleInjector->species_name_;
        newParticleInjector->species_number_  = particleInjector->species_number_;
        newParticleInjector->box_side_        = particleInjector->box_side_;
        newParticleInjector->position_initialization_                   = particleInjector->position_initialization_;
        newParticleInjector->position_initialization_on_injector_       = particleInjector->position_initialization_on_injector_;
        newParticleInjector->position_initialization_on_injector_index_ = particleInjector->position_initialization_on_injector_index_;
        
        return newParticleInjector;
        
    }

    //! Method to clone the whole vector
    static std::vector<ParticleInjector *> cloneVector(std::vector<ParticleInjector *> vecParticleInjector, Params &params, Patch *patch )
    {
        
        std::vector<ParticleInjector *> newVecParticleInjector;
        newVecParticleInjector.resize( 0 );
        
        for( unsigned int i_inj = 0; i_inj < vecParticleInjector.size(); i_inj++ ) {
            ParticleInjector *newParticleInjector = ParticleInjectorFactory::clone( vecParticleInjector[i_inj], params, patch );
            newVecParticleInjector.push_back( newParticleInjector );
        }
        
        return newVecParticleInjector;
    }
};

#endif
