#ifndef SPECIESFACTORY_H
#define SPECIESFACTORY_H

#include "Species.h"
#include "Species_norm.h"
#include "Species_rrll.h"

#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "PartBoundCond.h"

#include "Params.h"
#include "Patch.h"

#include "Tools.h"


class SpeciesFactory {
public:
    static Species* create(Params& params, int ispec, Patch* patch) {
        
        std::string species_name("");
        
        PyTools::extract("name",species_name,"Species",ispec);
        if (patch->isMaster()) MESSAGE(1, "Creating Species : " << species_name );
        
        unsigned int tot_species_number = PyTools::nComponents("Species");
        if(species_name.empty()) {
            std::ostringstream name("");
            name << "species" << std::setfill('0') << std::setw(log10(tot_species_number)+1) << ispec;
            species_name=name.str();
            if (patch->isMaster() ) MESSAGE("For species #" << ispec << ", name will be " << species_name);
        }
        
        // Extract type of species dynamics from namelist
        std::string pusher = "boris"; // default value
        if (!PyTools::extract("pusher", pusher ,"Species",ispec) )
            if ( patch->isMaster() ) WARNING("For species '" << species_name << "' pusher not defined: assumed = 'boris'.");
        
        // Create species object
        Species * thisSpecies=NULL;
        if (pusher=="boris" || pusher == "borisnr") {
             // Species with Boris (relativistic =='boris', nonrelativistic=='borisnr') dynamics
             thisSpecies = new Species_norm(params, patch);
        } else if (pusher=="vay") {
             // Species with J.L. Vay dynamics
             thisSpecies = new Species_norm(params, patch);
        } else if (pusher=="higueracary") {
             // Species with Higuary Cary dynamics
             thisSpecies = new Species_norm(params, patch);
        } else if (pusher=="rrll") {
             // Species with Boris dynamics + Radiation Back-Reaction (using the Landau-Lifshitz formula)
             thisSpecies = new Species_rrll(params, patch);
        } else {
            ERROR("For species `" << species_name << " pusher must be 'boris', 'borisnr', 'vay', 'higueracary' or 'rrll'")
        }
        
        thisSpecies->name = species_name;
        thisSpecies->pusher = pusher;
        thisSpecies->speciesNumber = ispec;
        
        // Extract various parameters from the namelist
        
        PyTools::extract("position_initialization",thisSpecies->position_initialization ,"Species",ispec);
        if (thisSpecies->position_initialization.empty()) {
            ERROR("For species '" << species_name << "' empty position_initialization");
        } else if ( (thisSpecies->position_initialization!="regular" )
                  &&(thisSpecies->position_initialization!="random"  )
                  &&(thisSpecies->position_initialization!="centered") ) {
            ERROR("For species '" << species_name << "' unknown position_initialization: " << thisSpecies->position_initialization);
        }
        
        PyTools::extract("momentum_initialization",thisSpecies->momentum_initialization ,"Species",ispec);
        if ( (thisSpecies->momentum_initialization=="mj") || (thisSpecies->momentum_initialization=="maxj") ) {
            thisSpecies->momentum_initialization="maxwell-juettner";
        }
        if (   (thisSpecies->momentum_initialization!="cold")
               && (thisSpecies->momentum_initialization!="maxwell-juettner")
               && (thisSpecies->momentum_initialization!="rectangular") ) {
            ERROR("For species '" << species_name << "' unknown momentum_initialization: "<<thisSpecies->momentum_initialization);
        }
        
        PyTools::extract("c_part_max",thisSpecies->c_part_max,"Species",ispec);
        
        if( !PyTools::extract("mass",thisSpecies->mass ,"Species",ispec) ) {
            ERROR("For species '" << species_name << "' mass not defined.");
        }
        
        PyTools::extract("time_frozen",thisSpecies->time_frozen ,"Species",ispec);
        if (thisSpecies->time_frozen > 0 && thisSpecies->momentum_initialization!="cold") {
            if ( patch->isMaster() ) WARNING("For species '" << species_name << "' possible conflict between time-frozen & not cold initialization");
        }
        
        PyTools::extract("radiating",thisSpecies->radiating ,"Species",ispec);
        if (thisSpecies->pusher=="rrll" && (!thisSpecies->radiating)) {
            if ( patch->isMaster() ) WARNING("For species '" << species_name << "', pusher='rrll' forcing radiating=True");
            thisSpecies->radiating=true;
        }
        
        if( !PyTools::extract("boundary_conditions", thisSpecies->boundary_conditions, "Species", ispec)  )
            ERROR("For species '" << species_name << "', boundary_conditions not defined" );
        
        if( thisSpecies->boundary_conditions.size() == 0 ) {
            ERROR("For species '" << species_name << "', boundary_conditions cannot be empty");
        } else if( thisSpecies->boundary_conditions.size() == 1 ) {
            while( thisSpecies->boundary_conditions.size() < params.nDim_particle )
                thisSpecies->boundary_conditions.push_back( thisSpecies->boundary_conditions[0] );
        } else if( thisSpecies->boundary_conditions.size() != params.nDim_particle ) {
            ERROR("For species '" << species_name << "', boundary_conditions must be the same size as the number of dimensions");
        }
        
        bool has_thermalize = false;
        for( unsigned int iDim=0; iDim<params.nDim_particle; iDim++ ) {
            if( thisSpecies->boundary_conditions[iDim].size() == 1 )
                thisSpecies->boundary_conditions[iDim].push_back( thisSpecies->boundary_conditions[iDim][0] );
            if( thisSpecies->boundary_conditions[iDim].size() != 2 )
                ERROR("For species '" << species_name << "', boundary_conditions["<<iDim<<"] must have one or two arguments")
            if( thisSpecies->boundary_conditions[iDim][0] == "thermalize"
             || thisSpecies->boundary_conditions[iDim][1] == "thermalize" )
                has_thermalize = true;
        }
        
        // for thermalizing BCs on particles check if thermal_boundary_temperature is correctly defined
        bool has_temperature = PyTools::extract("thermal_boundary_temperature",thisSpecies->thermal_boundary_temperature,"Species",ispec);
        bool has_velocity    = PyTools::extract("thermal_boundary_velocity",thisSpecies->thermal_boundary_velocity,"Species",ispec);
        if ( has_thermalize ) {
            if (!has_temperature)
                ERROR("For species '" << species_name << "' thermal_boundary_temperature needs to be defined due to thermalizing BC");
            if (!has_velocity)
                ERROR("For species '" << species_name << "' thermal_boundary_velocity needs to be defined due to thermalizing BC");
            
            if (thisSpecies->thermal_boundary_temperature.size()==1) {
                WARNING("For species '" << species_name << "' Using thermal_boundary_temperature[0] in all directions");
                thisSpecies->thermal_boundary_temperature.resize(3);
                thisSpecies->thermal_boundary_temperature[1] = thisSpecies->thermal_boundary_temperature[0];
                thisSpecies->thermal_boundary_temperature[2] = thisSpecies->thermal_boundary_temperature[0];
            }
            
            // Compute the thermalVelocity & Momentum for thermalizing bcs
            thisSpecies->thermalVelocity.resize(3);
            thisSpecies->thermalMomentum.resize(3);
            for (unsigned int i=0; i<3; i++) {
                thisSpecies->thermalVelocity[i] = sqrt(2.*thisSpecies->thermal_boundary_temperature[i]/thisSpecies->mass);
                thisSpecies->thermalMomentum[i] = thisSpecies->thermalVelocity[i];
                // Caution: momentum in SMILEI actually correspond to p/m
                if (thisSpecies->thermalVelocity[i]>0.3)
                    ERROR("For species '" << species_name << "' Thermalizing BCs require non-relativistic thermal_boundary_temperature");
            }
        }
        
        // Manage the ionization parameters
        thisSpecies->atomic_number = 0;
        PyTools::extract("atomic_number", thisSpecies->atomic_number, "Species",ispec);
        
        std::string model;
        if( PyTools::extract("ionization_model", model, "Species",ispec) && model!="none" ) {
        
            thisSpecies->ionization_model = model;
            
            if( ! PyTools::extract("ionization_electrons", thisSpecies->ionization_electrons, "Species",ispec) ) {
                ERROR("For species '" << species_name << "' undefined ionization_electrons (required for ionization)");
            }
            
            if( thisSpecies->atomic_number==0 ) {
                ERROR("For species '" << species_name << "' undefined atomic_number (required for ionization)");
            }
        }
        
        // Species geometry
        // ----------------
        
        // Density
        bool ok1, ok2;
        PyObject *profile1, *profile2, *profile3;
        ok1 = PyTools::extract_pyProfile("nb_density"    , profile1, "Species", ispec);
        ok2 = PyTools::extract_pyProfile("charge_density", profile1, "Species", ispec);
        if(  ok1 &&  ok2 ) ERROR("For species '" << species_name << "', cannot define both `nb_density` and `charge_density`.");
        if( !ok1 && !ok2 ) ERROR("For species '" << species_name << "', must define `nb_density` or `charge_density`.");
        if( ok1 ) thisSpecies->densityProfileType = "nb";
        if( ok2 ) thisSpecies->densityProfileType = "charge";
        
        thisSpecies->densityProfile = new Profile(profile1, params.nDim_particle, thisSpecies->densityProfileType+"_density "+species_name, true);
        
        // Number of particles per cell
        if( !PyTools::extract_pyProfile("n_part_per_cell", profile1, "Species", ispec))
            ERROR("For species '" << species_name << "', n_part_per_cell not found or not understood");
        thisSpecies->ppcProfile = new Profile(profile1, params.nDim_particle, "n_part_per_cell "+species_name, true);
        
        // Charge
        if( !PyTools::extract_pyProfile("charge", profile1, "Species", ispec))
            ERROR("For species '" << species_name << "', charge not found or not understood");
        thisSpecies->chargeProfile = new Profile(profile1, params.nDim_particle, "charge "+species_name, true);
        
        // Mean velocity
        PyTools::extract3Profiles("mean_velocity", ispec, profile1, profile2, profile3);
        thisSpecies->velocityProfile[0] = new Profile(profile1, params.nDim_particle, "mean_velocity[0] "+species_name, true);
        thisSpecies->velocityProfile[1] = new Profile(profile2, params.nDim_particle, "mean_velocity[1] "+species_name, true);
        thisSpecies->velocityProfile[2] = new Profile(profile3, params.nDim_particle, "mean_velocity[2] "+species_name, true);
        
        // Temperature
        PyTools::extract3Profiles("temperature", ispec, profile1, profile2, profile3);
        thisSpecies->temperatureProfile[0] = new Profile(profile1, params.nDim_particle, "temperature[0] "+species_name, true);
        thisSpecies->temperatureProfile[1] = new Profile(profile2, params.nDim_particle, "temperature[1] "+species_name, true);
        thisSpecies->temperatureProfile[2] = new Profile(profile3, params.nDim_particle, "temperature[2] "+species_name, true);
        
        // Get info about tracking
        unsigned int ntrack = PyTools::nComponents("DiagTrackParticles");
        thisSpecies->particles->tracked = false;
        for( unsigned int itrack=0; itrack<ntrack; itrack++ ) {
            std::string track_species;
            if( PyTools::extract("species", track_species, "DiagTrackParticles", itrack) && track_species==species_name ) {
                if( thisSpecies->particles->tracked )
                    ERROR("In this version, species '" << species_name << "' cannot be tracked by two DiagTrackParticles");
                thisSpecies->particles->tracked  = true;
            }
        }
        
        // Extract test Species flag
        PyTools::extract("is_test", thisSpecies->particles->is_test, "Species", ispec);
        
        // Verify they don't ionize
        if (thisSpecies->ionization_model!="none" && thisSpecies->particles->is_test) {
            ERROR("For species '" << species_name << "' test & ionized is currently impossible");
        }
        
        // Create the particles
        if (!params.restart) {
            // does a loop over all cells in the simulation
            // considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]
            thisSpecies->createParticles(params.n_space, params, patch, 0 );
            
        }
        else
            thisSpecies->particles->initialize( 0, params.nDim_particle );
        
        thisSpecies->initOperators(params, patch);
        
        return thisSpecies;
    } // End Species* create()

    
    // Method to clone a species from an existing one
    // Note that this must be only called from cloneVector, because additional init is needed
    static Species* clone(Species* species, Params &params, Patch* patch, bool with_particles = true) {
        // Create new species object
        Species * newSpecies = NULL;
        if (species->pusher=="boris" 
           || species->pusher=="higueracary"
           || species->pusher=="vay"
           || species->pusher=="borisnr") {
            newSpecies = new Species_norm(params, patch); // Boris
        } else if (species->pusher=="rrll") {
            newSpecies = new Species_rrll(params, patch); // Boris + Radiation Reaction
        }
        // Copy members
        newSpecies->name          = species->name;
        newSpecies->pusher         = species->pusher;
        newSpecies->speciesNumber         = species->speciesNumber;
        newSpecies->position_initialization     = species->position_initialization;
        newSpecies->momentum_initialization     = species->momentum_initialization;
        newSpecies->c_part_max            = species->c_part_max;
        newSpecies->mass                  = species->mass;
        newSpecies->time_frozen           = species->time_frozen;
        newSpecies->radiating             = species->radiating;
        newSpecies->boundary_conditions   = species->boundary_conditions;
        newSpecies->thermal_boundary_temperature  = species->thermal_boundary_temperature;
        newSpecies->thermal_boundary_velocity     = species->thermal_boundary_velocity;
        newSpecies->thermalVelocity       = species->thermalVelocity;
        newSpecies->thermalMomentum       = species->thermalMomentum;
        newSpecies->atomic_number         = species->atomic_number;
        newSpecies->ionization_model      = species->ionization_model;
        newSpecies->densityProfileType    = species->densityProfileType;
        newSpecies->densityProfile        = new Profile(species->densityProfile);
        newSpecies->ppcProfile            = new Profile(species->ppcProfile);
        newSpecies->chargeProfile         = new Profile(species->chargeProfile);
        newSpecies->velocityProfile.resize(3);
        newSpecies->velocityProfile[0]    = new Profile(species->velocityProfile[0]);
        newSpecies->velocityProfile[1]    = new Profile(species->velocityProfile[1]);
        newSpecies->velocityProfile[2]    = new Profile(species->velocityProfile[2]);
        newSpecies->temperatureProfile.resize(3);
        newSpecies->temperatureProfile[0] = new Profile(species->temperatureProfile[0]);
        newSpecies->temperatureProfile[1] = new Profile(species->temperatureProfile[1]);
        newSpecies->temperatureProfile[2] = new Profile(species->temperatureProfile[2]);
        newSpecies->max_charge            = species->max_charge;
        newSpecies->tracking_diagnostic   = species->tracking_diagnostic;
        
        newSpecies->particles->is_test    = species->particles->is_test;
        newSpecies->particles->tracked    = species->particles->tracked;
        
        // \todo : NOT SURE HOW THIS BEHAVES WITH RESTART
        if ( (!params.restart) && (with_particles) ) {
            newSpecies->createParticles(params.n_space, params, patch, 0 );
        }
        else
            newSpecies->particles->initialize( 0, (*species->particles) );

        newSpecies->initOperators(params, patch);
        
        return newSpecies;
    } // End Species* clone()
    
    
    static std::vector<Species*> createVector(Params& params, Patch* patch) {
        // this will be returned
        std::vector<Species*> retSpecies;
        retSpecies.resize(0);
        
        // read from python namelist
        unsigned int tot_species_number = PyTools::nComponents("Species");
        for (unsigned int ispec = 0; ispec < tot_species_number; ispec++) {
            Species* thisSpecies = SpeciesFactory::create(params, ispec, patch);
            
            // Put the newly created species in the vector of species
            retSpecies.push_back(thisSpecies);
            
        }
        
        // Loop species to find the electron species for ionizable species
        for (unsigned int ispec1 = 0; ispec1<retSpecies.size(); ispec1++) {
            if( ! retSpecies[ispec1]->Ionize ) continue;
            
            // Loop all other species
            for (unsigned int ispec2 = 0; ispec2<retSpecies.size(); ispec2++) {
                if( retSpecies[ispec1]->ionization_electrons == retSpecies[ispec2]->name) {
                    if( ispec1==ispec2 )
                        ERROR("For species '"<<retSpecies[ispec1]->name<<"' ionization_electrons must be a distinct species");
                    if (retSpecies[ispec2]->mass!=1)
                        ERROR("For species '"<<retSpecies[ispec1]->name<<"' ionization_electrons must be a species with mass==1");
                    retSpecies[ispec1]->electron_species_index = ispec2;
                    retSpecies[ispec1]->electron_species = retSpecies[ispec2];
                    retSpecies[ispec1]->Ionize->new_electrons.tracked = retSpecies[ispec1]->electron_species->particles->tracked;
                    retSpecies[ispec1]->Ionize->new_electrons.initialize(0, params.nDim_particle );
                    if ( ( !retSpecies[ispec1]->getNbrOfParticles() ) && ( !retSpecies[ispec2]->getNbrOfParticles() ) ) {
                        int max_eon_number = retSpecies[ispec1]->getNbrOfParticles() * retSpecies[ispec1]->atomic_number;
                        retSpecies[ispec2]->particles->reserve( max_eon_number, retSpecies[ispec2]->particles->dimension() );
                    }
                    break;
                }
            }
        }
        
        return retSpecies;
    }
    
    // Method to clone the whole vector of species
    static std::vector<Species*> cloneVector(std::vector<Species*> vecSpecies, Params& params, Patch* patch, bool with_particles = true)
    {
        std::vector<Species*> retSpecies;
        retSpecies.resize(0);
        
        for (unsigned int ispec = 0; ispec < vecSpecies.size(); ispec++) {
            Species* newSpecies = SpeciesFactory::clone(vecSpecies[ispec], params, patch, with_particles);
            retSpecies.push_back( newSpecies );
        }
        
        for (unsigned int i=0; i<retSpecies.size(); i++) {
            if (retSpecies[i]->Ionize) {
                retSpecies[i]->electron_species_index = vecSpecies[i]->electron_species_index;
                retSpecies[i]->electron_species = retSpecies[retSpecies[i]->electron_species_index];
                retSpecies[i]->Ionize->new_electrons.tracked = retSpecies[i]->electron_species->particles->tracked;
                retSpecies[i]->Ionize->new_electrons.initialize(0, params.nDim_particle );
            }
        }
        
        return retSpecies;
    }

};

#endif
