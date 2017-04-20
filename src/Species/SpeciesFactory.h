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
        
        std::string species_type("");
        
        PyTools::extract("species_type",species_type,"Species",ispec);
        
        unsigned int tot_species_number = PyTools::nComponents("Species");
        if(species_type.empty()) {
            std::ostringstream name("");
            name << "species" << std::setfill('0') << std::setw(log10(tot_species_number)+1) << ispec;
            species_type=name.str();
            if (patch->isMaster() ) MESSAGE("For species #" << ispec << ", parameter species_type will be " << species_type);
        }
        
        // Extract type of species dynamics from namelist
        std::string dynamics_type = "norm"; // default value
        if (!PyTools::extract("dynamics_type", dynamics_type ,"Species",ispec) )
            if ( patch->isMaster() ) WARNING("For species '" << species_type << "' dynamics_type not defined: assumed = 'norm'.");
        
        // Create species object
        Species * thisSpecies=NULL;
        if (dynamics_type=="norm" || dynamics_type == "borisnr") {
             // Species with Boris (relativistic =='norm', nonrelativistic=='borisnr') dynamics
             thisSpecies = new Species_norm(params, patch);
        } else if (dynamics_type=="vay") {
             // Species with J.L. Vay dynamics
             thisSpecies = new Species_norm(params, patch);
        } else if (dynamics_type=="higueracary") {
             // Species with Higuary Cary dynamics
             thisSpecies = new Species_norm(params, patch);
        } else if (dynamics_type=="rrll") {
             // Species with Boris dynamics + Radiation Back-Reaction (using the Landau-Lifshitz formula)
             thisSpecies = new Species_rrll(params, patch);
        } else {
            ERROR("For species `" << species_type << " dynamics_type must be 'norm', 'borisnr', 'vay', 'higueracary' or 'rrll'")
        }
        
        thisSpecies->species_type = species_type;
        thisSpecies->dynamics_type = dynamics_type;
        thisSpecies->speciesNumber = ispec;
        
        // Extract various parameters from the namelist
        
        PyTools::extract("initPosition_type",thisSpecies->initPosition_type ,"Species",ispec);
        if (thisSpecies->initPosition_type.empty()) {
            ERROR("For species '" << species_type << "' empty initPosition_type");
        } else if ( (thisSpecies->initPosition_type!="regular" )
                  &&(thisSpecies->initPosition_type!="random"  )
                  &&(thisSpecies->initPosition_type!="centered") ) {
            ERROR("For species '" << species_type << "' unknown initPosition_type: " << thisSpecies->initPosition_type);
        }
        
        PyTools::extract("initMomentum_type",thisSpecies->initMomentum_type ,"Species",ispec);
        if ( (thisSpecies->initMomentum_type=="mj") || (thisSpecies->initMomentum_type=="maxj") ) {
            thisSpecies->initMomentum_type="maxwell-juettner";
        }
        if (   (thisSpecies->initMomentum_type!="cold")
               && (thisSpecies->initMomentum_type!="maxwell-juettner")
               && (thisSpecies->initMomentum_type!="rectangular") ) {
            ERROR("For species '" << species_type << "' unknown initMomentum_type: "<<thisSpecies->initMomentum_type);
        }
        
        PyTools::extract("c_part_max",thisSpecies->c_part_max,"Species",ispec);
        
        if( !PyTools::extract("mass",thisSpecies->mass ,"Species",ispec) ) {
            ERROR("For species '" << species_type << "' mass not defined.");
        }
        
        PyTools::extract("time_frozen",thisSpecies->time_frozen ,"Species",ispec);
        if (thisSpecies->time_frozen > 0 && thisSpecies->initMomentum_type!="cold") {
            if ( patch->isMaster() ) WARNING("For species '" << species_type << "' possible conflict between time-frozen & not cold initialization");
        }
        
        PyTools::extract("radiating",thisSpecies->radiating ,"Species",ispec);
        if (thisSpecies->dynamics_type=="rrll" && (!thisSpecies->radiating)) {
            if ( patch->isMaster() ) WARNING("For species '" << species_type << "', dynamics_type='rrll' forcing radiating=True");
            thisSpecies->radiating=true;
        }
        
        if (!PyTools::extract("bc_part_type_xmin",thisSpecies->bc_part_type_xmin,"Species",ispec) )
            ERROR("For species '" << species_type << "', bc_part_type_xmin not defined");
        if (!PyTools::extract("bc_part_type_xmax",thisSpecies->bc_part_type_xmax,"Species",ispec) )
            ERROR("For species '" << species_type << "', bc_part_type_xmax not defined");
        
        if (params.nDim_particle>1) {
            if (!PyTools::extract("bc_part_type_ymin",thisSpecies->bc_part_type_ymin,"Species",ispec) )
                ERROR("For species '" << species_type << "', bc_part_type_ymin not defined");
            if (!PyTools::extract("bc_part_type_ymax",thisSpecies->bc_part_type_ymax,"Species",ispec) )
                ERROR("For species '" << species_type << "', bc_part_type_ymax not defined");
        }
        
        if (params.nDim_particle>2) {
            if (!PyTools::extract("bc_part_type_zmin",thisSpecies->bc_part_type_zmin,"Species",ispec) )
                ERROR("For species '" << species_type << "', bc_part_type_zmin not defined");
            if (!PyTools::extract("bc_part_type_zmax",thisSpecies->bc_part_type_zmax,"Species",ispec) )
                ERROR("For species '" << species_type << "', bc_part_type_zmax not defined");
        }
        
        // for thermalizing BCs on particles check if thermT is correctly defined
        bool thermTisDefined=false;
        bool thermVisDefined=false;
        if ( (thisSpecies->bc_part_type_xmin=="thermalize") || (thisSpecies->bc_part_type_xmax=="thermalize") ){
            thermTisDefined=PyTools::extract("thermT",thisSpecies->thermT,"Species",ispec);
            if (!thermTisDefined) ERROR("For species '" << species_type << "' thermT needs to be defined due to x-BC thermalize");
            thermVisDefined=PyTools::extract("thermVelocity",thisSpecies->thermVelocity,"Species",ispec);
            if (!thermVisDefined) ERROR("For species '" << species_type << "' thermVelocity needs to be defined due to x-BC thermalize");
        }
        if ( (params.nDim_particle==2) && (!thermTisDefined) && (!thermVisDefined) &&
             (thisSpecies->bc_part_type_ymin=="thermalize" || thisSpecies->bc_part_type_ymax=="thermalize") ) {
            thermTisDefined=PyTools::extract("thermT",thisSpecies->thermT,"Species",ispec);
            if (!thermTisDefined) ERROR("For species '" << species_type << "' thermT needs to be defined due to y-BC thermalize");
            thermVisDefined=PyTools::extract("thermVelocity",thisSpecies->thermVelocity,"Species",ispec);
            if (!thermTisDefined) ERROR("For species '" << species_type << "' thermVelocity needs to be defined due to y-BC thermalize");
        }
        if (thermTisDefined) {
            if (thisSpecies->thermT.size()==1) {
                WARNING("For species '" << species_type << "' Using thermT[0] in all directions");
                thisSpecies->thermT.resize(3);
                for (unsigned int i=1; i<3;i++)
                    thisSpecies->thermT[i]=thisSpecies->thermT[0];
            }
        } else {
            thisSpecies->thermT.resize(3);
            for (unsigned int i=0; i<3;i++)
                thisSpecies->thermT[i]=0.0;
            thisSpecies->thermVelocity.resize(3);
            for (unsigned int i=0; i<3;i++)
                thisSpecies->thermVelocity[i]=0.0;
        }
        
        // Compute the thermalVelocity & Momentum for thermalizing bcs
        thisSpecies->thermalVelocity.resize(3);
        thisSpecies->thermalMomentum.resize(3);
        
        for (unsigned int i=0; i<3; i++) {
            thisSpecies->thermalVelocity[i] = sqrt(2.*thisSpecies->thermT[i]/thisSpecies->mass);
            thisSpecies->thermalMomentum[i] = thisSpecies->thermalVelocity[i];
            // Caution: momentum in SMILEI actually correspond to p/m
            if (thisSpecies->thermalVelocity[i]>0.3) ERROR("For species '" << species_type << "' Thermalizing BCs require non-relativistic thermT");
        }
        
        
        
        // Manage the ionization parameters
        thisSpecies->atomic_number = 0;
        PyTools::extract("atomic_number", thisSpecies->atomic_number, "Species",ispec);
        
        std::string model;
        if( PyTools::extract("ionization_model", model, "Species",ispec) && model!="none" ) {
        
            thisSpecies->ionization_model = model;
            
            if( ! PyTools::extract("ionization_electrons", thisSpecies->ionization_electrons, "Species",ispec) ) {
                ERROR("For species '" << species_type << "' undefined ionization_electrons (required for ionization)");
            }
            
            if( thisSpecies->atomic_number==0 ) {
                ERROR("For species '" << species_type << "' undefined atomic_number (required for ionization)");
            }
        }
        
        // Species geometry
        // ----------------
        
        // Density
        bool ok1, ok2;
        PyObject *profile1, *profile2, *profile3;
        ok1 = PyTools::extract_pyProfile("nb_density"    , profile1, "Species", ispec);
        ok2 = PyTools::extract_pyProfile("charge_density", profile1, "Species", ispec);
        if(  ok1 &&  ok2 ) ERROR("For species '" << species_type << "', cannot define both `nb_density` and `charge_density`.");
        if( !ok1 && !ok2 ) ERROR("For species '" << species_type << "', must define `nb_density` or `charge_density`.");
        if( ok1 ) thisSpecies->densityProfileType = "nb";
        if( ok2 ) thisSpecies->densityProfileType = "charge";
        
        thisSpecies->densityProfile = new Profile(profile1, params.nDim_particle, thisSpecies->densityProfileType+"_density "+species_type);
        
        // Number of particles per cell
        if( !PyTools::extract_pyProfile("n_part_per_cell", profile1, "Species", ispec))
            ERROR("For species '" << species_type << "', n_part_per_cell not found or not understood");
        thisSpecies->ppcProfile = new Profile(profile1, params.nDim_particle, "n_part_per_cell "+species_type);
        
        // Charge
        if( !PyTools::extract_pyProfile("charge", profile1, "Species", ispec))
            ERROR("For species '" << species_type << "', charge not found or not understood");
        thisSpecies->chargeProfile = new Profile(profile1, params.nDim_particle, "charge "+species_type);
        
        // Mean velocity
        PyTools::extract3Profiles("mean_velocity", ispec, profile1, profile2, profile3);
        thisSpecies->velocityProfile[0] = new Profile(profile1, params.nDim_particle, "mean_velocity[0] "+species_type);
        thisSpecies->velocityProfile[1] = new Profile(profile2, params.nDim_particle, "mean_velocity[1] "+species_type);
        thisSpecies->velocityProfile[2] = new Profile(profile3, params.nDim_particle, "mean_velocity[2] "+species_type);
        
        // Temperature
        PyTools::extract3Profiles("temperature", ispec, profile1, profile2, profile3);
        thisSpecies->temperatureProfile[0] = new Profile(profile1, params.nDim_particle, "temperature[0] "+species_type);
        thisSpecies->temperatureProfile[1] = new Profile(profile2, params.nDim_particle, "temperature[1] "+species_type);
        thisSpecies->temperatureProfile[2] = new Profile(profile3, params.nDim_particle, "temperature[2] "+species_type);
        
        
        // CALCULATE USEFUL VALUES
        
        /*        double gamma=1.+thisSpecies->thermT[0]/thisSpecies->mass;
            
                  for (unsigned int i=0; i<3; i++) {
                  thisSpecies->thermalVelocity[i] = sqrt( 1.-1./gamma*gamma );
                  thisSpecies->thermalMomentum[i] = gamma*thisSpecies->thermalVelocity[i];
                  }
            
                  double gamma=1.+thisSpecies->thermT[0]/thisSpecies->mass;
        */
        
        thisSpecies->thermalVelocity.resize(3);
        thisSpecies->thermalMomentum.resize(3);
        
        if (thermTisDefined) {
            if ( patch->isMaster() ) WARNING("\tFor species '" << species_type << "' Using thermT[0] in all directions");
            if (thisSpecies->thermalVelocity[0]>0.3) {
                ERROR("For species '" << species_type << "' thermalising BCs require ThermT[0]="<<thisSpecies->thermT[0]<<"<<"<<thisSpecies->mass);
            }
            for (unsigned int i=0; i<3; i++) {
                thisSpecies->thermalVelocity[i] = sqrt(2.*thisSpecies->thermT[0]/thisSpecies->mass);
                thisSpecies->thermalMomentum[i] = thisSpecies->thermalVelocity[i];
            }
        }
        
        
        // Extract test Species flag
        PyTools::extract("isTest", thisSpecies->particles->isTest, "Species", ispec);
        
        // Verify they don't ionize
        if (thisSpecies->ionization_model!="none" && thisSpecies->particles->isTest) {
            ERROR("For species '" << species_type << "' test & ionized is currently impossible");
        }
        
        // Find out whether this species is tracked
        TimeSelection track_timeSelection( PyTools::extract_py("track_every", "Species", ispec), "Track" );
        thisSpecies->particles->tracked = ! track_timeSelection.isEmpty();
        
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
        if (species->dynamics_type=="norm" 
           || species->dynamics_type=="higueracary"
           || species->dynamics_type=="vay"
           || species->dynamics_type=="borisnr") {
            newSpecies = new Species_norm(params, patch); // Boris
        } else if (species->dynamics_type=="rrll") {
            newSpecies = new Species_rrll(params, patch); // Boris + Radiation Reaction
        }
        // Copy members
        newSpecies->species_type          = species->species_type;
        newSpecies->dynamics_type         = species->dynamics_type;
        newSpecies->speciesNumber         = species->speciesNumber;
        newSpecies->initPosition_type     = species->initPosition_type;
        newSpecies->initMomentum_type     = species->initMomentum_type;
        newSpecies->c_part_max            = species->c_part_max;
        newSpecies->mass                  = species->mass;
        newSpecies->time_frozen           = species->time_frozen;
        newSpecies->radiating             = species->radiating;
        newSpecies->bc_part_type_xmin     = species->bc_part_type_xmin;
        newSpecies->bc_part_type_xmax     = species->bc_part_type_xmax;
        newSpecies->bc_part_type_ymin    = species->bc_part_type_ymin;
        newSpecies->bc_part_type_ymax    = species->bc_part_type_ymax;
        newSpecies->bc_part_type_zmin   = species->bc_part_type_zmin;
        newSpecies->bc_part_type_zmax       = species->bc_part_type_zmax;
        newSpecies->thermT                = species->thermT;
        newSpecies->thermVelocity         = species->thermVelocity;
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
        
        newSpecies->particles->isTest              = species->particles->isTest;
        newSpecies->particles->tracked             = species->particles->tracked;
        
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
        
        if (patch->isMaster()) MESSAGE(1, "Creating Species :" );
        
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
                if( retSpecies[ispec1]->ionization_electrons == retSpecies[ispec2]->species_type) {
                    if( ispec1==ispec2 )
                        ERROR("For species '"<<retSpecies[ispec1]->species_type<<"' ionization_electrons must be a distinct species");
                    if (retSpecies[ispec2]->mass!=1)
                        ERROR("For species '"<<retSpecies[ispec1]->species_type<<"' ionization_electrons must be a species with mass==1");
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
