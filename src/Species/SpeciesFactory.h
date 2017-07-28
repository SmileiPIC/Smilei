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

        if (patch->isMaster()) MESSAGE(1, "Creating Species : " << species_type );

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


        // Extract type of species radiation from namelist
        std::string radiation_model = "none"; // default value
        if (!PyTools::extract("radiation_model", radiation_model ,"Species",ispec) )
            if ( patch->isMaster() ) WARNING("For species '" << species_type << "' radiation_model not defined: assumed = 'none'.");

        // Extract mass from namelist
        double mass;
        if (!PyTools::extract("mass", mass ,"Species",ispec) )
        {
            if ( patch->isMaster() ) ERROR("For species '" << species_type << "` mass is not defined.");
        }

        // Create species object
        Species * thisSpecies = NULL;

        // Particles
        if (mass > 0.)
        {
            // Dynamics of the species
            if (dynamics_type=="norm"
             || dynamics_type == "borisnr"
             || dynamics_type == "vay"
             || dynamics_type=="higueracary") {
                 // Species with relativistic Boris dynamics if  =='norm'
                 // Species with nonrelativistic Boris dynamics == 'borisnr'
                 // Species with J.L. Vay dynamics if == "vay"
                 // Species with Higuary Cary dynamics if == "higueracary"
                 thisSpecies = new SpeciesNorm(params, patch);
            } else {
                ERROR("For species `" << species_type << "` dynamics_type must be 'norm', 'borisnr', 'vay', 'higueracary'");
            }
            thisSpecies->dynamics_type = dynamics_type;

            // Radiation model of the species
            // Species with a Monte-Carlo process for the radiation loss
            if (radiation_model=="Monte-Carlo") {
                 thisSpecies->particles->isQuantumParameter = true;
                 thisSpecies->particles->isMonteCarlo = true;
                 thisSpecies->radiating = true;
            }
            // Species with another radiation loss model
            else if (radiation_model=="Landau-Lifshitz"
                 ||  radiation_model=="corrected-Landau-Lifshitz"
                 ||  radiation_model=="Niel")
            {
                 thisSpecies->particles->isQuantumParameter = true;
                 thisSpecies->radiating = true;
            }
            else if (radiation_model != "none")
            {
                ERROR("For species `" << species_type
                                      << " radiation_model must be 'none',"
                                      << " 'Landau-Lifshitz',"
                                      << " 'corrected-Landau-Lifshitz',"
                                      << " 'Niel' or 'Monte-Carlo'");
            }

            thisSpecies->radiation_model = radiation_model;

            if (radiation_model != "none")
            MESSAGE(2,"> Radiating species with model: `" << radiation_model << "`");

            // Non compatibility
            if ((dynamics_type=="borisnr")
            && (radiation_model=="Monte-Carlo"
            || radiation_model=="Landau-Lifshitz"
            || radiation_model=="corrected-Landau-Lifshitz"
            || radiation_model=="Niel"))
            {
                ERROR("For species `" << species_type
                                       << "` radiation_model `"
                                       << radiation_model
                                       << "` is not compatible with dynamics_type "
                                       << dynamics_type);

            }

        }
        // Photon species
        else if (mass == 0)
        {
            thisSpecies = new SpeciesNorm(params, patch);
            // Photon can not radiate
            thisSpecies->radiation_model = "none";
            thisSpecies->dynamics_type = "norm";

            MESSAGE(2,"> " <<species_type <<" is a photon species (mass==0).");
            MESSAGE(2,"> Radiation model set to none.");
            MESSAGE(2,"> Dynamic model set to norm.");
        }

        thisSpecies->species_type = species_type;
        thisSpecies->mass = mass;
        thisSpecies->speciesNumber = ispec;

        // Extract various parameters from the namelist

        // Radiation photon
        if (mass > 0.)
        {
            if (thisSpecies->radiation_model == "Monte-Carlo")
            {
                PyTools::extract("radiation_photons", thisSpecies->radiation_photons, "Species",ispec);
                if (thisSpecies->radiation_photons != "none")
                {
                    MESSAGE(2,"> radiation_photons set to the species `" << thisSpecies->radiation_photons << "`");
                }
                PyTools::extract("radiation_photon_sampling",
                                 thisSpecies->radiation_photon_sampling, "Species",ispec);
                if (thisSpecies->radiation_photon_sampling < 1)
                {
                    ERROR("For species '" << species_type << "' radiation_photon_sampling should be > 1");
                }
            }
        }

        // Multiphoton Breit-Wheeler
        if (mass == 0)
        {
            PyTools::extract("multiphoton_Breit_Wheeler", thisSpecies->multiphoton_Breit_Wheeler, "Species",ispec);
            // If the first species is not empty
            if (thisSpecies->multiphoton_Breit_Wheeler[0] == "" ||
               (thisSpecies->multiphoton_Breit_Wheeler[1] == "" && thisSpecies->multiphoton_Breit_Wheeler[0] != "none"))
            {
                ERROR("For species '" << species_type << "' multiphoton_Breit_Wheeler can not be empty, select electron and positron species.");
            }
            // Else If the first species is not none which means no process
            else if (thisSpecies->multiphoton_Breit_Wheeler[0] != "none")
            {
                if (thisSpecies->multiphoton_Breit_Wheeler[1] != "none")
                {
                    // Activation of the additional variables
                    thisSpecies->particles->isQuantumParameter = true;
                    thisSpecies->particles->isMonteCarlo = true;

                    thisSpecies->mBW_pair_creation_sampling.resize(2);
                    PyTools::extract("multiphoton_Breit_Wheeler_sampling",
                                     thisSpecies->mBW_pair_creation_sampling, "Species",ispec);

                    MESSAGE(2,"> Decay into pair via the multiphoton Breit-Wheeler activated");
                    MESSAGE(2,"> Generated electrons and positrons go to species: "
                    << thisSpecies->multiphoton_Breit_Wheeler[0]
                    << " & " << thisSpecies->multiphoton_Breit_Wheeler[1]);
                }
                else
                {
                    thisSpecies->multiphoton_Breit_Wheeler[0] = "none";
                    WARNING("For species '" << species_type << "' Positron species is none, Multiphoton Breit-Wheeler disabled");
                }
            }
        }

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
        // Matter particles
        if (thisSpecies->mass > 0) {
            if (   (thisSpecies->initMomentum_type!="cold")
                && (thisSpecies->initMomentum_type!="maxwell-juettner")
                && (thisSpecies->initMomentum_type!="rectangular") ) {
                    ERROR("For particle species '" << species_type
                                                << "' unknown initMomentum_type: "
                                                <<thisSpecies->initMomentum_type);
            }
        }
        // Photons
        else if (thisSpecies->mass == 0)
        {
            if (   (thisSpecies->initMomentum_type!="cold")
                && (thisSpecies->initMomentum_type!="rectangular") ) {
                    ERROR("For photon species '" << species_type
                                                << "' unknown initMomentum_type: "
                                                <<thisSpecies->initMomentum_type);
            }
        }

        PyTools::extract("c_part_max",thisSpecies->c_part_max,"Species",ispec);

        PyTools::extract("time_frozen",thisSpecies->time_frozen ,"Species",ispec);
        if (thisSpecies->time_frozen > 0 && thisSpecies->initMomentum_type!="cold") {
            if ( patch->isMaster() ) WARNING("For species '" << species_type << "' possible conflict between time-frozen & not cold initialization");
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
        // Matter particles
        if (thisSpecies->mass > 0) {
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

        }
        // Photons
        else if (thisSpecies->mass == 0)
        {
            if ( (thisSpecies->bc_part_type_xmin=="thermalize") ||
                 (thisSpecies->bc_part_type_xmax=="thermalize") ||
                 (thisSpecies->bc_part_type_ymin=="thermalize") ||
                 (thisSpecies->bc_part_type_ymax=="thermalize") ||
                 (thisSpecies->bc_part_type_zmin=="thermalize") ||
                 (thisSpecies->bc_part_type_zmax=="thermalize"))
            {
                ERROR("For photon species '" << species_type
                       << "' Thermalizing BCs are not available.");
            }
            if ( (thisSpecies->bc_part_type_xmin=="stop") ||
                 (thisSpecies->bc_part_type_xmax=="stop") ||
                 (thisSpecies->bc_part_type_ymin=="stop") ||
                 (thisSpecies->bc_part_type_ymax=="stop") ||
                 (thisSpecies->bc_part_type_zmin=="stop") ||
                 (thisSpecies->bc_part_type_zmax=="stop"))
            {
                ERROR("For photon species '" << species_type
                       << "' stop BCs are not physical.");
            }
        }



        // Manage the ionization parameters
        if (thisSpecies->mass > 0)
        {
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
        }

        // Species geometry
        // ----------------

        // Density
        bool ok1, ok2;
        PyObject *profile1, *profile2, *profile3;
        // Matter particles
        if (thisSpecies->mass > 0)
        {
            ok1 = PyTools::extract_pyProfile("nb_density"    , profile1, "Species", ispec);
            ok2 = PyTools::extract_pyProfile("charge_density", profile1, "Species", ispec);
            if(  ok1 &&  ok2 ) ERROR("For species '" << species_type << "', cannot define both `nb_density` and `charge_density`.");
            if( !ok1 && !ok2 ) ERROR("For species '" << species_type << "', must define `nb_density` or `charge_density`.");
            if( ok1 ) thisSpecies->densityProfileType = "nb";
            if( ok2 ) thisSpecies->densityProfileType = "charge";
        }
        // Photons
        else if (thisSpecies->mass == 0)
        {
            ok1 = PyTools::extract_pyProfile("nb_density"    , profile1, "Species", ispec);
            ok2 = PyTools::extract_pyProfile("charge_density", profile1, "Species", ispec);
            if( !ok1 && ok2 ) ERROR("For photon species '" << species_type
                         << "', must define `nb_density`, `charge_density has no meaning for photons`.");
            if (!ok1) ERROR("For photon species '" << species_type
                         << "', must define `nb_density`.");
            if( ok1 ) thisSpecies->densityProfileType = "nb";
        }

        thisSpecies->densityProfile = new Profile(profile1, params.nDim_particle, thisSpecies->densityProfileType+"_density "+species_type, true);

        // Number of particles per cell
        if( !PyTools::extract_pyProfile("n_part_per_cell", profile1, "Species", ispec))
            ERROR("For species '" << species_type << "', n_part_per_cell not found or not understood");
        thisSpecies->ppcProfile = new Profile(profile1, params.nDim_particle, "n_part_per_cell "+species_type, true);

        // Charge
        if( !PyTools::extract_pyProfile("charge", profile1, "Species", ispec))
            ERROR("For species '" << species_type << "', charge not found or not understood");
        thisSpecies->chargeProfile = new Profile(profile1, params.nDim_particle, "charge "+species_type, true);

        // Mean velocity
        PyTools::extract3Profiles("mean_velocity", ispec, profile1, profile2, profile3);
        thisSpecies->velocityProfile[0] = new Profile(profile1, params.nDim_particle, "mean_velocity[0] "+species_type, true);
        thisSpecies->velocityProfile[1] = new Profile(profile2, params.nDim_particle, "mean_velocity[1] "+species_type, true);
        thisSpecies->velocityProfile[2] = new Profile(profile3, params.nDim_particle, "mean_velocity[2] "+species_type, true);

        // Temperature
        PyTools::extract3Profiles("temperature", ispec, profile1, profile2, profile3);
        thisSpecies->temperatureProfile[0] = new Profile(profile1, params.nDim_particle, "temperature[0] "+species_type, true);
        thisSpecies->temperatureProfile[1] = new Profile(profile2, params.nDim_particle, "temperature[1] "+species_type, true);
        thisSpecies->temperatureProfile[2] = new Profile(profile3, params.nDim_particle, "temperature[2] "+species_type, true);


        // CALCULATE USEFUL VALUES

        /*        double gamma=1.+thisSpecies->thermT[0]/thisSpecies->mass;

                  for (unsigned int i=0; i<3; i++) {
                  thisSpecies->thermalVelocity[i] = sqrt( 1.-1./gamma*gamma );
                  thisSpecies->thermalMomentum[i] = gamma*thisSpecies->thermalVelocity[i];
                  }

                  double gamma=1.+thisSpecies->thermT[0]/thisSpecies->mass;
        */

        // Matter particles
        if (thisSpecies->mass > 0) {
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
        {
            thisSpecies->particles->initialize( 0, params.nDim_particle );
        }

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
        || species->dynamics_type=="borisnr")
        {
            // Boris, Vay or Higuera-Cary
            newSpecies = new SpeciesNorm(params, patch);
        }

        // Copy members
        newSpecies->species_type          = species->species_type;
        newSpecies->dynamics_type         = species->dynamics_type;
        newSpecies->radiation_model       = species->radiation_model;
        newSpecies->radiation_photons     = species->radiation_photons;
        newSpecies->radiation_photon_sampling = species->radiation_photon_sampling;
        newSpecies->photon_species        = species->photon_species;
        newSpecies->speciesNumber         = species->speciesNumber;
        newSpecies->initPosition_type     = species->initPosition_type;
        newSpecies->initMomentum_type     = species->initMomentum_type;
        newSpecies->c_part_max            = species->c_part_max;
        newSpecies->mass                  = species->mass;
        newSpecies->time_frozen           = species->time_frozen;
        newSpecies->radiating             = species->radiating;
        newSpecies->bc_part_type_xmin     = species->bc_part_type_xmin;
        newSpecies->bc_part_type_xmax     = species->bc_part_type_xmax;
        newSpecies->bc_part_type_ymin     = species->bc_part_type_ymin;
        newSpecies->bc_part_type_ymax     = species->bc_part_type_ymax;
        newSpecies->bc_part_type_zmin     = species->bc_part_type_zmin;
        newSpecies->bc_part_type_zmax     = species->bc_part_type_zmax;
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
        if (newSpecies->mass==0) {
            newSpecies->multiphoton_Breit_Wheeler[0]   = species->multiphoton_Breit_Wheeler[0];
            newSpecies->multiphoton_Breit_Wheeler[1]   = species->multiphoton_Breit_Wheeler[1];
            newSpecies->mBW_pair_creation_sampling[0] = species->mBW_pair_creation_sampling[0];
            newSpecies->mBW_pair_creation_sampling[1] = species->mBW_pair_creation_sampling[1];
        }

        newSpecies->particles->isTest              = species->particles->isTest;
        newSpecies->particles->tracked             = species->particles->tracked;
        newSpecies->particles->isQuantumParameter  = species->particles->isQuantumParameter;
        newSpecies->particles->isMonteCarlo        = species->particles->isMonteCarlo;

        // \todo : NOT SURE HOW THIS BEHAVES WITH RESTART
        if ( (!params.restart) && (with_particles) ) {
            newSpecies->createParticles(params.n_space, params, patch, 0 );
        }
        else
        {
            newSpecies->particles->initialize( 0, (*species->particles) );
        }

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
            if (retSpecies[ispec1]->electron_species_index==-1) {
                ERROR("For species '"<<retSpecies[ispec1]->species_type<<"' ionization_electrons named " << retSpecies[ispec1]->ionization_electrons << " could not be found");
            }
        }

        // Loop species to find the photon species for radiation species
        for (unsigned int ispec1 = 0; ispec1<retSpecies.size(); ispec1++) {
            if( ! retSpecies[ispec1]->Radiate )
            {
                continue;
            }

            // No emission of discrete photon, only scalar diagnostics are updated
            if( retSpecies[ispec1]->radiation_photons == "none")
            {
                retSpecies[ispec1]->photon_species_index = -1;
                retSpecies[ispec1]->photon_species = NULL;
            }
            // Else, there will be emission of macro-photons.
            else
            {
                for (unsigned int ispec2 = 0; ispec2<retSpecies.size(); ispec2++) {
                    if( retSpecies[ispec1]->radiation_photons == retSpecies[ispec2]->species_type) {
                        if( ispec1==ispec2 )
                            ERROR("For species '"<<retSpecies[ispec1]->species_type<<"' radiation_photons must be a distinct photon species");
                        if (retSpecies[ispec2]->mass!=0)
                        {
                            ERROR("For species '"<<retSpecies[ispec1]->species_type<<"' radiation_photons must be a photon species with mass==0");
                        }
                        retSpecies[ispec1]->photon_species_index = ispec2;
                        retSpecies[ispec1]->photon_species = retSpecies[ispec2];
                        retSpecies[ispec1]->Radiate->new_photons.tracked = retSpecies[ispec1]->photon_species->particles->tracked;
                        retSpecies[ispec1]->Radiate->new_photons.isQuantumParameter = retSpecies[ispec1]->photon_species->particles->isQuantumParameter;
                        retSpecies[ispec1]->Radiate->new_photons.isMonteCarlo = retSpecies[ispec1]->photon_species->particles->isMonteCarlo;
                        retSpecies[ispec1]->Radiate->new_photons.initialize(0,
                                                                            params.nDim_particle );
                        //retSpecies[ispec1]->Radiate->new_photons.initialize(retSpecies[ispec1]->getNbrOfParticles(),
                        //                                                    params.nDim_particle );
                        retSpecies[ispec2]->particles->reserve(retSpecies[ispec1]->getNbrOfParticles(),
                                                               retSpecies[ispec2]->particles->dimension() );
                    }
                }
            }
        }

        // Loop species to find the electron and positron
        // species for multiphoton Breit-wheeler
        for (unsigned int ispec1 = 0; ispec1<retSpecies.size(); ispec1++) {
            if(!retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process)
            {
                continue;
            }
            else
            {
                for (unsigned int ispec2 = 0; ispec2<retSpecies.size(); ispec2++)
                {
                    for (int k=0;k<2;k++)
                    {
                        if( retSpecies[ispec1]->multiphoton_Breit_Wheeler[k] == retSpecies[ispec2]->species_type)
                        {
                            if( ispec1==ispec2 )
                            {
                                ERROR("For species '" << retSpecies[ispec1]->species_type
                                                      << "' pair species must be a distinct particle species");
                            }
                            if (retSpecies[ispec2]->mass != 1)
                            {
                                ERROR("For species '"<<retSpecies[ispec1]->species_type<<"' pair species must be an electron and positron species");
                            }
                            retSpecies[ispec1]->mBW_pair_species_index[k] = ispec2;
                            retSpecies[ispec1]->mBW_pair_species[k] = retSpecies[ispec2];
                            retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].tracked = retSpecies[ispec1]->mBW_pair_species[k]->particles->tracked;
                            retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].isQuantumParameter = retSpecies[ispec1]->mBW_pair_species[k]->particles->isQuantumParameter;
                            retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].isMonteCarlo = retSpecies[ispec1]->mBW_pair_species[k]->particles->isMonteCarlo;
                            retSpecies[ispec1]->Multiphoton_Breit_Wheeler_process->new_pair[k].initialize(0,
                                                                                params.nDim_particle );
                            retSpecies[ispec2]->particles->reserve(retSpecies[ispec1]->getNbrOfParticles(),
                                                                   retSpecies[ispec2]->particles->dimension() );
                        }
                    }
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

        // Ionization
        for (unsigned int i=0; i<retSpecies.size(); i++) {
            if (retSpecies[i]->Ionize) {
                retSpecies[i]->electron_species_index = vecSpecies[i]->electron_species_index;
                retSpecies[i]->electron_species = retSpecies[retSpecies[i]->electron_species_index];
                retSpecies[i]->Ionize->new_electrons.tracked = retSpecies[i]->electron_species->particles->tracked;
                retSpecies[i]->Ionize->new_electrons.initialize(0, params.nDim_particle );
            }
        }

        // Synchortron-like radiation
        for (unsigned int i=0; i<retSpecies.size(); i++) {
            if (retSpecies[i]->Radiate) {
                retSpecies[i]->radiation_photons = vecSpecies[i]->radiation_photons;
                retSpecies[i]->photon_species_index = vecSpecies[i]->photon_species_index;
                if (vecSpecies[i]->photon_species)
                {
                    retSpecies[i]->photon_species = retSpecies[retSpecies[i]->photon_species_index];
                    retSpecies[i]->Radiate->new_photons.tracked = retSpecies[i]->photon_species->particles->tracked;
                    retSpecies[i]->Radiate->new_photons.isQuantumParameter = retSpecies[i]->photon_species->particles->isQuantumParameter;
                    retSpecies[i]->Radiate->new_photons.isMonteCarlo = retSpecies[i]->photon_species->particles->isMonteCarlo;
                    //retSpecies[i]->Radiate->new_photons.initialize(retSpecies[i]->getNbrOfParticles(),
                    //                                               params.nDim_particle );
                    retSpecies[i]->Radiate->new_photons.initialize(0,
                                                                  params.nDim_particle );
                }
                else
                {
                    retSpecies[i]->photon_species = NULL;
                }
            }
        }

        // multiphoton Breit-Wheeler
        for (unsigned int i=0; i<retSpecies.size(); i++)
        {
            if (retSpecies[i]->Multiphoton_Breit_Wheeler_process) {
                // Loop on pairs
                for (int k=0;k<2;k++)
                {
                    retSpecies[i]->multiphoton_Breit_Wheeler[k] = vecSpecies[i]->multiphoton_Breit_Wheeler[k];
                    retSpecies[i]->mBW_pair_species_index[k] = vecSpecies[i]->mBW_pair_species_index[k];
                    retSpecies[i]->mBW_pair_species[k] = retSpecies[retSpecies[i]->mBW_pair_species_index[k]];
                    retSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].tracked = retSpecies[i]->mBW_pair_species[k]->particles->tracked;
                    retSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].isQuantumParameter = retSpecies[i]->mBW_pair_species[k]->particles->isQuantumParameter;
                    retSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].isMonteCarlo = retSpecies[i]->mBW_pair_species[k]->particles->isMonteCarlo;
                    retSpecies[i]->Multiphoton_Breit_Wheeler_process->new_pair[k].initialize(
                                        0,params.nDim_particle );
                }
            }
            else
            {
                retSpecies[i]->mBW_pair_species[0] = NULL;
                retSpecies[i]->mBW_pair_species[1] = NULL;
            }
        }

        return retSpecies;
    }

};

#endif
