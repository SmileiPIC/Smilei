#ifndef SPECIESFACTORY_H
#define SPECIESFACTORY_H

#include "Species.h"
#include "Species_norm.h"
#include "Species_rrll.h"

#include "Params.h"
#include "SmileiMPI.h"

#include "Tools.h"

class SpeciesFactory {
public:

    static std::vector<Species*> createVector(Params& params, SmileiMPI* smpi) {
        // this will be returned
        std::vector<Species*> retSpecies;
        
        // at some point we will find it... it is needed for ionization
        Species *electron_species=NULL;

        // read from python namelist
        bool ok;
        for (unsigned int ispec = 0; ispec < (unsigned int) PyTools::nComponents("Species"); ispec++) {
            SpeciesStructure spec_struct;
            
            spec_struct.speciesNumber=ispec;
            
            PyTools::extract("species_type",spec_struct.species_type,"Species",ispec);
            if(spec_struct.species_type.empty()) {
                ERROR("For species #" << ispec << " empty type");
            }
            PyTools::extract("initPosition_type",spec_struct.initPosition_type ,"Species",ispec);
            if (spec_struct.initPosition_type.empty()) {
                ERROR("For species #" << ispec << " empty initPosition_type");
            } else if ( (spec_struct.initPosition_type!="regular")&&(spec_struct.initPosition_type!="random") ) {
                ERROR("For species #" << ispec << " bad definition of initPosition_type " << spec_struct.initPosition_type);
            }
            
            PyTools::extract("initMomentum_type",spec_struct.initMomentum_type ,"Species",ispec);
            if ( (spec_struct.initMomentum_type=="mj") || (spec_struct.initMomentum_type=="maxj") ) {
                spec_struct.initMomentum_type="maxwell-juettner";
            }
            if (   (spec_struct.initMomentum_type!="cold")
                && (spec_struct.initMomentum_type!="maxwell-juettner")
                && (spec_struct.initMomentum_type!="rectangular") ) {
                ERROR("For species #" << ispec << " bad definition of initMomentum_type");
            }
            
            spec_struct.c_part_max = 1.0;// default value
            PyTools::extract("c_part_max",spec_struct.c_part_max,"Species",ispec);
            
            if( !PyTools::extract("mass",spec_struct.mass ,"Species",ispec) ) {
                ERROR("For species #" << ispec << ", mass not defined.");
            }
            
            spec_struct.dynamics_type = "norm"; // default value
            if (!PyTools::extract("dynamics_type",spec_struct.dynamics_type ,"Species",ispec) )
                WARNING("For species #" << ispec << ", dynamics_type not defined: assumed = 'norm'.");
            if (spec_struct.dynamics_type!="norm"){
                ERROR("dynamics_type different than norm not yet implemented");
            }
            
            spec_struct.time_frozen = 0.0; // default value
            PyTools::extract("time_frozen",spec_struct.time_frozen ,"Species",ispec);
            if (spec_struct.time_frozen > 0 && \
                spec_struct.initMomentum_type!="cold") {
                WARNING("For species #" << ispec << " possible conflict between time-frozen & not cold initialization");
            }
            
            spec_struct.radiating = false; // default value
            PyTools::extract("radiating",spec_struct.radiating ,"Species",ispec);
            if (spec_struct.dynamics_type=="rrll" && (!spec_struct.radiating)) {
                WARNING("For species #" << ispec << ", dynamics_type='rrll' forcing radiating=True");
                spec_struct.radiating=true;
            }
            
            if (!PyTools::extract("bc_part_type_west",spec_struct.bc_part_type_west,"Species",ispec) )
                ERROR("For species #" << ispec << ", bc_part_type_west not defined");
            if (!PyTools::extract("bc_part_type_east",spec_struct.bc_part_type_east,"Species",ispec) )
                ERROR("For species #" << ispec << ", bc_part_type_east not defined");
            
            if (params.nDim_particle>1) {
                if (!PyTools::extract("bc_part_type_south",spec_struct.bc_part_type_south,"Species",ispec) )
                    ERROR("For species #" << ispec << ", bc_part_type_south not defined");
                if (!PyTools::extract("bc_part_type_north",spec_struct.bc_part_type_north,"Species",ispec) )
                    ERROR("For species #" << ispec << ", bc_part_type_north not defined");
            }
            
            // for thermalizing BCs on particles check if thermT is correctly defined
            bool thermTisDefined=false;
            if ( (spec_struct.bc_part_type_west=="thermalize") || (spec_struct.bc_part_type_east=="thermalize") ){
                thermTisDefined=PyTools::extract("thermT",spec_struct.thermT,"Species",ispec);
                if (!thermTisDefined) ERROR("thermT needs to be defined for species " <<ispec<< " due to x-BC thermalize");
            }
            if ( (params.nDim_particle==2) && (!thermTisDefined) &&
                (spec_struct.bc_part_type_south=="thermalize" || spec_struct.bc_part_type_north=="thermalize") ) {
                thermTisDefined=PyTools::extract("thermT",spec_struct.thermT,"Species",ispec);
                if (!thermTisDefined) ERROR("thermT needs to be defined for species " <<ispec<< " due to y-BC thermalize");
            }
            if (thermTisDefined) {
                if (spec_struct.thermT.size()==1) {
                    spec_struct.thermT.resize(3);
                    for (unsigned int i=1; i<3;i++)
                        spec_struct.thermT[i]=spec_struct.thermT[0];
                }
            } else {
                spec_struct.thermT.resize(3);
                for (unsigned int i=0; i<3;i++)
                    spec_struct.thermT[i]=0.0;
            }
            
            
            spec_struct.ionization_model = "none"; // default value
            PyTools::extract("ionization_model", spec_struct.ionization_model, "Species",ispec);
            
            ok = PyTools::extract("atomic_number", spec_struct.atomic_number, "Species",ispec);
            if( !ok && spec_struct.ionization_model!="none" ) {
                ERROR("For species #" << ispec << ", `atomic_number` not found => required for the ionization model .");
            }
            
            spec_struct.isTest = false; // default value
            PyTools::extract("isTest",spec_struct.isTest ,"Species",ispec);
            if (spec_struct.ionization_model!="none" && (!spec_struct.isTest)) {
                ERROR("For species #" << ispec << ", disabled for now : test & ionized");
            }
            // Define the number of timesteps for dumping test particles
            spec_struct.test_dump_every = 1;
            if (PyTools::extract("dump_every",spec_struct.test_dump_every ,"Species",ispec)) {
                if (spec_struct.test_dump_every>1 && !spec_struct.isTest)
                    WARNING("For species #" << ispec << ", dump_every discarded because not test particles");
            }
            
            
            // Species geometry
            // ----------------
            
            // Density
            bool ok1, ok2;
            ok1 = PyTools::extract_pyProfile("nb_density"    , spec_struct.dens_profile, "Species", ispec);
            ok2 = PyTools::extract_pyProfile("charge_density", spec_struct.dens_profile, "Species", ispec);
            
            if(  ok1 &&  ok2 ) ERROR("For species #" << ispec << ", cannot define both `nb_density` and `charge_density`.");
            if( !ok1 && !ok2 ) ERROR("For species #" << ispec << ", must define `nb_density` or `charge_density`.");
            if( ok1 ) spec_struct.density_type = "nb";
            if( ok2 ) spec_struct.density_type = "charge";
            
            // Number of particles per cell
            if( !PyTools::extract_pyProfile("n_part_per_cell"    , spec_struct.ppc_profile, "Species", ispec))
                ERROR("For species #" << ispec << ", n_part_per_cell not found or not understood");
            
            // Charge
            if( !PyTools::extract_pyProfile("charge"    , spec_struct.charge_profile, "Species", ispec))
                ERROR("For species #" << ispec << ", charge not found or not understood");
            
            // Mean velocity
            PyTools::extract3Profiles("mean_velocity", ispec, spec_struct.mvel_x_profile, spec_struct.mvel_y_profile, spec_struct.mvel_z_profile);
            
            // Temperature
            PyTools::extract3Profiles("temperature", ispec, spec_struct.temp_x_profile, spec_struct.temp_y_profile, spec_struct.temp_z_profile);
            
            
            // CALCULATE USEFUL VALUES
            
            // define thermal velocity as \sqrt{T/m}
            spec_struct.thermalVelocity.resize(3);
            spec_struct.thermalMomentum.resize(3);
            for (unsigned int i=0; i<3; i++) {
                spec_struct.thermalVelocity[i]=sqrt(2.0*spec_struct.thermT[i]/spec_struct.mass);
                spec_struct.thermalMomentum[i]=spec_struct.mass * spec_struct.thermalVelocity[i];
            }

            // create species
            Species *thisSpecies=NULL;
            
            if (spec_struct.dynamics_type=="norm") {
                // Species with Boris dynamics
                thisSpecies = new Species_norm(params, spec_struct, smpi);
            } else if (spec_struct.dynamics_type=="rrll") {
                // Species with Boris dynamics + Radiation Back-Reaction (using the Landau-Lifshitz formula)
                thisSpecies = new Species_rrll(params, spec_struct, smpi);
            } else {
                ERROR("Species " << spec_struct.species_type << " dynamics_type must be either norm or rrll")
            }
            
            if (spec_struct.isTest) {
                int locNbrParticles = thisSpecies->getNbrOfParticles();
                int* allNbrParticles = new int[smpi->smilei_sz];
                MPI_Gather( &locNbrParticles, 1, MPI_INTEGER, allNbrParticles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
                int nParticles(0);
                if (smpi->isMaster()) {
                    nParticles =  allNbrParticles[0];
                    for (int irk=1 ; irk<smpi->getSize() ; irk++){
                        allNbrParticles[irk] += nParticles;
                        nParticles = allNbrParticles[irk];
                    }
                    for (int irk=smpi->getSize()-1 ; irk>0 ; irk--){
                        allNbrParticles[irk] = allNbrParticles[irk-1];
                    }
                    allNbrParticles[0] = 0;
                    
                }
                int offset(0);
                MPI_Scatter(allNbrParticles, 1 , MPI_INTEGER, &offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
                thisSpecies->particles.addIdOffsets(offset);
            }
            
            retSpecies.push_back(thisSpecies);
            
            unsigned int nPart = thisSpecies->getNbrOfParticles();
            MPI_Reduce(smpi->isMaster()?MPI_IN_PLACE:&nPart, &nPart, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
            MESSAGE(1,"Species " << ispec << " (" << spec_struct.species_type << ") created with " << nPart << " particles" );
        
            if (thisSpecies->sparams.species_type=="electron") {
                electron_species=thisSpecies;
            }
        }
        
        // we cycle again to fix electron species for ionizable species
        for (unsigned int i=0; i<retSpecies.size(); i++) {
            if (retSpecies[i]->Ionize)  {
                if (electron_species) {
                    retSpecies[i]->electron_species=electron_species;
                    PMESSAGE(2,smpi->getRank(),"Added electron species to species " << retSpecies[i]->sparams.species_type);
                } else {
                    ERROR("Ionization needs a species called \"electron\" to be defined");
                }
            }
        }
        
        
        
        // Plasma related parameters
        // -------------------------
        MESSAGE("Plasma related parameters");
        MESSAGE(1,"n_species       : " << retSpecies.size());
        for ( unsigned int i=0 ; i<retSpecies.size() ; i++ ) {
            MESSAGE(1,"            type : "<< retSpecies[i]->sparams.species_type);
        }
        
        return retSpecies;
    }

};

#endif
