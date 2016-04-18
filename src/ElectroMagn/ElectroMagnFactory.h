#ifndef ELECTROMAGNFACTORY_H
#define ELECTROMAGNFACTORY_H

#include <sstream>
#include "ElectroMagn.h"
#include "ElectroMagn1D.h"
#include "ElectroMagn2D.h"
#include "ElectroMagnBC.h"

#include "Patch.h"
#include "Params.h"
#include "Laser.h"
#include "Tools.h"

class ElectroMagnFactory {
public:
    static ElectroMagn* create(Params& params, std::vector<Species*>& vecSpecies,  Patch* patch) {
        ElectroMagn* EMfields = NULL;
        if ( params.geometry == "1d3v" ) {
            EMfields = new ElectroMagn1D(params, vecSpecies, patch);
        }
        else if ( params.geometry == "2d3v" ) {
            EMfields = new ElectroMagn2D(params, vecSpecies, patch);
        }
        else {
            ERROR( "Unknown geometry : " << params.geometry );
        }
        
        // -----------------
        // Lasers properties
        // -----------------
        if( patch->isMaster() ) MESSAGE(1, "Laser parameters :");
        int nlaser = PyTools::nComponents("Laser");
        for (int ilaser = 0; ilaser < nlaser; ilaser++) {
            Laser * laser = new Laser(params, ilaser, patch);
            if     ( laser->boxSide == "west" && EMfields->emBoundCond[0])
                EMfields->emBoundCond[0]->vecLaser.push_back( laser );
            else if( laser->boxSide == "east" && EMfields->emBoundCond[1])
                EMfields->emBoundCond[1]->vecLaser.push_back( laser );
        }
        
        // -----------------
        // ExtFields properties
        // -----------------
        unsigned int numExtFields=PyTools::nComponents("ExtField");
        for (unsigned int n_extfield = 0; n_extfield < numExtFields; n_extfield++) {
            MESSAGE("ExtField " << n_extfield);
            ExtField extField;
            PyObject * profile;
            std::ostringstream name;
            if( !PyTools::extract("field",extField.fields,"ExtField",n_extfield))
                ERROR("ExtField #"<<n_extfield<<": parameter 'field' not provided'");
            
            // Now import the profile
            name.str("");
            name << "ExtField[" << n_extfield <<"].profile";
            if (!PyTools::extract_pyProfile("profile",profile,"ExtField",n_extfield))
                ERROR(" ExtField #"<<n_extfield<<": parameter 'profile' not understood");
            extField.profile = new Profile(profile, params.nDim_field, name.str());
            
            EMfields->extFields.push_back(extField);
        }
        
        
        // -----------------
        // Antenna properties
        // -----------------
        unsigned int numAntenna=PyTools::nComponents("Antenna");
        for (unsigned int n_antenna = 0; n_antenna < numAntenna; n_antenna++) {
            Antenna antenna;
            PyObject * profile;
            std::ostringstream name;
            antenna.field = NULL;
            if( !PyTools::extract("field",antenna.fieldName,"Antenna",n_antenna))
                ERROR("Antenna #"<<n_antenna<<": parameter 'field' not provided'");
            if (antenna.fieldName != "Jx" && antenna.fieldName != "Jy" && antenna.fieldName != "Jz")
                ERROR("Antenna #"<<n_antenna<<": parameter 'field' must be one of Jx, Jy, Jz");
            
            // Extract the space profile
            name.str("");
            name << "Antenna[" << n_antenna <<"].space_profile";
            if (!PyTools::extract_pyProfile("space_profile",profile,"Antenna",n_antenna))
                ERROR(" Antenna #"<<n_antenna<<": parameter 'space_profile' not understood");
            antenna.space_profile = new Profile(profile, params.nDim_field, name.str());
            
            // Extract the time profile
            name.str("");
            name << "Antenna[" << n_antenna <<"].time_profile";
            if (!PyTools::extract_pyProfile("time_profile" ,profile,"Antenna",n_antenna))
                ERROR(" Antenna #"<<n_antenna<<": parameter 'time_profile' not understood");
            antenna.time_profile =  new Profile(profile, 1, name.str());
            
            EMfields->antennas.push_back(antenna);
        }
        
        
        EMfields->finishInitialization(vecSpecies.size(), patch);
        
        // Some output
        std::stringstream ss;
        for (std::vector<Field*>::iterator iterField=EMfields->allFields.begin(); iterField!=EMfields->allFields.end(); iterField++) {
            ss << (*iterField)->name << " ";
        }
        if (patch->isMaster()) {
            MESSAGE(1,"EM fields dump      :");
            MESSAGE(2, ss.str() );
        }
        ss.str("");
        for (std::vector<Field*>::iterator iterField=EMfields->allFields_avg.begin(); iterField!=EMfields->allFields_avg.end(); iterField++) {
            ss << (*iterField)->name << " ";
        }
        if (patch->isMaster()) {
            MESSAGE(1,"EM avg. fields dump :");
            MESSAGE(2, ss.str() );
        }
        
        return EMfields;
    }
    
    
    static ElectroMagn* clone(ElectroMagn* EMfields, Params& params, std::vector<Species*>& vecSpecies,  Patch* patch)
    {
        ElectroMagn* newEMfields = NULL;
        if ( params.geometry == "1d3v" ) {
            newEMfields = new ElectroMagn1D(params, vecSpecies, patch);
        } else if ( params.geometry == "2d3v" ) {
            newEMfields = new ElectroMagn2D(params, vecSpecies, patch);
        }
        
        // -----------------
        // Clone Lasers properties
        // -----------------
        int nlaser;
        newEMfields->emBoundCond[0]->vecLaser.resize(0);
        nlaser = EMfields->emBoundCond[0]->vecLaser.size();
        for (int ilaser = 0; ilaser < nlaser; ilaser++) {
            newEMfields->emBoundCond[0]->vecLaser.push_back(
                new Laser(EMfields->emBoundCond[0]->vecLaser[ilaser])
            );
        }
        newEMfields->emBoundCond[1]->vecLaser.resize(0);
        nlaser = EMfields->emBoundCond[1]->vecLaser.size();
        for (int ilaser = 0; ilaser < nlaser; ilaser++) {
            newEMfields->emBoundCond[1]->vecLaser.push_back(
                new Laser(EMfields->emBoundCond[1]->vecLaser[ilaser])
            );
        }
        
        // -----------------
        // Clone ExtFields properties
        // -----------------
        for (unsigned int n_extfield = 0; n_extfield < EMfields->extFields.size(); n_extfield++) {
            ExtField extField;
            extField.fields  = EMfields->extFields[n_extfield].fields;
            extField.profile = EMfields->extFields[n_extfield].profile;
            newEMfields->extFields.push_back(extField);
        }
        
        // -----------------
        // Clone Antenna properties
        // -----------------
        for (unsigned int n_antenna = 0; n_antenna < EMfields->antennas.size(); n_antenna++) {
            Antenna antenna;
            antenna.fieldName     = EMfields->antennas[n_antenna].fieldName    ;
            antenna.field         = EMfields->antennas[n_antenna].field        ;
            antenna.space_profile = EMfields->antennas[n_antenna].space_profile;
            antenna.time_profile  = EMfields->antennas[n_antenna].time_profile ;
            newEMfields->antennas.push_back(antenna);
        }
        
        
        newEMfields->finishInitialization(vecSpecies.size(), patch);
        
        return newEMfields;
    }

};

#endif

