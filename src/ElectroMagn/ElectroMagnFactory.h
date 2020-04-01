#ifndef ELECTROMAGNFACTORY_H
#define ELECTROMAGNFACTORY_H

#include <sstream>
#include "ElectroMagn.h"
#include "ElectroMagn1D.h"
#include "ElectroMagn2D.h"
#include "ElectroMagn3D.h"
#include "ElectroMagnAM.h"
#include "ElectroMagnBC.h"
#include "EnvelopeFactory.h"

#include "Patch.h"
#include "Params.h"
#include "Laser.h"
#include "Tools.h"


class ElectroMagnFactory
{
public:
    static ElectroMagn *create( Params &params, DomainDecomposition *domain_decomposition, std::vector<Species *> &vecSpecies,  Patch *patch )
    {
        ElectroMagn *EMfields = NULL;
        if( params.geometry == "1Dcartesian" ) {
            EMfields = new ElectroMagn1D( params, domain_decomposition, vecSpecies, patch );
        } else if( params.geometry == "2Dcartesian" ) {
            EMfields = new ElectroMagn2D( params, domain_decomposition, vecSpecies, patch );
        } else if( params.geometry == "3Dcartesian" ) {
            EMfields = new ElectroMagn3D( params, domain_decomposition, vecSpecies, patch );
        } else if( params.geometry == "AMcylindrical" ) {
            EMfields = new ElectroMagnAM( params, domain_decomposition, vecSpecies, patch );
        } else {
            ERROR( "Unknown geometry : " << params.geometry << "!" );
        }
        
        EMfields->finishInitialization( vecSpecies.size(), patch );
        
        // initialize the envelope if used
        if( params.Laser_Envelope_model ) { // for the moment it works only with one envelope
            EMfields->envelope = EnvelopeFactory::create( params, patch, EMfields );
        }
        
        
        // -----------------
        // Lasers properties
        // -----------------
        int nlaser = PyTools::nComponents( "Laser" );
        if( patch->isMaster() && nlaser > 0) {
            TITLE("Initializing laser parameters" );
        }
        for( int ilaser = 0; ilaser < nlaser; ilaser++ ) {
            Laser *laser = new Laser( params, ilaser, patch );
            if( laser->box_side == "xmin" && EMfields->emBoundCond[0] ) {
                if( patch->isXmin() ) {
                    laser->createFields( params, patch );
                }
                EMfields->emBoundCond[0]->vecLaser.push_back( laser );
            } else if( laser->box_side == "xmax" && EMfields->emBoundCond[1] ) {
                if( patch->isXmax() ) {
                    laser->createFields( params, patch );
                }
                EMfields->emBoundCond[1]->vecLaser.push_back( laser );
            } else {
                delete laser;
            }
        }
        
        // -----------------
        // ExtFields properties
        // -----------------
        unsigned int numExtFields=PyTools::nComponents( "ExternalField" );
        if( patch->isMaster() && numExtFields > 0) {
            TITLE("Initializing External fields" );
        }
        for( unsigned int n_extfield = 0; n_extfield < numExtFields; n_extfield++ ) {
            ExtField extField;
            PyObject *profile;
            if( !PyTools::extract( "field", extField.field, "ExternalField", n_extfield ) ) {
                ERROR( "ExternalField #"<<n_extfield<<": parameter 'field' not provided'" );
            }
            // Now import the profile
            std::ostringstream name( "" );
            name << "ExternalField[" << n_extfield <<"].profile";
            if( !PyTools::extract_pyProfile( "profile", profile, "ExternalField", n_extfield ) ) {
                ERROR( "ExternalField #"<<n_extfield<<": parameter 'profile' not understood" );
            }
            extField.profile = new Profile( profile, params.nDim_field, name.str(), true );
            // Find which index the field is in the allFields vector
            extField.index = 1000;
            for( unsigned int ifield=0; ifield<EMfields->allFields.size(); ifield++ ) {
                if( EMfields->allFields[ifield]
                        && extField.field==EMfields->allFields[ifield]->name ) {
                    extField.index = ifield;
                    break;
                }
            }
            if( extField.index > EMfields->allFields.size()-1 ) {
                ERROR( "ExternalField #"<<n_extfield<<": field "<<extField.field<<" not found" );
            }
            
            MESSAGE( 1, "External field " << extField.field << ": " << extField.profile->getInfo() );
            EMfields->extFields.push_back( extField );
        }

        // -----------------
        // PrescribedFields properties
        // -----------------
        unsigned int external_time_field_number =PyTools::nComponents( "PrescribedField" );
        if( patch->isMaster() && external_time_field_number > 0) {
            TITLE("Initializing Prescribed (external) Fields" );
        }
        for( unsigned int n_extfield = 0; n_extfield < PyTools::nComponents( "PrescribedField" ); n_extfield++ ) {
            ExtTimeField extField;
            PyObject *profile;
            std::string fieldName("");
            if( !PyTools::extract( "field", fieldName, "PrescribedField", n_extfield ) ) {
                ERROR( "PrescribedField #"<<n_extfield<<": parameter 'field' not provided'" );
            }
            // Now import the profile
            std::ostringstream name( "" );
            name << "PrescribedField[" << n_extfield <<"].profile";
            if( !PyTools::extract_pyProfile( "profile", profile, "PrescribedField", n_extfield ) ) {
                ERROR( "PrescribedField #"<<n_extfield<<": parameter 'profile' not understood" );
            }
            extField.profile = new Profile( profile, params.nDim_field+1, name.str(), true );
            // Find which index the field is in the allFields vector
            extField.index = 1000;
            for( unsigned int ifield=0; ifield<EMfields->allFields.size(); ifield++ ) {
                if( EMfields->allFields[ifield]
                        && fieldName==EMfields->allFields[ifield]->name ) {
                	
                	if (params.nDim_field == 1) {
                		extField.savedField = new Field1D(EMfields->allFields[ifield]->dims());
                	} else if (params.nDim_field == 2){
                		extField.savedField = new Field2D(EMfields->allFields[ifield]->dims());
                	} else if (params.nDim_field == 3){
                		extField.savedField = new Field3D(EMfields->allFields[ifield]->dims());
                	}
                    extField.savedField->copyFrom(EMfields->allFields[ifield]);
                    extField.savedField->name = EMfields->allFields[ifield]->name;
                    extField.index =  ifield;
                    break;
                }
            }
            if( extField.index > EMfields->allFields.size()-1 ) {
                ERROR( "PrescribedField #"<<n_extfield<<": field "<<fieldName<<" not found" );
            }
            
            MESSAGE(1, "Prescribed field " << fieldName << ": " << extField.profile->getInfo());
            EMfields->extTimeFields.push_back( extField );
        }
        
        
        // -----------------
        // Antenna properties
        // -----------------
        unsigned int antenna_number=PyTools::nComponents( "Antenna" );
        if( patch->isMaster() && antenna_number > 0) {
            TITLE("Initializing Antenna" );
        }
        for( unsigned int n_antenna = 0; n_antenna < antenna_number; n_antenna++ ) {
            Antenna antenna;
            PyObject *profile;
            std::ostringstream name;
            antenna.field = NULL;
            if( !PyTools::extract( "field", antenna.fieldName, "Antenna", n_antenna ) ) {
                ERROR( "Antenna #"<<n_antenna<<": parameter 'field' not provided'" );
            }
            if( antenna.fieldName != "Jx" && antenna.fieldName != "Jy" && antenna.fieldName != "Jz" ) {
                ERROR( "Antenna #"<<n_antenna<<": parameter 'field' must be one of Jx, Jy, Jz" );
            }
            
            // Extract the space profile
            name.str( "" );
            name << "Antenna[" << n_antenna <<"].space_profile";
            if( !PyTools::extract_pyProfile( "space_profile", profile, "Antenna", n_antenna ) ) {
                ERROR( " Antenna #"<<n_antenna<<": parameter 'space_profile' not understood" );
            }
            antenna.space_profile = new Profile( profile, params.nDim_field, name.str() );
            
            // Extract the time profile
            name.str( "" );
            name << "Antenna[" << n_antenna <<"].time_profile";
            if( !PyTools::extract_pyProfile( "time_profile", profile, "Antenna", n_antenna ) ) {
                ERROR( " Antenna #"<<n_antenna<<": parameter 'time_profile' not understood" );
            }
            antenna.time_profile =  new Profile( profile, 1, name.str() );
            
            // Find the index of the field in allFields
            antenna.index = 1000;
            for( unsigned int ifield=0; ifield<EMfields->allFields.size(); ifield++ ) {
                if( EMfields->allFields[ifield]
                        && antenna.fieldName==EMfields->allFields[ifield]->name ) {
                    antenna.index = ifield;
                    break;
                }
            }
            if( antenna.index > EMfields->allFields.size()-1 ) {
                ERROR( "Antenna #"<<n_antenna<<": field "<<antenna.fieldName<<" not found" );
            }
            
            EMfields->antennas.push_back( antenna );
        }
        
        
        return EMfields;
    }
    
    
    static ElectroMagn *clone( ElectroMagn *EMfields, Params &params, std::vector<Species *> &vecSpecies,  Patch *patch, unsigned int n_moved )
    {
        // Workaround for a Laser bug
        // count laser for later
        int nlaser_tot( 0 );
        for( int iBC=0; iBC<2; iBC++ ) { // xmax and xmin
            if( ! EMfields->emBoundCond[iBC] ) {
                continue;
            }
            nlaser_tot += EMfields->emBoundCond[iBC]->vecLaser.size();
        }
        
        ElectroMagn *newEMfields = NULL;
        if( params.geometry == "1Dcartesian" ) {
            newEMfields = new ElectroMagn1D( static_cast<ElectroMagn1D *>( EMfields ), params, patch );
        } else if( params.geometry == "2Dcartesian" ) {
            newEMfields = new ElectroMagn2D( static_cast<ElectroMagn2D *>( EMfields ), params, patch );
        } else if( params.geometry == "3Dcartesian" ) {
            newEMfields = new ElectroMagn3D( static_cast<ElectroMagn3D *>( EMfields ), params, patch );
        } else if( params.geometry == "AMcylindrical" ) {
            newEMfields = new ElectroMagnAM( static_cast<ElectroMagnAM *>( EMfields ), params, patch );
        }
        
        newEMfields->finishInitialization( vecSpecies.size(), patch );
        
        // initialize the envelope if used
        if( EMfields->envelope != NULL ) {
            newEMfields->envelope = EnvelopeFactory::clone( EMfields->envelope, patch, EMfields, params, n_moved );
        }
        
        // -----------------
        // Clone time-average fields
        // -----------------
        newEMfields->allFields_avg.resize( EMfields->allFields_avg.size() );
        for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++ ) {
            for( unsigned int ifield=0; ifield<EMfields->allFields_avg[idiag].size(); ifield++ )
                newEMfields->allFields_avg[idiag].push_back(
                    newEMfields->createField( EMfields->allFields_avg[idiag][ifield]->name )
                );
        }
        
        // -----------------
        // Clone Lasers properties
        // -----------------
        if( nlaser_tot>0 ) {
            int nlaser;
            for( int iBC=0; iBC<2; iBC++ ) { // xmax and xmin
                if( ! newEMfields->emBoundCond[iBC] ) {
                    continue;
                }
                
                newEMfields->emBoundCond[iBC]->vecLaser.resize( 0 );
                nlaser = EMfields->emBoundCond[iBC]->vecLaser.size();
                // Create lasers one by one
                for( int ilaser = 0; ilaser < nlaser; ilaser++ ) {
                    // Create laser
                    Laser *laser = new Laser( EMfields->emBoundCond[iBC]->vecLaser[ilaser], params );
                    // If patch is on border, then fill the fields arrays
                    if( ( iBC==0 && patch->isXmin() )
                            || ( iBC==1 && patch->isXmax() ) ) {
                        laser->createFields( params, patch );
                    }
                    // Append the laser to the vector
                    newEMfields->emBoundCond[iBC]->vecLaser.push_back( laser );
                }
            }
        }
        
        // -----------------
        // Clone ExternalFields properties
        // -----------------
        for( unsigned int n_extfield = 0; n_extfield < EMfields->extFields.size(); n_extfield++ ) {
            ExtField extField;
            extField.field   = EMfields->extFields[n_extfield].field;
            extField.profile = EMfields->extFields[n_extfield].profile;
            extField.index   = EMfields->extFields[n_extfield].index;
            newEMfields->extFields.push_back( extField );
        }
        
        // -----------------
        // Clone PrescribedFields properties
        // -----------------
        for( unsigned int n_extfield = 0; n_extfield < EMfields->extTimeFields.size(); n_extfield++ ) {
            ExtTimeField extField;
            extField.savedField   = EMfields->extTimeFields[n_extfield].savedField;
            extField.profile = EMfields->extTimeFields[n_extfield].profile;
            extField.index   = EMfields->extTimeFields[n_extfield].index;
            newEMfields->extTimeFields.push_back( extField );
        }
        
        // -----------------
        // Clone Antenna properties
        // -----------------
        for( unsigned int n_antenna = 0; n_antenna < EMfields->antennas.size(); n_antenna++ ) {
            Antenna antenna;
            antenna.field = NULL;
            antenna.fieldName     = EMfields->antennas[n_antenna].fieldName    ;
            antenna.space_profile = EMfields->antennas[n_antenna].space_profile;
            antenna.time_profile  = EMfields->antennas[n_antenna].time_profile ;
            antenna.index         = EMfields->antennas[n_antenna].index        ;
            newEMfields->antennas.push_back( antenna );
        }
        
        //newEMfields->finishInitialization(vecSpecies.size(), patch);
        
        return newEMfields;
    }
    
};

#endif
