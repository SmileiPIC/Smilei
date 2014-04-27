#include "DiagnosticPhaseSpace.h"

#include <iomanip>
#include <string>
#include <iomanip>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"
#include "DiagnosticPhase2DxP.h"

using namespace std;

DiagnosticPhaseSpace::~DiagnosticPhaseSpace() {
}

void DiagnosticPhaseSpace::close() {
    //! check if we're on the master (the only one that opened the file)
	if (fileId != 0) {
        //!close all hdf5  groups (in principle this is done also with H5Fclose below...)
        for (std::map<DiagnosticPhase*, std::map<std::string,hid_t> >::iterator iterMap=mapGroupId.begin(); iterMap!= mapGroupId.end(); iterMap++) {
            for (map<string, hid_t>::iterator iter = iterMap->second.begin(); iter != iterMap->second.end(); iter++) {
                H5Gclose(iter->second);
            }
        }
        //! close the hdf5 file
        H5Fclose(fileId);
	}
}



DiagnosticPhaseSpace::DiagnosticPhaseSpace(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi) : fileId(0), ndim(params->nDim_particle) {
	for (unsigned int i =0 ; i < diagParams->vecPhase.size(); i++) {
		DiagnosticPhase *diagPhase=NULL;
		
        hid_t gidParent=0;
        
		if (smpi->isMaster()) {
            //!open hdf5 file
			if (i==0) {
				ostringstream file_name("");
				file_name<<"PhaseSpace.h5";
				fileId = H5Fcreate( file_name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                string ver(__VERSION);
                
                // write version
                hid_t aid3  = H5Screate(H5S_SCALAR);
                hid_t atype = H5Tcopy(H5T_C_S1);
                H5Tset_size(atype, ver.size());
                H5Tset_strpad(atype,H5T_STR_NULLTERM);
                hid_t attr3 = H5Acreate2(fileId, "Version", atype, aid3, H5P_DEFAULT, H5P_DEFAULT);
                
                H5Awrite(attr3, atype, ver.c_str());
                
                H5Aclose(attr3);
                H5Sclose(aid3);
                H5Tclose(atype);
			}
			ostringstream groupName("");
			groupName << "ps_" << i << "_" << diagParams->vecPhase[i].kind;
            gidParent = H5Gcreate(fileId, groupName.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);			
		}
        
        // create DiagnosticPhase
		if (params->geometry == "1d3v") {
			if (diagParams->vecPhase[i].kind == "xpx") {
				diagPhase =  new DiagnosticPhase2DxP(diagParams->vecPhase[i],0);
            } else if (diagParams->vecPhase[i].kind == "xpy") {
				diagPhase =  new DiagnosticPhase2DxP(diagParams->vecPhase[i],1);
            } else if (diagParams->vecPhase[i].kind == "xpz") {
				diagPhase =  new DiagnosticPhase2DxP(diagParams->vecPhase[i],2);
            } else {
                ERROR("kind " << diagParams->vecPhase[i].kind << " not implemented for geometry " << params->geometry);
            }
		} else {
			ERROR("DiagnosticPhase not implemented for geometry " << params->geometry);
		}
		
		if (diagPhase) {
            if (smpi->isMaster()) {
                //! create a group for each species of this diag and keep track of its ID.
                map<string,hid_t> localmap;
                for (unsigned int k=0; k<diagParams->vecPhase[i].species.size(); k++) {
                    hid_t gid = H5Gcreate(gidParent, diagParams->vecPhase[i].species[k].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    localmap[diagParams->vecPhase[i].species[k]]=gid;
                }
                mapGroupId[diagPhase]=localmap;
            }
			vecDiagPhase.push_back(diagPhase);	
		}
	}
}

void DiagnosticPhaseSpace::run(int timestep, std::vector<Species*>& vecSpecies) {
	//! check which diagnosticPhase to run at this timestep
	vector<DiagnosticPhase*> vecDiagPhaseActiveTimestep;	
	for (unsigned int i =0 ; i < vecDiagPhase.size(); i++) {
		if (timestep%vecDiagPhase[i]->every==0) vecDiagPhaseActiveTimestep.push_back(vecDiagPhase[i]);
	}
	
	if (vecDiagPhaseActiveTimestep.size()>0) {
		for (unsigned int j=0; j < vecSpecies.size(); j++) {
			
			//! check which diagnosticPhase to run for the species 
			vector<DiagnosticPhase*> vecDiagPhaseToRun;
			for (unsigned int i =0 ; i < vecDiagPhaseActiveTimestep.size(); i++) {
				if(find(vecDiagPhaseActiveTimestep[i]->my_species.begin(), vecDiagPhaseActiveTimestep[i]->my_species.end(), vecSpecies[j]->name_str) != vecDiagPhaseActiveTimestep[i]->my_species.end()) { 
					vecDiagPhaseToRun.push_back(vecDiagPhaseActiveTimestep[i]);
				}
			}
			
			if (vecDiagPhaseToRun.size()>0) {
                partStruct my_part;
                my_part.pos.resize(ndim);
                my_part.mom.resize(3);
                
				//! cycle over all the particles
				for (unsigned int ibin = 0 ; ibin < vecSpecies[j]->bmin.size() ; ibin++) {
					for (int iPart=vecSpecies[j]->bmin[ibin] ; iPart<vecSpecies[j]->bmax[ibin]; iPart++ ) {
						for (unsigned int i =0 ; i < vecDiagPhaseToRun.size(); i++) {
                            //! fill the my_part structure
							for(unsigned int k=0;k<ndim;k++) {
								my_part.pos[k]=vecSpecies[j]->particles.position(k,iPart);
							}
							for(unsigned int k=0;k<3;k++) {
								my_part.mom[k]=vecSpecies[j]->particles.momentum(k,iPart);
							}
							my_part.weight=vecSpecies[j]->particles.weight(iPart);
							my_part.charge=vecSpecies[j]->particles.charge(iPart);
                            //! do something with each partcle
							vecDiagPhaseToRun[i]->doSomething(my_part);
						}						
					}
				}
                //! and finally write the data (reduce data on 1 proc, write it and clear memory for future usage)
				for (unsigned int i =0 ; i < vecDiagPhaseToRun.size(); i++) {
					vecDiagPhaseToRun[i]->writeData(timestep, mapGroupId[vecDiagPhaseToRun[i]][vecSpecies[j]->name_str]);
				}
                if(fileId>0) H5Fflush(fileId, H5F_SCOPE_GLOBAL);
			}
			
		}
	}
}
