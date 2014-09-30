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
#include "DiagnosticPhasePosMom.h"
#include "DiagnosticPhasePosLor.h"
#include "DiagnosticPhaseMomMom.h"

using namespace std;

DiagnosticPhaseSpace::~DiagnosticPhaseSpace() {
}

void DiagnosticPhaseSpace::close() {
    //! check if we're on the master (the only one that opened the file)
	if (fileId != 0) {
        //!close all hdf5  groups (in principle this is done also with H5Fclose below...)
        for (std::map<DiagnosticPhase*, std::map<std::string,hid_t> >::iterator iterMap=mapDataId.begin(); iterMap!= mapDataId.end(); iterMap++) {
            for (map<string, hid_t>::iterator iter = iterMap->second.begin(); iter != iterMap->second.end(); iter++) {
                H5Dclose(iter->second);
            }
        }
        //! close the hdf5 file
        H5Fclose(fileId);
	}
}



DiagnosticPhaseSpace::DiagnosticPhaseSpace(PicParams &params, DiagParams &diagParams, SmileiMPI* smpi) : fileId(0), ndim(params.nDim_particle) {
	for (unsigned int i = 0 ; i < diagParams.vecPhase.size(); i++) {
        hid_t gidParent=0;
        if (smpi->isMaster() ) {
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
                hid_t attr3 = H5Acreate(fileId, "Version", atype, aid3, H5P_DEFAULT, H5P_DEFAULT);
                
                H5Awrite(attr3, atype, ver.c_str());
                
                H5Aclose(attr3);
                H5Sclose(aid3);
                H5Tclose(atype);
            }
            ostringstream groupName("");
            groupName << "ps" << setw(4) << setfill('0') << i;
            gidParent = H5Gcreate(fileId, groupName.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 

            hid_t sid = H5Screate(H5S_SCALAR);	
            hid_t aid = H5Acreate(gidParent, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(aid, H5T_NATIVE_UINT, &diagParams.vecPhase[i].every);
            H5Sclose(sid);
            H5Aclose(aid);
            
        }
        
        for (unsigned int ii=0 ; ii < diagParams.vecPhase[i].kind.size(); ii++) {
            DiagnosticPhase *diagPhase=NULL;
                        
            // create DiagnosticPhase
            if (params.geometry == "1d3v") {
                if (diagParams.vecPhase[i].kind[ii] == "xpx") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,0);
                } else if (diagParams.vecPhase[i].kind[ii] == "xpy") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "xpz") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "xlor") {
                    diagPhase =  new DiagnosticPhasePosLor(diagParams.vecPhase[i],0);
                } else if (diagParams.vecPhase[i].kind[ii] == "pxpy") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],0,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "pxpz") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],0,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "pypz") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],1,2);
                } else {
                    ERROR("kind " << diagParams.vecPhase[i].kind[ii] << " not implemented for geometry " << params.geometry);
                }
            } else if (params.geometry == "2d3v") {
                if (diagParams.vecPhase[i].kind[ii] == "xpx") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,0);
                } else if (diagParams.vecPhase[i].kind[ii] == "xpy") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "xpz") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "ypx") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],1,0);
                } else if (diagParams.vecPhase[i].kind[ii] == "ypy") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],1,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "ypz") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],1,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "pxpy") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],0,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "pxpz") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],0,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "pypz") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],1,2);                    
                } else if (diagParams.vecPhase[i].kind[ii] == "xlor") {
                    diagPhase =  new DiagnosticPhasePosLor(diagParams.vecPhase[i],0);
                } else if (diagParams.vecPhase[i].kind[ii] == "ylor") {
                    diagPhase =  new DiagnosticPhasePosLor(diagParams.vecPhase[i],1);
                } else {
                    ERROR("kind " << diagParams.vecPhase[i].kind[ii] << " not implemented for geometry " << params.geometry);
                }
            } else if (params.geometry == "3d3v") {
                if (diagParams.vecPhase[i].kind[ii] == "xpx") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,0);
                } else if (diagParams.vecPhase[i].kind[ii] == "xpy") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "xpz") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],0,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "ypx") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],1,0);
                } else if (diagParams.vecPhase[i].kind[ii] == "ypy") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],1,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "ypz") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],1,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "zpx") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],2,0);
                } else if (diagParams.vecPhase[i].kind[ii] == "zpy") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],2,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "zpz") {
                    diagPhase =  new DiagnosticPhasePosMom(diagParams.vecPhase[i],2,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "pxpy") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],0,1);
                } else if (diagParams.vecPhase[i].kind[ii] == "pxpz") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],0,2);
                } else if (diagParams.vecPhase[i].kind[ii] == "pypz") {
                    diagPhase =  new DiagnosticPhaseMomMom(diagParams.vecPhase[i],1,2);                    
                } else if (diagParams.vecPhase[i].kind[ii] == "xlor") {
                    diagPhase =  new DiagnosticPhasePosLor(diagParams.vecPhase[i],0);
                } else if (diagParams.vecPhase[i].kind[ii] == "ylor") {
                    diagPhase =  new DiagnosticPhasePosLor(diagParams.vecPhase[i],1);
                } else if (diagParams.vecPhase[i].kind[ii] == "zlor") {
                    diagPhase =  new DiagnosticPhasePosLor(diagParams.vecPhase[i],2);
                } else {
                    ERROR("kind " << diagParams.vecPhase[i].kind[ii] << " not implemented for geometry " << params.geometry);
                }                
            } else {
                ERROR("DiagnosticPhase not implemented for geometry " << params.geometry);
            }
            
            if (diagPhase) {
                if (smpi->isMaster()) {
                    //! create a group for each species of this diag and keep track of its ID.
                    map<string,hid_t> localmap;
                    hid_t gid = H5Gcreate(gidParent, diagParams.vecPhase[i].kind[ii].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);			
                    for (unsigned int k=0; k<diagParams.vecPhase[i].species.size(); k++) {
                        
                        hsize_t dims[3] = {0,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                        hsize_t max_dims[3] = {H5S_UNLIMITED,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                        hsize_t chunk_dims[3] = {1,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                        
                        hid_t sid = H5Screate_simple (3, dims, max_dims);	
                        hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
                        H5Pset_layout(pid, H5D_CHUNKED);
                        H5Pset_chunk(pid, 3, chunk_dims);

                        H5Pset_deflate (pid, std::min((unsigned int)9,diagParams.vecPhase[i].deflate));
                        
                        hid_t did = H5Dcreate (gid, diagParams.vecPhase[i].species[k].c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid,H5P_DEFAULT);
                        H5Pclose (pid);	
                        H5Sclose (sid);
                        
                        localmap[diagParams.vecPhase[i].species[k]]=did;
                    }
                    
                    hsize_t dimsPos[2] = {2,2};
                    hid_t sid = H5Screate_simple(2, dimsPos, NULL);
                    hid_t aid = H5Acreate (gid, "extents", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
                    double tmp[4] = {diagPhase->firstmin, diagPhase->firstmax, diagPhase->secondmin, diagPhase->secondmax};
                    H5Awrite(aid, H5T_NATIVE_DOUBLE, tmp);
                    H5Aclose(aid);
                    H5Sclose(sid);
                    
                    H5Gclose(gid);
                    mapDataId[diagPhase]=localmap;
                }
                vecDiagPhase.push_back(diagPhase);	
            }
        }
        if (smpi->isMaster() ) {
            H5Gclose(gidParent);
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
        //! create the particle structure
        partStruct my_part;
        my_part.pos.resize(ndim);
        my_part.mom.resize(3);

		for (unsigned int j=0; j < vecSpecies.size(); j++) {
			
			//! check which diagnosticPhase to run for the species 
			vector<DiagnosticPhase*> vecDiagPhaseToRun;
			for (unsigned int i =0 ; i < vecDiagPhaseActiveTimestep.size(); i++) {
				if(find(vecDiagPhaseActiveTimestep[i]->my_species.begin(), vecDiagPhaseActiveTimestep[i]->my_species.end(), vecSpecies[j]->name_str) != vecDiagPhaseActiveTimestep[i]->my_species.end()) { 
					vecDiagPhaseToRun.push_back(vecDiagPhaseActiveTimestep[i]);
				}
			}
			
			if (vecDiagPhaseToRun.size()>0) {
                
				//! cycle over all the particles
				for (unsigned int ibin = 0 ; ibin < vecSpecies[j]->bmin.size() ; ibin++) {
					for (int iPart=vecSpecies[j]->bmin[ibin] ; iPart<vecSpecies[j]->bmax[ibin]; iPart++ ) {
                        //! fill the my_part structure
                        for(unsigned int k=0;k<ndim;k++) {
                            my_part.pos[k]=vecSpecies[j]->particles.position(k,iPart);
                        }
                        for(unsigned int k=0;k<3;k++) {
                            my_part.mom[k]=vecSpecies[j]->particles.momentum(k,iPart);
                        }
						for (unsigned int i =0 ; i < vecDiagPhaseToRun.size(); i++) {
							my_part.weight=vecSpecies[j]->particles.weight(iPart);
							my_part.charge=vecSpecies[j]->particles.charge(iPart);
                            //! do something with each partcle
							vecDiagPhaseToRun[i]->run(my_part);
						}						
					}
				}
                //! and finally write the data (reduce data on 1 proc, write it and clear memory for future usage)
				for (unsigned int i =0 ; i < vecDiagPhaseToRun.size(); i++) {
					vecDiagPhaseToRun[i]->writeData(mapDataId[vecDiagPhaseToRun[i]][vecSpecies[j]->name_str]);
				}
			}
		}
        if(fileId>0) H5Fflush(fileId, H5F_SCOPE_GLOBAL);
	}
}
