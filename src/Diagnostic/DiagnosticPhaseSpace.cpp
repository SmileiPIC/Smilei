#include "DiagnosticPhaseSpace.h"

#include <iomanip>
#include <string>
#include <iomanip>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"

#include "DiagnosticPhaseMomMom.h"
#include "DiagnosticPhasePosLor.h"
#include "DiagnosticPhasePosMom.h"

using namespace std;

DiagnosticPhaseSpace::~DiagnosticPhaseSpace() {
}

void DiagnosticPhaseSpace::close() {
    //! check if we're on the master (the only one that opened the file)
	if (fileId != 0) {
        H5Fclose(fileId);
	}
}



//DiagnosticPhaseSpace::DiagnosticPhaseSpace() : fileId(0) {}
DiagnosticPhaseSpace::DiagnosticPhaseSpace(Params &params, Patch* patch) : fileId(0), ndim(params.nDim_particle) {
    unsigned int numPhases=PyTools::nComponents("DiagPhase");

    //! create the particle structure
    my_part.pos.resize(ndim);
    my_part.mom.resize(3);

    for (unsigned int n_phase = 0 ; n_phase < numPhases ; n_phase++) {

	// ----- Start of extract data ------ //
	phaseStructure my_phase;
	bool ok;
        
        my_phase.every=0;
        ok=PyTools::extract("every",my_phase.every,"DiagPhase",n_phase);
        if (!ok) {
            //            if (n_probephase>0) {
            //                my_phase.every=phases.vecDiagPhase.end()->every;
            //            } else {
            my_phase.every=params.global_every;
            //            }
        }
        
        vector<string> kind;
        PyTools::extract("kind",kind,"DiagPhase",n_phase);        
        for (vector<string>::iterator it=kind.begin(); it!=kind.end();it++) {
            if (std::find(kind.begin(), it, *it) == it) {
                my_phase.kind.push_back(*it); 
            } else {
                WARNING("removed duplicate " << *it << " in \"DiagPhase\" " << n_phase);
            }
        }
        
        vector<double> time_range(2,0.);
        ok=PyTools::extract("time_range",time_range,"DiagPhase",n_phase);        
        
        if (!ok) { 
            my_phase.tmin = 0.;
            my_phase.tmax = params.sim_time;
        }
        else {
            my_phase.tmin = time_range[0];
            my_phase.tmax = time_range[1];
        }
        
        
        PyTools::extract("species",my_phase.species,"DiagPhase",n_phase);
        
        my_phase.deflate=0;
        PyTools::extract("deflate",my_phase.deflate,"DiagPhase",n_phase);
        
        if (my_phase.species.size()==0) {
            WARNING("adding all species to the \"DiagPhase\" " << n_phase);
            for (unsigned int i=0;i<params.species_param.size(); i++) {
                my_phase.species.push_back(params.species_param[i].species_type);
            }
        }
        
        PyTools::extract("pos_min",my_phase.pos_min,"DiagPhase",n_phase);
	// ----- End of extract data ------ //


	hid_t gidParent=0;
	if (patch->isMaster() ) {
	    if (n_phase==0) {
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
	    groupName << "ps" << setw(4) << setfill('0') << n_phase;
	    gidParent = H5Gcreate(fileId, groupName.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 

	    hid_t sid = H5Screate(H5S_SCALAR);	
	    hid_t aid = H5Acreate(gidParent, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
	    H5Awrite(aid, H5T_NATIVE_UINT, &my_phase.every);
	    H5Sclose(sid);
	    H5Aclose(aid);
            
	}
        
	for (unsigned int ii=0 ; ii < my_phase.kind.size(); ii++) {
	    DiagnosticPhase *diagPhase=NULL;
                        
	    // create DiagnosticPhase
	    if (params.geometry == "1d3v") {
		if (my_phase.kind[ii] == "xpx") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,0);
		} else if (my_phase.kind[ii] == "xpy") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,1);
		} else if (my_phase.kind[ii] == "xpz") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,2);
		} else if (my_phase.kind[ii] == "xlor") {
		    diagPhase =  new DiagnosticPhasePosLor(my_phase,0);
		} else if (my_phase.kind[ii] == "pxpy") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,1);
		} else if (my_phase.kind[ii] == "pxpz") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,2);
		} else if (my_phase.kind[ii] == "pypz") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,1,2);
		} else {
		    ERROR("kind " << my_phase.kind[ii] << " not implemented for geometry " << params.geometry);
		}
	    } else if (params.geometry == "2d3v") {
		if (my_phase.kind[ii] == "xpx") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,0);
		} else if (my_phase.kind[ii] == "xpy") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,1);
		} else if (my_phase.kind[ii] == "xpz") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,2);
		} else if (my_phase.kind[ii] == "ypx") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,0);
		} else if (my_phase.kind[ii] == "ypy") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,1);
		} else if (my_phase.kind[ii] == "ypz") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,2);
		} else if (my_phase.kind[ii] == "pxpy") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,1);
		} else if (my_phase.kind[ii] == "pxpz") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,2);
		} else if (my_phase.kind[ii] == "pypz") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,1,2);                    
		} else if (my_phase.kind[ii] == "xlor") {
		    diagPhase =  new DiagnosticPhasePosLor(my_phase,0);
		} else if (my_phase.kind[ii] == "ylor") {
		    diagPhase =  new DiagnosticPhasePosLor(my_phase,1);
		} else {
		    ERROR("kind " << my_phase.kind[ii] << " not implemented for geometry " << params.geometry);
		}
	    } else if (params.geometry == "3d3v") {
		if (my_phase.kind[ii] == "xpx") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,0);
		} else if (my_phase.kind[ii] == "xpy") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,1);
		} else if (my_phase.kind[ii] == "xpz") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,2);
		} else if (my_phase.kind[ii] == "ypx") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,0);
		} else if (my_phase.kind[ii] == "ypy") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,1);
		} else if (my_phase.kind[ii] == "ypz") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,2);
		} else if (my_phase.kind[ii] == "zpx") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,2,0);
		} else if (my_phase.kind[ii] == "zpy") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,2,1);
		} else if (my_phase.kind[ii] == "zpz") {
		    diagPhase =  new DiagnosticPhasePosMom(my_phase,2,2);
		} else if (my_phase.kind[ii] == "pxpy") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,1);
		} else if (my_phase.kind[ii] == "pxpz") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,2);
		} else if (my_phase.kind[ii] == "pypz") {
		    diagPhase =  new DiagnosticPhaseMomMom(my_phase,1,2);                    
		} else if (my_phase.kind[ii] == "xlor") {
		    diagPhase =  new DiagnosticPhasePosLor(my_phase,0);
		} else if (my_phase.kind[ii] == "ylor") {
		    diagPhase =  new DiagnosticPhasePosLor(my_phase,1);
		} else if (my_phase.kind[ii] == "zlor") {
		    diagPhase =  new DiagnosticPhasePosLor(my_phase,2);
		} else {
		    ERROR("kind " << my_phase.kind[ii] << " not implemented for geometry " << params.geometry);
		}                
	    } else {
		ERROR("DiagnosticPhase not implemented for geometry " << params.geometry);
	    }
            
	    if (diagPhase) {
		if (patch->isMaster()) {
		    //! create a group for each species of this diag and keep track of its ID.
                    
		    hsize_t dims[3] = {0,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
		    hsize_t max_dims[3] = {H5S_UNLIMITED,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
		    hsize_t chunk_dims[3] = {1,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                    
		    hid_t sid = H5Screate_simple (3, dims, max_dims);	
		    hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
		    H5Pset_layout(pid, H5D_CHUNKED);
		    H5Pset_chunk(pid, 3, chunk_dims);
                    
		    H5Pset_deflate (pid, std::min((unsigned int)9,my_phase.deflate));
                    
		    hid_t did = H5Dcreate (gidParent, my_phase.kind[ii].c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid,H5P_DEFAULT);
		    H5Pclose (pid);	
		    H5Sclose (sid);
                    
		    // write attribute of species present in the phaseSpace
		    string namediag;
		    for (unsigned int k=0; k<my_phase.species.size(); k++) {
			namediag+=my_phase.species[k]+" ";
		    }
		    namediag=namediag.substr(0, namediag.size()-1);
		    sid = H5Screate(H5S_SCALAR);
		    hid_t tid = H5Tcopy(H5T_C_S1);
		    H5Tset_size(tid, namediag.size());
		    H5Tset_strpad(tid,H5T_STR_NULLTERM);
		    hid_t aid = H5Acreate(gidParent, "species", tid, sid, H5P_DEFAULT, H5P_DEFAULT);
		    H5Awrite(aid, tid, namediag.c_str());
                    
                    
		    // write attribute extent of the phaseSpace
		    hsize_t dimsPos[2] = {2,2};
		    sid = H5Screate_simple(2, dimsPos, NULL);
		    aid = H5Acreate (gidParent, "extents", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
		    double tmp[4] = {diagPhase->firstmin, diagPhase->firstmax, diagPhase->secondmin, diagPhase->secondmax};
		    H5Awrite(aid, H5T_NATIVE_DOUBLE, tmp);
		    H5Aclose(aid);
		    H5Sclose(sid);
                    
		    diagPhase->dataId=did;
                    
		}
		vecDiagPhase.push_back(diagPhase);	
	    }
	}
	if (patch->isMaster() ) {
	    H5Gclose(gidParent);
	}
    }
}

void DiagnosticPhaseSpace::run(int timestep, std::vector<Species*>& vecSpecies) {
    //! check which diagnosticPhase to run at this timestep
    vector<DiagnosticPhase*> vecDiagPhaseActiveTimestep;	
    for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhase.begin() ; diag != vecDiagPhase.end(); diag++) {
	if (timestep % (*diag)->every==0) vecDiagPhaseActiveTimestep.push_back(*diag);
    }
	
    if (vecDiagPhaseActiveTimestep.size()>0) {
	vecDiagPhaseToRun.clear();

	for (vector<Species*>::const_iterator mySpec=vecSpecies.begin(); mySpec!= vecSpecies.end(); mySpec++) {
	    if (!(*mySpec)->particles->isTestParticles) {
			
		//! check which diagnosticPhase to run for the species 
		for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseActiveTimestep.begin() ; diag != vecDiagPhaseActiveTimestep.end(); diag++) {
		    if(find((*diag)->my_species.begin(), (*diag)->my_species.end(), (*mySpec)->species_param.species_type) != (*diag)->my_species.end()) { 
			vecDiagPhaseToRun.push_back(*diag);
		    }
		}

		if (vecDiagPhaseToRun.size()>0) {
                
		    //! cycle over all the particles
		    for (unsigned int ibin = 0 ; ibin < (*mySpec)->bmin.size() ; ibin++) {
			for (int iPart=(*mySpec)->bmin[ibin] ; iPart<(*mySpec)->bmax[ibin]; iPart++ ) {
			    //! fill the my_part structure
			    for(unsigned int k=0;k<ndim;k++) {
				my_part.pos[k]=(*mySpec)->particles->position(k,iPart);
			    }
			    for(unsigned int k=0;k<3;k++) {
				my_part.mom[k]=(*mySpec)->particles->momentum(k,iPart);
			    }
			    my_part.weight=(*mySpec)->particles->weight(iPart);
			    my_part.charge=(*mySpec)->particles->charge(iPart);
			    for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseToRun.begin() ; diag != vecDiagPhaseToRun.end(); diag++) {
				//! do something with each particle
				(*diag)->run(my_part);
			    }						
			}
		    }
		}
	    }
	} // END for (vector<Species*>::const_iterator mySpec
	/*
        //! and finally write the data (reduce data on 1 proc, write it and clear memory for future usage)
        for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseActiveTimestep.begin() ; diag != vecDiagPhaseActiveTimestep.end(); diag++) {
	(*diag)->writeData();
	}*/
    }
}
