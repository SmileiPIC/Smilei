#include "DiagnosticPhaseSpace.h"

#include <iomanip>
#include <string>
#include <iomanip>

#include "Params.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"
#include "H5.h"

#include "DiagnosticPhaseMomMom.h"
#include "DiagnosticPhasePosLor.h"
#include "DiagnosticPhasePosMom.h"


using namespace std;

DiagnosticPhaseSpace::DiagnosticPhaseSpace(Params& params, SmileiMPI *smpi) :
fileId(0)
{
    HEREIAM("");
    //! create the particle structure
    my_part.pos.resize(params.nDim_particle);
    my_part.mom.resize(3);
    
    unsigned int numPhases=PyTools::nComponents("DiagPhase");
    for (unsigned int n_phase = 0; n_phase < numPhases; n_phase++) {
        HEREIAM("");
        MESSAGE(1,"Activating DiagPhase " << n_phase);
        
        vector<string> kind;
        if (!PyTools::extract("kind",kind,"DiagPhase",n_phase)) {
            ERROR("For DiagPhase " <<n_phase << " missing kind");
        }
        sort( kind.begin(), kind.end() );
        kind.erase( unique( kind.begin(), kind.end() ), kind.end() );
        hid_t gidParent=0;
        if (smpi->isMaster()) {
            if (n_phase == 0) {
                ostringstream file_name("");
                file_name<<"PhaseSpace.h5";
                fileId = H5Fcreate( file_name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                // write version
                H5::attr(fileId, "Version", string(__VERSION));
            }
            
            ostringstream groupName("");
            groupName << "ps" << setw(4) << setfill('0') << n_phase;
            gidParent = H5::group(fileId, groupName.str());
        }
        
        for (unsigned int ii=0 ; ii < kind.size(); ii++) {
            DiagnosticPhase *diagPhase=NULL;
            
            // create DiagnosticPhase
            if (kind[ii] == "xpx") {
                diagPhase =  new DiagnosticPhasePosMom(params,n_phase,0,0);
            } else if (kind[ii] == "xpy") {
                diagPhase =  new DiagnosticPhasePosMom(params,n_phase,0,1);
            } else if (kind[ii] == "xpz") {
                diagPhase =  new DiagnosticPhasePosMom(params,n_phase,0,2);
            } else if (kind[ii] == "xlor") {
                diagPhase =  new DiagnosticPhasePosLor(params,n_phase,0);
            } else if (kind[ii] == "pxpy") {
                diagPhase =  new DiagnosticPhaseMomMom(params,n_phase,0,1);
            } else if (kind[ii] == "pxpz") {
                diagPhase =  new DiagnosticPhaseMomMom(params,n_phase,0,2);
            } else if (kind[ii] == "pypz") {
                diagPhase =  new DiagnosticPhaseMomMom(params,n_phase,1,2);
            }
            
            if (params.geometry == "2d3v" || params.geometry == "3d3v") {
                if (kind[ii] == "ypx") {
                    diagPhase =  new DiagnosticPhasePosMom(params,n_phase,1,0);
                } else if (kind[ii] == "ypy") {
                    diagPhase =  new DiagnosticPhasePosMom(params,n_phase,1,1);
                } else if (kind[ii] == "ypz") {
                    diagPhase =  new DiagnosticPhasePosMom(params,n_phase,1,2);
                } else if (kind[ii] == "ylor") {
                    diagPhase =  new DiagnosticPhasePosLor(params,n_phase,1);
                }
            }
            
            if (params.geometry == "3d3v") {
                if (kind[ii] == "zpx") {
                    diagPhase =  new DiagnosticPhasePosMom(params,n_phase,2,0);
                } else if (kind[ii] == "zpy") {
                    diagPhase =  new DiagnosticPhasePosMom(params,n_phase,2,1);
                } else if (kind[ii] == "zpz") {
                    diagPhase =  new DiagnosticPhasePosMom(params,n_phase,2,2);
                } else if (kind[ii] == "zlor") {
                    diagPhase =  new DiagnosticPhasePosLor(params,n_phase,2);
                }
            }
            
            if (diagPhase) {
                vecDiagPhase.push_back(diagPhase);
                
                if (smpi->isMaster()) {
                    //! create data for each species of this diag and keep track of its ID.
                    
                    hsize_t dims[3] = {0,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                    hsize_t max_dims[3] = {H5S_UNLIMITED,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                    hsize_t chunk_dims[3] = {1,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                    
                    hid_t sid = H5Screate_simple (3, dims, max_dims);
                    hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
                    H5Pset_layout(pid, H5D_CHUNKED);
                    H5Pset_chunk(pid, 3, chunk_dims);
                    
                    unsigned int deflate;
                    PyTools::extract("deflate",deflate,"DiagPhase",n_phase);
                    
                    H5Pset_deflate(pid, std::min((unsigned int)9,deflate));
                    diagPhase->dataId = H5Dcreate (gidParent, kind[ii].c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid,H5P_DEFAULT);
                    H5Pclose (pid);
                    H5Sclose (sid);
                    
                    // write attribute of species present in the phaseSpace
                    string namediag;
                    for (unsigned int k=0; k<diagPhase->my_species.size(); k++) {
                        namediag+=diagPhase->my_species[k]+" ";
                    }
                    namediag=namediag.substr(0, namediag.size()-1);
                    H5::attr(diagPhase->dataId,"species",namediag);
                    
                    // write attribute extent of the phaseSpace
                    hsize_t dimsPos[2] = {2,2};
                    sid = H5Screate_simple(2, dimsPos, NULL);
                    hid_t aid = H5Acreate (diagPhase->dataId, "extents", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
                    double tmp[4] = {diagPhase->firstmin, diagPhase->firstmax, diagPhase->secondmin, diagPhase->secondmax};
                    H5Awrite(aid, H5T_NATIVE_DOUBLE, tmp);
                    H5Aclose(aid);
                    H5Sclose(sid);
                    
                }
            } else {
                ERROR("DiagnosticPhase not implemented for geometry " << params.geometry);
            }
            
        }
        
        if (smpi->isMaster() ) {
            H5Gclose(gidParent);
        }
    }

}

void DiagnosticPhaseSpace::close() {
    //! check if we're on the master (the only one that opened the file)
    if (fileId != 0) {
        H5Fclose(fileId);
    }
}

void DiagnosticPhaseSpace::run(int timestep, std::vector<Species*>& vecSpecies) {
    //! check which diagnosticPhase to run at this timestep
    vector<DiagnosticPhase*> vecDiagPhaseActiveTimestep;
    for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhase.begin() ; diag != vecDiagPhase.end(); diag++) {
        if (timestep % (*diag)->every==0) vecDiagPhaseActiveTimestep.push_back(*diag);
    }
    
    if (vecDiagPhaseActiveTimestep.size()>0) {
        for (vector<Species*>::const_iterator mySpec=vecSpecies.begin(); mySpec!= vecSpecies.end(); mySpec++) {
            if (!(*mySpec)->particles.isTestParticles) {
                
                //! check which diagnosticPhase to run for the species
                vector<DiagnosticPhase*> vecDiagPhaseToRun;
                for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseActiveTimestep.begin() ; diag != vecDiagPhaseActiveTimestep.end(); diag++) {
                    if(find((*diag)->my_species.begin(), (*diag)->my_species.end(), (*mySpec)->species_type) != (*diag)->my_species.end()) {
                        vecDiagPhaseToRun.push_back(*diag);
                    }
                }
                if (vecDiagPhaseToRun.size()>0) {
                    
                    //! cycle over all the particles
                    for (unsigned int ibin = 0 ; ibin < (*mySpec)->bmin.size() ; ibin++) {
                        for (int iPart=(*mySpec)->bmin[ibin] ; iPart<(*mySpec)->bmax[ibin]; iPart++ ) {
                            //! fill the my_part structure
                            for(unsigned int k=0;k<my_part.pos.size();k++) {
                                my_part.pos[k]=(*mySpec)->particles.position(k,iPart);
                            }
                            for(unsigned int k=0;k<3;k++) {
                                my_part.mom[k]=(*mySpec)->particles.momentum(k,iPart);
                            }
                            my_part.weight=(*mySpec)->particles.weight(iPart);
                            my_part.charge=(*mySpec)->particles.charge(iPart);
                            for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseToRun.begin() ; diag != vecDiagPhaseToRun.end(); diag++) {
                                //! do something with each particle
                                (*diag)->run(my_part);
                            }						
                        }
                    }
                }
            }
        } // END for (vector<Species*>::const_iterator mySpec
        //! and finally write the data (reduce data on 1 proc, write it and clear memory for future usage)
        for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseActiveTimestep.begin() ; diag != vecDiagPhaseActiveTimestep.end(); diag++) {
            (*diag)->writeData();
        }
    }
}
