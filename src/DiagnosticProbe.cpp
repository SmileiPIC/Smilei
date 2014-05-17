#include "DiagnosticProbe.h"

#include <iomanip>
#include <string>
#include <iomanip>
#include <sstream>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"

using namespace std;

DiagnosticProbe::DiagnosticProbe(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi):
smpi_(smpi), probeSize(6), fileId(0) {
    
    if (diagParams->probeStruc.size() >0) {
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId = H5Fcreate( "Probes.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);
        
        string ver(__VERSION);
        hid_t aid3  = H5Screate(H5S_SCALAR);
        hid_t atype = H5Tcopy(H5T_C_S1);
        H5Tset_size(atype, ver.size());
        H5Tset_strpad(atype,H5T_STR_NULLTERM);
        hid_t attr3 = H5Acreate(fileId, "Version", atype, aid3, H5P_DEFAULT, H5P_DEFAULT);
        
        H5Awrite(attr3, atype, ver.c_str());
        
        H5Aclose(attr3);
        H5Sclose(aid3);
        H5Tclose(atype);    

        every.resize(diagParams->probeStruc.size());
        probeParticles.resize(diagParams->probeStruc.size());
        probeId.resize(diagParams->probeStruc.size());
        
        for (unsigned int np=0; np<diagParams->probeStruc.size(); np++) {
            every[np]=diagParams->probeStruc[np].every;
            unsigned int dimProbe=diagParams->probeStruc[np].dim+2;
            unsigned int ndim=params->nDim_particle;
            
            vector<unsigned int> vecNumber=diagParams->probeStruc[np].number;
            
            unsigned int totPart=1;
            for (unsigned int iDimProb=0; iDimProb<diagParams->probeStruc[np].dim; iDimProb++) {
                totPart *= vecNumber[iDimProb];
            }
            
            probeParticles[np].initialize(totPart, ndim);
            probeId[np].resize(totPart);
            
            vector<double> partPos(ndim*totPart,0.0);
            
            for(unsigned int ipart=0; ipart!=totPart; ++ipart) {
                
                int found=smpi_->getRank();
                
                for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
                    unsigned int k=iDim+ipart*ndim;
                    
                    partPos[k]=diagParams->probeStruc[np].pos[0][iDim];
                    
                    for (unsigned int iDimProb=0; iDimProb<diagParams->probeStruc[np].dim; iDimProb++) {
                        partPos[k] += (ipart%vecNumber[iDimProb])*(diagParams->probeStruc[np].pos[iDimProb+1][iDim]-diagParams->probeStruc[np].pos[iDimProb][iDim])/(vecNumber[iDimProb]-1);
                    }
                    probeParticles[np].position(iDim,ipart) = 2*M_PI*partPos[k];
                    
                    if(smpi->getDomainLocalMin(iDim) >  probeParticles[np].position(iDim,ipart) || smpi->getDomainLocalMax(iDim) <= probeParticles[np].position(iDim,ipart)) {
                        found=-1;
                    }
                }
                probeId[np][ipart] = found;
            }
            vector<hsize_t> dims(dimProbe);
            vector<hsize_t> max_dims(dimProbe);
            vector<hsize_t> chunk_dims(dimProbe);
            
            dims[0]=0;
            max_dims[0]=H5S_UNLIMITED;
            chunk_dims[0]=1;
            
            for (unsigned int iDimProb=0; iDimProb<vecNumber.size(); iDimProb++) {
                dims[iDimProb+1]=vecNumber[iDimProb];
                max_dims[iDimProb+1]=vecNumber[iDimProb];
                chunk_dims[iDimProb+1]=1;
            }
            dims.back()=probeSize;
            max_dims.back()=probeSize;
            chunk_dims.back()=probeSize;
            
            hid_t sid = H5Screate_simple(dimProbe, &dims[0], &max_dims[0]);
            hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_layout(plist, H5D_CHUNKED);
            H5Pset_chunk(plist, dimProbe, &chunk_dims[0]);
            
            hid_t did = H5Dcreate(fileId, probeName(np).c_str(), H5T_NATIVE_FLOAT, sid, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Pclose(plist);
            H5Sclose(sid);
            
            unsigned int vecNumberProd=1;
            for (unsigned int iDimProb=0; iDimProb<vecNumber.size(); iDimProb++) {
                vecNumberProd*=vecNumber[iDimProb];
            }    
            
            vector<hsize_t> dimsPos(1+vecNumber.size());
            for (unsigned int iDimProb=0; iDimProb<vecNumber.size(); iDimProb++) {
                dimsPos[iDimProb]=vecNumber[iDimProb];
            }
            dimsPos[vecNumber.size()]=ndim;
            
            hid_t dataspace_id = H5Screate_simple(dimsPos.size(), &dimsPos[0], NULL);
            
            hid_t attribute_id = H5Acreate (did, "position", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &partPos[0]);
            H5Aclose(attribute_id);
            H5Sclose(dataspace_id);
            
            sid = H5Screate(H5S_SCALAR);	
            hid_t aid = H5Acreate(did, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(aid, H5T_NATIVE_UINT, &every[np]);
            H5Sclose(sid);
            H5Aclose(aid);
            
            sid = H5Screate(H5S_SCALAR);	
            aid = H5Acreate(did, "dimension", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(aid, H5T_NATIVE_UINT, &diagParams->probeStruc[np].dim);
            H5Sclose(sid);
            H5Aclose(aid);
            
            H5Dclose(did);
        }
        
    }
}

void DiagnosticProbe::close() {
    if (fileId>0) {
        H5Fclose(fileId);
    }
}

string DiagnosticProbe::probeName(int p) {
    ostringstream prob_name("");
    prob_name << "p" << setfill('0') << setw(4) << p;
    return prob_name.str();
}

void DiagnosticProbe::run(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp) {
    for (unsigned int np=0; np<every.size(); np++) {
        if (every[np] && timestep % every[np] == 0) {
            vector<double> data(probeSize);
            
            hid_t did = H5Dopen2(fileId, probeName(np).c_str(), H5P_DEFAULT);
            
            // All rank open all probes dataset
            hid_t sid = H5Dget_space(did);
            
            hsize_t dimProbe = H5Sget_simple_extent_dims(sid, NULL, NULL);
            
            vector<hsize_t> dims(dimProbe);
            vector<hsize_t> nulldims(dimProbe);
            for (unsigned int iDimProb=0; iDimProb<dimProbe-1; iDimProb++) {
                dims[iDimProb]=1;
                nulldims[iDimProb]=0;
            }
            dims.back()=probeSize;
            nulldims.back()=0;
            
            hid_t  partMemSpace = H5Screate_simple(dimProbe, &dims[0], NULL);
            hid_t  partMemSpaceNull = H5Screate_simple(dimProbe, &nulldims[0], NULL);
            
            // Get dataset existing dims
            vector<hsize_t> dimsO(dimProbe);
            
            H5Sget_simple_extent_dims(sid, &dimsO[0], NULL);
            
            // Increment dataset size
            dimsO[0]++;
            
            H5Dset_extent(did, &dimsO[0]);
            
            sid = H5Dget_space(did);
            
            for (unsigned int iprob=0; iprob <probeParticles[np].size(); iprob++) {
                
                vector<hsize_t>  count(dimProbe);
                if  (probeId[np][iprob]==smpi_->getRank()) {
                    for (unsigned int iDimProb=0; iDimProb< dimProbe-1; iDimProb++) {
                        count[iDimProb]=1;
                    }
                    count.back()=probeSize;
                } else {
                    for (unsigned int iDimProb=0; iDimProb< dimProbe; iDimProb++) {
                        count[iDimProb]=0;
                    }
                }
                
                vector<hsize_t> start(dimProbe);
                start[0]=dimsO[0]-1;
                for (unsigned int iDimProb=1; iDimProb<dimProbe-1; iDimProb++) {
                    start[iDimProb]=iprob % dimsO[iDimProb];
                }
                start.back()=0;
                
                H5Sselect_hyperslab(sid, H5S_SELECT_SET, &start[0], NULL, &count[0], NULL);
                
                hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
                
                if  (probeId[np][iprob]==smpi_->getRank()) {
                    (*interp)(EMfields,probeParticles[np],iprob,&Eloc_fields,&Bloc_fields);
                    
                    //! here we fill the probe data!!!
                    data[0]=Eloc_fields.x;
                    data[1]=Eloc_fields.y;
                    data[2]=Eloc_fields.z;
                    data[3]=Bloc_fields.x;
                    data[4]=Bloc_fields.y;
                    data[5]=Bloc_fields.z;
                    
                    H5Dwrite(did, H5T_NATIVE_DOUBLE, partMemSpace, sid, write_plist,&data[0]);
                } else {
                    H5Dwrite(did, H5T_NATIVE_DOUBLE, partMemSpaceNull, sid, write_plist,&data[0]);
                }
                
                H5Pclose( write_plist );
                
            }
            H5Sclose(sid);
            H5Dclose(did);
            H5Sclose(partMemSpaceNull);
            H5Sclose(partMemSpace);
        }
    }
    H5Fflush(fileId, H5F_SCOPE_GLOBAL );
}
