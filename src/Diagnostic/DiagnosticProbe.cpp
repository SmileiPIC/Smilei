#include "DiagnosticProbe.h"

#include <iomanip>
#include <string>
#include <iomanip>
#include <sstream>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"
#include "SmileiMPI_Cart2D.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"
#include "DiagParams.h"

using namespace std;

DiagnosticProbe::DiagnosticProbe(PicParams &params, DiagParams &diagParams, SmileiMPI* smpi):
cpuRank((int)smpi->getRank()),
probeSize(10), 
fileId(0) {
    
    if (diagParams.probeStruc.size() >0) {
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId = H5Fcreate( "Probes.h5", H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        string ver(__VERSION);
        hid_t sid = H5Screate(H5S_SCALAR);
        hid_t tid = H5Tcopy(H5T_C_S1);
        H5Tset_size(tid, ver.size());
        H5Tset_strpad(tid,H5T_STR_NULLTERM);
        hid_t aid = H5Acreate(fileId, "Version", tid, sid, H5P_DEFAULT, H5P_DEFAULT);
        
        H5Awrite(aid, tid, ver.c_str());
        
        H5Aclose(aid);
        H5Sclose(sid);
        H5Tclose(tid);
        
        every.resize(diagParams.probeStruc.size());
        probeParticles.resize(diagParams.probeStruc.size());
        probeId.resize(diagParams.probeStruc.size());
        
        for (unsigned int np=0; np<diagParams.probeStruc.size(); np++) {
            every[np]=diagParams.probeStruc[np].every;
            unsigned int dimProbe=diagParams.probeStruc[np].dim+2;
            unsigned int ndim=params.nDim_particle;
            
            vector<unsigned int> vecNumber=diagParams.probeStruc[np].number;
            
            unsigned int totPart=1;
            for (unsigned int iDimProbe=0; iDimProbe<diagParams.probeStruc[np].dim; iDimProbe++) {
                totPart *= vecNumber[iDimProbe];
            }
            
            probeParticles[np].initialize(totPart, ndim);
            probeId[np].resize(totPart);
            
            vector<double> partPos(ndim*totPart,0.0);
            
            for(unsigned int ipart=0; ipart!=totPart; ++ipart) {
                int found=cpuRank;
                for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
                    partPos[iDim+ipart*ndim]=diagParams.probeStruc[np].pos[0][iDim];
                    // the particle position is a linear combiantion of the point pos with posFirst or posSecond or posThird
                    for (unsigned int iDimProbe=0; iDimProbe<diagParams.probeStruc[np].dim; iDimProbe++) {
                        partPos[iDim+ipart*ndim] += (ipart%vecNumber[iDimProbe])*(diagParams.probeStruc[np].pos[iDimProbe+1][iDim]-diagParams.probeStruc[np].pos[0][iDim])/(vecNumber[iDimProbe]-1);
                    }
                    probeParticles[np].position(iDim,ipart) = 2*M_PI*partPos[iDim+ipart*ndim];
                    
                    //!fixme this is awful: we add one cell if we're on the upper border
                    double maxToCheck=smpi->getDomainLocalMax(iDim);                    
                    if (ndim==1) {
                        if ((static_cast<SmileiMPI_Cart1D*>(smpi))->isEastern()) {
                            maxToCheck+=params.cell_length[iDim];
                        }
                    } else if (ndim==2) {
                        if ((iDim == 0 && (static_cast<SmileiMPI_Cart2D*>(smpi))->isEastern()) ||
                            (iDim == 1 && (static_cast<SmileiMPI_Cart2D*>(smpi))->isNorthern())) {
                            maxToCheck+=params.cell_length[iDim];
                        }                        
                    } else {
                        ERROR("implement here");
                    }

                    if (probeParticles[np].position(iDim,ipart) < smpi->getDomainLocalMin(iDim) ||
                        probeParticles[np].position(iDim,ipart) >= maxToCheck) {
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
            
            for (unsigned int iDimProbe=1; iDimProbe<dimProbe-1; iDimProbe++) {
                dims[iDimProbe]=vecNumber[iDimProbe-1];
                max_dims[iDimProbe]=vecNumber[iDimProbe-1];
                chunk_dims[iDimProbe]=1;
            }
            dims.back()=probeSize;
            max_dims.back()=probeSize;
            chunk_dims.back()=probeSize;
            
            sid = H5Screate_simple(dimProbe, &dims[0], &max_dims[0]);
            hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_layout(pid, H5D_CHUNKED);
            H5Pset_chunk(pid, dimProbe, &chunk_dims[0]);
            
            hid_t did = H5Dcreate(fileId, probeName(np).c_str(), H5T_NATIVE_FLOAT, sid, H5P_DEFAULT, pid, H5P_DEFAULT);
            H5Pclose(pid);
            H5Sclose(sid);
            
            unsigned int vecNumberProd=1;
            for (unsigned int iDimProbe=0; iDimProbe<vecNumber.size(); iDimProbe++) {
                vecNumberProd*=vecNumber[iDimProbe];
            }  
            
            vector<hsize_t> dimsPos(1+vecNumber.size());
            for (unsigned int iDimProbe=0; iDimProbe<vecNumber.size(); iDimProbe++) {
                dimsPos[iDimProbe]=vecNumber[iDimProbe];
            }
            dimsPos[vecNumber.size()]=ndim;
            
            hid_t sid = H5Screate_simple(dimsPos.size(), &dimsPos[0], NULL);
            
            hid_t aid = H5Acreate (did, "position", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(aid, H5T_NATIVE_DOUBLE, &partPos[0]);
            H5Aclose(aid);
            H5Sclose(sid);
            
            sid = H5Screate(H5S_SCALAR);	
            aid = H5Acreate(did, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(aid, H5T_NATIVE_UINT, &every[np]);
            H5Sclose(sid);
            H5Aclose(aid);
            
            sid = H5Screate(H5S_SCALAR);	
            aid = H5Acreate(did, "dimension", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(aid, H5T_NATIVE_UINT, &diagParams.probeStruc[np].dim);
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
            
            hid_t sidPart = H5Screate_simple(dimProbe, &dims[0], NULL);
            hid_t sidNull = H5Screate_simple(dimProbe, &nulldims[0], NULL);
            
            // Get dataset existing dims
            vector<hsize_t> dimsO(dimProbe);
            
            H5Sget_simple_extent_dims(sid, &dimsO[0], NULL);
            
            // Increment dataset size
            dimsO[0]++;
            
            H5Dset_extent(did, &dimsO[0]);
            
            H5Sclose(sid);
            
            sid = H5Dget_space(did);
            
            for (int iprob=0; iprob <probeParticles[np].size(); iprob++) {
                
                vector<hsize_t> count(dimProbe);
                if (probeId[np][iprob]==cpuRank) {
                    fill(count.begin(),count.end()-1,1);
                    count.back()=probeSize;
                } else {
                    fill(count.begin(),count.end(),0);
                }
                
                vector<hsize_t> start(dimProbe);
                start.front()=dimsO[0]-1;                
                for (unsigned int iDimProb=1; iDimProb<dimProbe-1; iDimProb++) {
                    start[iDimProb]=iprob % dimsO[iDimProb];
                }
                start.back()=0;
                
                H5Sselect_hyperslab(sid, H5S_SELECT_SET, &start[0], NULL, &count[0], NULL);
                
                hid_t pid = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(pid, H5FD_MPIO_INDEPENDENT);
                
                if (probeId[np][iprob]==cpuRank) {
                    (*interp)(EMfields,probeParticles[np],iprob,&Eloc_fields,&Bloc_fields,&Jloc_fields,&data[probeSize-1]);
                    
                    //! here we fill the probe data!!!
                    data[0]=Eloc_fields.x;
                    data[1]=Eloc_fields.y;
                    data[2]=Eloc_fields.z;
                    data[3]=Bloc_fields.x;
                    data[4]=Bloc_fields.y;
                    data[5]=Bloc_fields.z;
                    data[6]=Jloc_fields.x;
                    data[7]=Jloc_fields.y;
                    data[8]=Jloc_fields.z;
                    
                    H5Dwrite(did, H5T_NATIVE_DOUBLE, sidPart, sid, pid, &data[0]);
                } else {
                    H5Dwrite(did, H5T_NATIVE_DOUBLE, sidNull, sid, pid, &data[0]);
                }
                
                H5Pclose( pid );
                
            }
            H5Sclose(sidNull);
            H5Sclose(sidPart);
            H5Sclose(sid);
            H5Dclose(did);
        }
    }
    if (fileId) H5Fflush(fileId, H5F_SCOPE_GLOBAL );
}
