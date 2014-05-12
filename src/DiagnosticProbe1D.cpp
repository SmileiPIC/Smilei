#include "DiagnosticProbe1D.h"

#include <iomanip>
#include <string>
#include <iomanip>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"

using namespace std;

DiagnosticProbe1D::~DiagnosticProbe1D() {
}

void DiagnosticProbe1D::close() {
    if (fileId>0) {
        H5Fclose(fileId);
    }
}

DiagnosticProbe1D::DiagnosticProbe1D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi) :
smpi_(smpi),
probeSize(6),
fileId(0),
probeParticles(diagParams->probe1DStruc.size()),
probeId(diagParams->probe1DStruc.size()),
every(diagParams->probe1DStruc.size())
{
    if (diagParams->probe1DStruc.size() == 0) return;
    // Management of global IO file
    // All probe in a single file, a dataset per probe
    ostringstream file_name("");
    file_name<<"Probes1D.h5";

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    fileId = H5Fcreate( file_name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

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

    for (unsigned int np=0; np<diagParams->probe1DStruc.size(); np++) {
        
        every[np]=diagParams->probe1DStruc[np].every;
        unsigned int nprob=diagParams->probe1DStruc[np].number;
        
        hsize_t dims[3] = {0, probeSize, nprob};
        hsize_t max_dims[3] = {H5S_UNLIMITED, probeSize, nprob};
        hid_t file_space = H5Screate_simple(3, dims, max_dims);

        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        hsize_t chunk_dims[3] = {1, probeSize, 1};
        H5Pset_chunk(plist, 3, chunk_dims);
        

        unsigned int ndim=params->nDim_particle;

        
        probeParticles[np].initialize(diagParams->probe1DStruc[np].number, ndim);
        probeId[np].resize(diagParams->probe1DStruc[np].number);

        vector<double> partPos(ndim*nprob);

        for(unsigned int count=0; count!=nprob; ++count) {
            int found=smpi_->getRank();
            for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
                if (diagParams->probe1DStruc[np].number>1) {
                    partPos[iDim+count*ndim]=diagParams->probe1DStruc[np].posStart[iDim]+count*(diagParams->probe1DStruc[np].posEnd[iDim]-diagParams->probe1DStruc[np].posStart[iDim])/(diagParams->probe1DStruc[np].number-1);
                } else {
                    partPos[iDim+count*ndim]=0.5*(diagParams->probe1DStruc[np].posStart[iDim]+diagParams->probe1DStruc[np].posEnd[iDim]);
                }
                if(smpi_->getDomainLocalMin(iDim) >  partPos[iDim+count*ndim] || smpi_->getDomainLocalMax(iDim) <= partPos[iDim+count*ndim]) {
                    found=-1;
                }
                probeParticles[np].position(iDim,count)=partPos[iDim+count*ndim];
                partPos[iDim+count*ndim]/=2*M_PI;
            }
            probeId[np][count] = found;
        }
        
        //! write probe positions \todo check with 2D the row major order
        
        hid_t probeDataset_id = H5Dcreate(fileId, probeName(np).c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Pclose(plist);
        H5Sclose(file_space);
        
        hsize_t dimsPos[2] = {ndim, nprob};
        
        hid_t dataspace_id = H5Screate_simple(2, dimsPos, NULL);
        
        hid_t attribute_id = H5Acreate2 (probeDataset_id, "position", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &partPos[0]);
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);

        
        hsize_t dims1D[1] = {1};
        hid_t sid = H5Screate_simple(1, dims1D, NULL);	
        hid_t aid = H5Acreate(probeDataset_id, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, H5T_NATIVE_UINT, &diagParams->probe1DStruc[np].every);
        H5Sclose(sid);
        H5Aclose(aid);
        
        H5Dclose(probeDataset_id);        
        
    }
}

string DiagnosticProbe1D::probeName(int p) {
    ostringstream prob_name("");
    prob_name << "p"<< setfill('0') << setw(4) << p;
    return prob_name.str();
}

void DiagnosticProbe1D::run(int timestep, unsigned int np, ElectroMagn* EMfields, Interpolator* interp) {
    
    hsize_t dims[3] = {1, probeSize, 1};
    hid_t  partMemSpace = H5Screate_simple(3, dims, NULL);
    hsize_t nulldims[3] = {0, 0, 0};
    hid_t  partMemSpaceNull = H5Screate_simple(3, nulldims, NULL);

    vector<double> data(probeSize);
    
    unsigned int nprob=probeParticles[np].size();

    hid_t dataset_id = H5Dopen2(fileId, probeName(np).c_str(), H5P_DEFAULT);

    // All rank open all probes dataset
    hid_t file_space = H5Dget_space(dataset_id);

    // Get dataset existing dims
    hsize_t dimsO[3];
    H5Sget_simple_extent_dims(file_space, dimsO, NULL);
    
    // Increment dataset size
    hsize_t newDims[3] = { dimsO[0]+1, dimsO[1], dimsO[2]};
    H5Dset_extent(dataset_id, newDims);
    
    file_space = H5Dget_space(dataset_id);
    
    
    
    
//    hid_t attribute_idTime = H5Aopen(dataset_id, "Timesteps", H5P_DEFAULT);
    
    
    
    
    
    
    
    
    
    
    
    
    
    for (int count=0; count <nprob; count++) {
        if (probeId[np][count]==smpi_->getRank())
            (*interp)(EMfields,probeParticles[np],count,&Eloc_fields,&Bloc_fields);

        hsize_t count2[3];
        if  (probeId[np][count]==smpi_->getRank()) {
            count2[0] = 1;
            count2[1] = probeSize;
            count2[2] = 1;
        } else {
            count2[0] = 0;
            count2[1] = 0;
            count2[2] = 0;
        }
        hsize_t start[3] = { dimsO[0], 0, count};
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count2, NULL);

        //! here we fill the probe data!!!
        data[0]=Eloc_fields.x;
        data[1]=Eloc_fields.y;
        data[2]=Eloc_fields.z;
        data[3]=Bloc_fields.x;
        data[4]=Bloc_fields.y;
        data[5]=Bloc_fields.z;

        hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);

        if  (probeId[np][count]==smpi_->getRank()) {
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, partMemSpace, file_space, write_plist,&data[0]);
        } else {
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, partMemSpaceNull, file_space, write_plist,&data[0]);
        }

        H5Pclose( write_plist );

    }
    H5Sclose(file_space);
    
    H5Dclose(dataset_id);

    H5Fflush(fileId, H5F_SCOPE_GLOBAL );

    H5Sclose(partMemSpaceNull);
    H5Sclose(partMemSpace);
    
    
    
    
}
