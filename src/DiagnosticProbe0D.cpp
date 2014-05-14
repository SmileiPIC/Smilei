#include "DiagnosticProbe0D.h"

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

DiagnosticProbe0D::DiagnosticProbe0D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi) : DiagnosticProbe(params,diagParams,smpi) {
    every.resize(diagParams->probe0DStruc.size());
    probeParticles.resize(diagParams->probe0DStruc.size());
    probeId.resize(diagParams->probe0DStruc.size());
    
    if (diagParams->probe0DStruc.size() == 0) return;
    
    open("Probes0D.h5");
    
    for (unsigned int np=0; np<diagParams->probe0DStruc.size(); np++) {
        
        every[np]=diagParams->probe0DStruc[np].every;

        hsize_t dims[3] = {0, probeSize, 1};
        hsize_t max_dims[3] = {H5S_UNLIMITED, probeSize, 1};
        hid_t file_space = H5Screate_simple(3, dims, max_dims);
        
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        hsize_t chunk_dims[3] = {1, probeSize, 1};
        H5Pset_chunk(plist, 3, chunk_dims);
        
        
        unsigned int ndim=params->nDim_particle;
        
        probeParticles[np].initialize(1, ndim);
        probeId[np].resize(1);
        
        vector<double> partPos(ndim);
        
        int found=smpi->getRank();
        for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
            partPos[iDim]=diagParams->probe0DStruc[np].pos[iDim];
            probeParticles[np].position(iDim,0) = 2*M_PI*partPos[iDim];
            if(smpi->getDomainLocalMin(iDim) >  probeParticles[np].position(iDim,0) || smpi->getDomainLocalMax(iDim) <= probeParticles[np].position(iDim,0)) {
                found=-1;
            }
        }
        probeId[np][0] = found;
        
        //! write probe positions \todo check with 2D the row major order
        
        hid_t probeDataset_id = H5Dcreate(fileId, probeName(np).c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Pclose(plist);
        H5Sclose(file_space);
        
        hsize_t dimsPos[2] = {ndim, 1};
        
        hid_t dataspace_id = H5Screate_simple(2, dimsPos, NULL);
        
        hid_t attribute_id = H5Acreate2 (probeDataset_id, "position", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &partPos[0]);
        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);
        
        
        hsize_t dims0D[1] = {1};
        hid_t sid = H5Screate_simple(1, dims0D, NULL);	
        hid_t aid = H5Acreate(probeDataset_id, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, H5T_NATIVE_UINT, &diagParams->probe0DStruc[np].every);
        H5Sclose(sid);
        H5Aclose(aid);
        
        H5Dclose(probeDataset_id);        
        
    }
}
