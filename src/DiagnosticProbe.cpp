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

DiagnosticProbe::~DiagnosticProbe() {
}

void DiagnosticProbe::close() {
    if (fileId>0) {
        H5Fclose(fileId);
    }
}

void DiagnosticProbe::open() {
    ostringstream file_name("");
    file_name << "Probes" << dimProbe-2 << "D.h5";
    
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
}

string DiagnosticProbe::probeName(int p) {
    ostringstream prob_name("");
    prob_name << "p" << dimProbe-2 << "_" << setfill('0') << setw(4) << p;
    return prob_name.str();
}

void DiagnosticProbe::addProbe(unsigned int np, vector<double> partPos, vector<unsigned int> vecNumber) {
    
    DEBUG(dimProbe << " " << vecNumber.size());
    vector<hsize_t> dims(dimProbe);
    vector<hsize_t> max_dims(dimProbe);
    vector<hsize_t> chunk_dims(dimProbe);
    
    dims[0]=0;
    max_dims[0]=H5S_UNLIMITED;
    chunk_dims[0]=1;
    
    for (unsigned int i=0; i<vecNumber.size(); i++) {
        dims[i+1]=vecNumber[i];
        max_dims[i+1]=vecNumber[i];
        chunk_dims[i+1]=1;
    }
    dims.back()=probeSize;
    max_dims.back()=probeSize;
    chunk_dims.back()=probeSize;
    
    hid_t file_space = H5Screate_simple(dimProbe, &dims[0], &max_dims[0]);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    H5Pset_chunk(plist, dimProbe, &chunk_dims[0]);
    
    hid_t probeDataset_id = H5Dcreate(fileId, probeName(np).c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(file_space);
    
    unsigned int vecNumberProd=1;
    for (unsigned int i=0; i<vecNumber.size(); i++) {
        vecNumberProd*=vecNumber[i];
    }    
    unsigned int ndim=partPos.size()/vecNumberProd;
    
    vector<hsize_t> dimsPos(1+vecNumber.size());
    for (unsigned int i=0; i<vecNumber.size(); i++) {
        dimsPos[i]=vecNumber[i];
    }
    dimsPos[vecNumber.size()]=ndim;
    
    hid_t dataspace_id = H5Screate_simple(dimsPos.size(), &dimsPos[0], NULL);
    
    hid_t attribute_id = H5Acreate2 (probeDataset_id, "position", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &partPos[0]);
    H5Aclose(attribute_id);
    H5Sclose(dataspace_id);
    
    hsize_t dims0D = 1;
    hid_t sid = H5Screate_simple(1, &dims0D, NULL);	
    hid_t aid = H5Acreate(probeDataset_id, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, H5T_NATIVE_UINT, &every[np]);
    H5Sclose(sid);
    H5Aclose(aid);
    
    H5Dclose(probeDataset_id);
    
}

void DiagnosticProbe::runAll(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp) {
    for (unsigned int i=0; i<every.size(); i++) {
        if (every[i] && timestep % every[i] == 0) {
            run(i, EMfields, interp);
        }
    }
}

void DiagnosticProbe::run(unsigned int np, ElectroMagn* EMfields, Interpolator* interp) {
        
    vector<hsize_t> dims(dimProbe);
    vector<hsize_t> nulldims(dimProbe);
    for (unsigned int i=0; i<dimProbe-1; i++) {
        dims[i]=1;
        nulldims[i]=0;
    }
    dims.back()=probeSize;
    nulldims.back()=0;
    
    hid_t  partMemSpace = H5Screate_simple(dimProbe, &dims[0], NULL);
    hid_t  partMemSpaceNull = H5Screate_simple(dimProbe, &nulldims[0], NULL);
    
    vector<double> data(probeSize);
    
    unsigned int nprob=probeParticles[np].size();
    
    hid_t dataset_id = H5Dopen2(fileId, probeName(np).c_str(), H5P_DEFAULT);
    
    // All rank open all probes dataset
    hid_t file_space = H5Dget_space(dataset_id);
    
    // Get dataset existing dims
    
    vector<hsize_t> dimsO(dimProbe);
    
    H5Sget_simple_extent_dims(file_space, &dimsO[0], NULL);
    
    // Increment dataset size
    vector<hsize_t>  newDims=dimsO;
    newDims[0]++;

    H5Dset_extent(dataset_id, &newDims[0]);
    
    file_space = H5Dget_space(dataset_id);
    
    for (unsigned int iprob=0; iprob <nprob; iprob++) {
        if (probeId[np][iprob]==smpi_->getRank())
            (*interp)(EMfields,probeParticles[np],iprob,&Eloc_fields,&Bloc_fields);
        
        vector<hsize_t>  count2(dimProbe);
        if  (probeId[np][iprob]==smpi_->getRank()) {
            for (unsigned int i=0; i< dimProbe-1; i++) {
                count2[i]=1;
            }
            count2.back()=probeSize;
        } else {
            for (unsigned int i=0; i< dimProbe; i++) {
                count2[i]=0;
            }
        }
        
        vector<hsize_t> start(dimProbe);
        start[0]=dimsO[0];
    
        for (unsigned int i=1; i<dimProbe-1; i++) {
            start[i]=iprob % dimsO[i];
        }
        start.back()=0;
        
        
//        DEBUG("<><><><><><><><><><><><><><><> " << setw(4) << iprob);
//        for (unsigned int i=0; i<dimProbe; i++) {
//            DEBUG(dimProbe-2 << "D dimsO " << setw(4) << dimsO[i] << ": start " << setw(4)  << start[i]);
//        }
        
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, &start[0], NULL, &count2[0], NULL);
        
        //! here we fill the probe data!!!
        data[0]=Eloc_fields.x;
        data[1]=Eloc_fields.y;
        data[2]=Eloc_fields.z;
        data[3]=Bloc_fields.x;
        data[4]=Bloc_fields.y;
        data[5]=Bloc_fields.z;
        
        hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
        
        if  (probeId[np][iprob]==smpi_->getRank()) {
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
