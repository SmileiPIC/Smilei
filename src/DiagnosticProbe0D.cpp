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

DiagnosticProbe0D::~DiagnosticProbe0D() {
}

void DiagnosticProbe0D::close() {
    H5Fclose(fileId);
}

DiagnosticProbe0D::DiagnosticProbe0D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi) :
    smpi_(smpi),
    probeSize(7)
{
    // Management of global IO file
    // All probe in a single file, a dataset per probe
    ostringstream file_name("");
    file_name<<"Probes0D.h5";

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

    hsize_t dims[2] = {0, probeSize};
    hsize_t max_dims[2] = {H5S_UNLIMITED, probeSize};
    hid_t file_space = H5Screate_simple(2, dims, max_dims);

    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    hsize_t chunk_dims[2] = {1, probeSize};
    H5Pset_chunk(plist, 2, chunk_dims);

    hsize_t dimsPos = diagParams->ps_0d_coord.size();
    DEBUG("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << dimsPos);
    hid_t dataspace_id = H5Screate_simple(1, &dimsPos, NULL);

    probeParticles.initialize(diagParams->ps_0d_coord[0].size(), diagParams->ps_0d_coord.size());
    DEBUG("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << diagParams->ps_0d_coord[0].size());
    
    double* tmp = new double[diagParams->ps_0d_coord.size()];
    for(unsigned int count=0; count!=diagParams->ps_0d_coord[0].size(); ++count) {
        int found=smpi_->getRank();
        for(unsigned int iDim=0; iDim!=diagParams->ps_0d_coord.size(); ++iDim) {
            if(smpi_->getDomainLocalMin(iDim)>diagParams->ps_0d_coord[iDim][count] || smpi_->getDomainLocalMax(iDim)<=diagParams->ps_0d_coord[iDim][count]) {
                found=-1;
            }
        }
        probeId.push_back(found);

        for(unsigned int iDim=0; iDim!=diagParams->ps_0d_coord.size(); ++iDim) {
            probeParticles.position(iDim,count)=diagParams->ps_0d_coord[iDim][count];
        }

        hid_t probeDataset_id = H5Dcreate(fileId, probeName(count).c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        hid_t attribute_id = H5Acreate2 (probeDataset_id, "Position", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        for (unsigned int iDim=0; iDim!=diagParams->ps_0d_coord.size(); ++iDim)
            tmp[iDim] = probeParticles.position(iDim,count)/(2.0*M_PI);
        H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, tmp);
        H5Aclose(attribute_id);
        H5Dclose(probeDataset_id);
    }
    delete [] tmp;
    H5Sclose(dataspace_id);

    // Create 1 dataset per probe

    H5Pclose(plist);
    H5Sclose(file_space);
}

string DiagnosticProbe0D::probeName(int p) {
    ostringstream prob_name("");
    prob_name << "p"<< setfill('0') << setw(4) << p;
    return prob_name.str();
}

void DiagnosticProbe0D::run(int timestep, ElectroMagn* EMfields, Interpolator* interp) {
    hsize_t dims[2];
    dims[0] = 1;
    dims[1] = probeSize;
    hid_t  partMemSpace = H5Screate_simple(2, dims, NULL);
    hsize_t nulldims[2];
    nulldims[0] = 0;
    nulldims[1] = 0;
    hid_t  partMemSpaceNull = H5Screate_simple(2, nulldims, NULL);

    vector<double> data(probeSize);

    for (int count=0; count <probeParticles.size(); count++) {
        if (probeId[count]==smpi_->getRank())
            (*interp)(EMfields,probeParticles,count,&Eloc_fields,&Bloc_fields);

        // All rank open all probes dataset
        hid_t dataset_id = H5Dopen2(fileId, probeName(count).c_str(), H5P_DEFAULT);
        hid_t file_space = H5Dget_space(dataset_id);

        // Get dataset existing dims
        hsize_t dimsO[2];
        H5Sget_simple_extent_dims(file_space, dimsO, NULL);
        H5Sclose(file_space);

        // Increment dataset size
        hsize_t dims[2];
        dims[0] = dimsO[0]+1;
        dims[1] = dimsO[1];
        H5Dset_extent(dataset_id, dims);
        //
        file_space = H5Dget_space(dataset_id);
        hsize_t start[2];
        hsize_t count2[2];
        if  (probeId[count]==smpi_->getRank()) {
            count2[0] = 1;
            count2[1] = probeSize;
        } else {
            count2[0] = 0;
            count2[1] = 0;
        }
        start[0] = dimsO[0];
        start[1] = 0;
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count2, NULL);

        //! here we fill the probe data!!!
        data[0]=timestep;
        data[1]=Eloc_fields.x;
        data[2]=Eloc_fields.y;
        data[3]=Eloc_fields.z;
        data[4]=Bloc_fields.x;
        data[5]=Bloc_fields.y;
        data[6]=Bloc_fields.z;

        hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);

        if  (probeId[count]==smpi_->getRank()) {
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, partMemSpace, file_space, write_plist,&data[0]);
        } else {
            // Write 0 data
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, partMemSpaceNull, file_space, write_plist,&data[0]);
        }

        H5Pclose( write_plist );

        H5Sclose(file_space);

        H5Dclose(dataset_id);
    }
    H5Fflush(fileId, H5F_SCOPE_GLOBAL );

    H5Sclose(partMemSpaceNull);
    H5Sclose(partMemSpace);

}
