#include "DiagnosticProbe0D.h"

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"
#include "hdf5.h"
#include <iomanip>
#include <string>
#include <iomanip>

using namespace std;


DiagnosticProbe0D::~DiagnosticProbe0D() {
}

void DiagnosticProbe0D::close() {
    H5Fclose(probeFileGlobal_id);
}


DiagnosticProbe0D::DiagnosticProbe0D(PicParams* params, SmileiMPI* smpi, vector<vector<double> > ps_coord) :
smpi_(smpi)
{
	// Management of global IO file
	// All probe in a single file, a dataset per probe
	ostringstream file_name(""); file_name<<"probe_global.h5";
	
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
	probeFileGlobal_id = H5Fcreate( file_name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	//>>>>>>>>>>>
	hid_t group_id = H5Gcreate2(probeFileGlobal_id, "Positions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//<<<<<<<<<<<

	
	for(unsigned int count=0;count!=ps_coord[0].size();++count){
		bool found=true;
		for(unsigned int iDim=0; iDim!=ps_coord.size();++iDim){
			DEBUG(smpi_->getRank() << " " << iDim << " " <<  count << " " << smpi_->getDomainLocalMin(iDim) << " " << ps_coord[iDim][count]  << " " << smpi_->getDomainLocalMax(iDim));
			if(smpi_->getDomainLocalMin(iDim)>ps_coord[iDim][count] || smpi_->getDomainLocalMax(iDim)<=ps_coord[iDim][count]) {
				found=false;
			}
		}	
		
		
		probeId.push_back(found);
		
		Particle *part=new Particle(params->nDim_field);
		for(unsigned int iDim=0;iDim!=ps_coord.size();++iDim){
			part->position(iDim)=ps_coord[iDim][count];
		}
		probeParticles.push_back(part);
		
		
		hsize_t dims[1] = {ps_coord.size()};
		hid_t file_space = H5Screate_simple(1, dims, NULL);
		hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
		
		hid_t probeDataset_id = H5Dcreate(group_id, ("/Positions/"+probeName(count)).c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
		H5Dwrite(probeDataset_id, H5T_NATIVE_DOUBLE, H5P_DEFAULT, file_space, H5P_DEFAULT,part->position());

		H5Dclose(probeDataset_id);
		H5Pclose(plist);
		H5Sclose(file_space);
		

	}
	H5Gclose(group_id);

	//    DEBUG("--------------------------------------------------------------");
	
	Eloc_fields.resize(probeParticles.size());
	Bloc_fields.resize(probeParticles.size());
	
	
	//    DEBUG("--------------------------------------------------------------");
	
	
	

	
	// 7 = timestep + Exyz + Bxyz
	probeSize = 7;
	
	hsize_t dims[2] = {0, probeSize};
	hsize_t max_dims[2] = {H5S_UNLIMITED, probeSize};
	hid_t file_space = H5Screate_simple(2, dims, max_dims);
	
	hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_layout(plist, H5D_CHUNKED);
	hsize_t chunk_dims[2] = {1, probeSize};
	H5Pset_chunk(plist, 2, chunk_dims);
	// Create 1 dataset per probe
	for(unsigned int count=0;count!=ps_coord[0].size();++count){
		hid_t probeDataset_id = H5Dcreate(probeFileGlobal_id, probeName(count).c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
		H5Dclose(probeDataset_id);
	}
	
	H5Pclose(plist);
	H5Sclose(file_space);
	
}

string DiagnosticProbe0D::probeName(int p) {
	ostringstream prob_name(""); prob_name<<"/probe_"<< setfill('0') << setw(4) << p;
	return prob_name.str();
}

void DiagnosticProbe0D::run(int timestep, ElectroMagn* EMfields, Interpolator* interp){
	hsize_t dims[2];
	dims[0] = 1;
	dims[1] = probeSize;
	hid_t  partMemSpace = H5Screate_simple(2, dims, NULL);
	hsize_t nulldims[2];
	nulldims[0] = 0;
	nulldims[1] = 0;
	hid_t  partMemSpaceNull = H5Screate_simple(2, nulldims, NULL);
	
	for (unsigned int count=0; count <probeParticles.size(); count++) {
		if (probeId[count])
			(*interp)(EMfields,probeParticles[count],&Eloc_fields[count],&Bloc_fields[count]);
		
		// All rank open all probes dataset
		hid_t dataset_id = H5Dopen2(probeFileGlobal_id, probeName(count).c_str(), H5P_DEFAULT);
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
		if  (probeId[count]) {
			count2[0] = 1;
			count2[1] = probeSize;
		} else {
			count2[0] = 0;
			count2[1] = 0;		
		}
        start[0] = dimsO[0];
        start[1] = 0;
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count2, NULL);
		
		vector<double> data(probeSize);
        data[0]=timestep;
        data[1]=Eloc_fields[count].x;
        data[2]=Eloc_fields[count].y;
        data[3]=Eloc_fields[count].z;
        data[4]=Bloc_fields[count].x;
        data[5]=Bloc_fields[count].y;
        data[6]=Bloc_fields[count].z;		
		
		hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
		
		if  (probeId[count]) {
			H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, partMemSpace, file_space, write_plist,&data[0]);
		} else {
			// Write 0 data
			H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, partMemSpaceNull, file_space, write_plist,&data[0]);
		}
		
        H5Pclose( write_plist );
		
		H5Sclose(file_space);
		
		H5Dclose(dataset_id);
    }
    H5Fflush( probeFileGlobal_id, H5F_SCOPE_GLOBAL );
    
    H5Sclose(partMemSpaceNull);
    H5Sclose(partMemSpace);
	
}
