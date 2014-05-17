#include "DiagnosticPhase2D.h"

#include <sstream>
#include <iomanip>

#include "SmileiMPI.h"

using namespace std;

DiagnosticPhase2D::DiagnosticPhase2D(phaseStructure phaseStruct): DiagnosticPhase(phaseStruct) {
}

void DiagnosticPhase2D::writeAttributes(hid_t gid) {

    const vector<hsize_t> dimsPos(2,2);
    hid_t sid = H5Screate_simple(2, &dimsPos[0], NULL);

    vector<double> tmp(dimsPos[0]*dimsPos[1]);

    hid_t aid = H5Acreate (gid, "extents", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
    
    tmp[0]=firstmin;
    tmp[1]=secondmin;
    tmp[2]=firstmax;
    tmp[3]=secondmax;

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &tmp[0]);
    H5Aclose(aid);
    H5Sclose(sid);
    
}

void DiagnosticPhase2D::writeData(unsigned int timestep, hid_t gid) {
	
	Field2D my_data_sum;
	
	my_data_sum.allocateDims(my_data.dims());
	MPI_Reduce(my_data.data_,my_data_sum.data_,my_data.globalDims_,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	if (gid>0) {
		ostringstream name("");
		name << "t" << setfill('0') << setw(8) << timestep;
		
		hsize_t dims[2]={my_data.dims()[0],my_data.dims()[1]};
		hid_t sid = H5Screate_simple (2, dims, NULL);	
		hid_t did = H5Dcreate (gid, name.str().c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, my_data_sum.data_);
		H5Dclose (did);	
		H5Sclose(sid);
	}
    
	//! we want to clean the Field back after ending for the next time
	for (unsigned int i=0; i<my_data.globalDims_; i++) {
		my_data.data_[i]=0.0;
	}
	
}

