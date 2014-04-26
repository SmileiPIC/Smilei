#include "DiagnosticPhase2D.h"

#include <sstream>
#include <iomanip>

#include "SmileiMPI.h"

using namespace std;


void DiagnosticPhase2D::writeData(unsigned int timestep, std::string species_name, SmileiMPI* smpi_) {
	
	Field2D my_data_sum;
	
	my_data_sum.allocateDims(my_data.dims());
	MPI_Reduce(my_data.data_,my_data_sum.data_,my_data.globalDims_,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	if (smpi_->isMaster()) {
		ostringstream name("");
		name << species_name << setfill('0') << setw(8) << timestep;
		
		hsize_t dims[2]={my_data.dims()[0],my_data.dims()[1]};
		hid_t sid = H5Screate_simple (2, dims, NULL);	
		hid_t did = H5Dcreate (groupID, name.str().c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, my_data_sum.data_);
		H5Dclose (did);	
		H5Sclose(sid);
	}
	
}