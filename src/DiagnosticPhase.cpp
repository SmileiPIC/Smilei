#include "DiagnosticPhase.h"

#include <sstream>
#include <iomanip>

#include "SmileiMPI.h"

using namespace std;

DiagnosticPhase::DiagnosticPhase(phaseStructure phaseStruct) :
my_species(phaseStruct.species)
{
	every=phaseStruct.every;
	if (every==0) ERROR("every cannot be zero");
}

void DiagnosticPhase::writeData(hid_t did) {
	
	Field2D my_data_sum;
	my_data_sum.allocateDims(my_data.dims());

	MPI_Reduce(my_data.data_,my_data_sum.data_,my_data.globalDims_,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	    
    if (did>0) {
        hid_t sid = H5Dget_space(did);
        hsize_t dims[3];
        H5Sget_simple_extent_dims(sid, dims, NULL);
        H5Sclose(sid);
        
        hsize_t start[3]={dims[0],0,0};
        hsize_t count[3]={1,dims[1],dims[2]};

        // Increment dataset size
        dims[0]++;
        H5Dset_extent(did, dims);
        sid = H5Dget_space(did);
        
        hid_t sidChunk = H5Screate_simple(3,count,NULL);
        
        H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, count, NULL);
        
        H5Dwrite(did, H5T_NATIVE_DOUBLE, sidChunk, sid, H5P_DEFAULT, my_data_sum.data_);
        
		H5Sclose(sidChunk);
		H5Sclose(sid);
	}
    
	//! we want to clean the Field back after ending for the next time    
    my_data.put_to(0.0);
	
}

