#include "DiagnosticPhase.h"

#include <sstream>
#include <iomanip>

#include "PyTools.h"

#include "SmileiMPI.h"

using namespace std;


DiagnosticPhase::~DiagnosticPhase(){
    if(dataId) H5Dclose(dataId);
};

DiagnosticPhase::DiagnosticPhase(Params &params, unsigned int n_phase) :
every(0),
dataId(0)
{
    
    if (!PyTools::extract("every",every,"DiagPhase",n_phase)) {
        every=params.global_every;
    }
        
    vector<double> time_range(2,0.);
    if (!PyTools::extract("time_range",time_range,"DiagPhase",n_phase)) {
        time_range[0]=0.;
        time_range[1]=params.sim_time;
    }
    tmin = time_range[0];
    tmax = time_range[1];

    
    PyTools::extract("species",my_species,"DiagPhase",n_phase);
    
    vector<double> min_max_n;
    
    PyTools::extract("first",min_max_n,"DiagPhase",n_phase);
    if (min_max_n.size()!=3) {
        ERROR("DiagPhase "<< n_phase << " first must have 3 components: [min, max, num]" );
    }
    if (min_max_n[2] != (int) min_max_n[2]) {
        ERROR("DiagPhase "<< n_phase << " first: third component must be integer" );
    }

    firstmin = min_max_n[0];
    firstmax = min_max_n[1];
    firstnum = (unsigned int) min_max_n[2];

    PyTools::extract("second",min_max_n,"DiagPhase",n_phase);
    if (min_max_n.size()!=3) {
        ERROR("DiagPhase "<< n_phase << " second must have 3 components: [min, max, num]" );
    }
    if (min_max_n[2] != (int) min_max_n[2]) {
        ERROR("DiagPhase "<< n_phase << " second : third component must be integer" );
    }

    secondmin = min_max_n[0];
    secondmax = min_max_n[1];
    secondnum = (unsigned int) min_max_n[2];
    
	if (every==0) ERROR("every cannot be zero");
    
    my_data.allocateDims(firstnum,secondnum);

}

void DiagnosticPhase::writeData() {
	
    MPI_Reduce(dataId>0?MPI_IN_PLACE:my_data.data_, my_data.data_,my_data.globalDims_,MPI_DOUBLE,  MPI_SUM, 0, MPI_COMM_WORLD);

    if (dataId>0) {
        hid_t sid = H5Dget_space(dataId);
        hsize_t dims[3];
        H5Sget_simple_extent_dims(sid, dims, NULL);
        H5Sclose(sid);
        
        hsize_t start[3]={dims[0],0,0};
        hsize_t count[3]={1,dims[1],dims[2]};

        // Increment dataset size
        dims[0]++;
        H5Dset_extent(dataId, dims);
        sid = H5Dget_space(dataId);
        
        hid_t sidChunk = H5Screate_simple(3,count,NULL);
        
        H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, count, NULL);
        
        H5Dwrite(dataId, H5T_NATIVE_DOUBLE, sidChunk, sid, H5P_DEFAULT, my_data.data_);
        
		H5Sclose(sidChunk);
		H5Sclose(sid);
	}
    
	//! we want to clean the Field back after ending for the next time    
    my_data.put_to(0.0);
	
}

