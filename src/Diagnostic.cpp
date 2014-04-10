#include "Diagnostic.h"

#include <string>

#include <hdf5.h>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "Interpolator.h"
#include "DiagnosticScalar.h"

using namespace std;

Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi, Interpolator *interp) :
    diagScal(params, smpi),
	probe0D(params, diagparams, smpi),
	interp_(interp)
{
    everyScalar = diagparams->scalar_every;
    everyMap = diagparams->map_every;
    everyProbe0D = diagparams->probe0d_every;
	
	if (diagparams->phase1D.size()>0) {
		
		ostringstream file_name("");
		file_name<<"Phase1D.h5";
		
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

		
		for (unsigned int i=0; i<diagparams->phase1D.size();i++) {
			DiagnosticPhase1D *p=new DiagnosticPhase1D(params, diagparams, smpi,i);
			phase1DVec.push_back(p);
		}
	}
}

void Diagnostic::closeAll () {
    probe0D.close();
	
	H5Fclose(fileId);

}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies) {
    if (everyScalar && timestep % everyScalar == 0) {
        diagScal.run(timestep, EMfields, vecSpecies);
    }
    if (everyProbe0D && timestep % everyProbe0D == 0) {
        probe0D.run(timestep, EMfields, interp_);
    }
	
	
	for (unsigned int i=0; i<phase1DVec.size();i++) {
		if (phase1DVec[i]->every && timestep % phase1DVec[i]->every==0) {
			phase1DVec[i]->run(timestep, vecSpecies);
		}
	}
	
	
}

