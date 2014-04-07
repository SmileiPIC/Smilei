#include "DiagnosticPhase1D.h"

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

DiagnosticPhase1D::~DiagnosticPhase1D() {
}

void DiagnosticPhase1D::close() {
    H5Fclose(fileId);
}

DiagnosticPhase1D::DiagnosticPhase1D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi) :
    smpi_(smpi)
{
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

	//! momentum min
	momentum_min.resize(params->n_species);
	momentum_max.resize(params->n_species);
	
	
	for (unsigned int i=0; i<params->n_species; i++) {
		momentum_min[i].resize(3);
		momentum_max[i].resize(3);
	}
	
	
	lorentz_factor_min.resize(params->n_species);
	lorentz_factor_max.resize(params->n_species);
	
}

void DiagnosticPhase1D::run(int timestep, std::vector<Species*>& vecSpecies) {
	
	for (unsigned int j=0; j < vecSpecies.size(); j++) {
		
		// initialize momentum min/max for phase space diags
		if (vecSpecies[j]->particles.size()>0) {
			for (unsigned int i=0; i<3; i++) {
				momentum_min[j][i] = momentum_max[j][i] = vecSpecies[j]->particles.momentum(i,0);
			}
			lorentz_factor_min[j]=lorentz_factor_max[j]=vecSpecies[j]->particles.lor_fac(0);
		}
		
		for (unsigned int ibin = 0 ; ibin < vecSpecies[j]->bmin.size() ; ibin++) {
			for (int iPart=vecSpecies[j]->bmin[ibin] ; iPart<vecSpecies[j]->bmax[ibin]; iPart++ ) {
				
				for (unsigned int i=0; i<3; i++) {
					if (vecSpecies[j]->particles.momentum(i,iPart) < momentum_min[j][i]) momentum_min[j][i] = vecSpecies[j]->particles.momentum(i,iPart);
					if (vecSpecies[j]->particles.momentum(i,iPart) > momentum_max[j][i]) momentum_max[j][i] = vecSpecies[j]->particles.momentum(i,iPart);					
				}
				if (vecSpecies[j]->particles.lor_fac(iPart) < lorentz_factor_min[j]) lorentz_factor_min[j] = vecSpecies[j]->particles.lor_fac(iPart);
				if (vecSpecies[j]->particles.lor_fac(iPart) > lorentz_factor_max[j]) lorentz_factor_max[j] = vecSpecies[j]->particles.lor_fac(iPart);
				
			}
		}
	}
	
	for (unsigned int j=0; j < vecSpecies.size(); j++) {
		for (unsigned int i=0; i<3; i++) {
			double tempmin=momentum_min[j][i];
			MPI_Allreduce(&tempmin,&momentum_min[j][i],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
			double tempmax=momentum_max[j][i];
			MPI_Allreduce(&tempmax,&momentum_max[j][i],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		}
		double gammamin=lorentz_factor_min[j];
		MPI_Allreduce(&gammamin,&lorentz_factor_min[j],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		double gammamax=lorentz_factor_max[j];
		MPI_Allreduce(&gammamax,&lorentz_factor_max[j],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

		if (smpi_->isMaster()) {
			for (unsigned int k=0; k <momentum_min[j].size(); k++) {				
				DEBUG(timestep << " spec. " << j << " min: " << momentum_min[j][k] << " max: "<< momentum_max[j][k]);
			}
			DEBUG(timestep << " spec. " << j << " Lmin: " << lorentz_factor_min[j] << " Lmax: "<< lorentz_factor_max[j]);
		}
	}

}
