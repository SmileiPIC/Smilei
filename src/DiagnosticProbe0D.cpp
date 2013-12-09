#include "DiagnosticProbe0D.h"

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"

#include <string>

using namespace std;


DiagnosticProbe0D::~DiagnosticProbe0D() {
//	H5Fclose(probe_file_id_);
}

DiagnosticProbe0D::DiagnosticProbe0D(PicParams* params, SmileiMPI* smpi, vector<vector<double> > ps_coord) :
smpi_(smpi)
{
	// Management of global IO file
	
//	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
//	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
//	
//    probe_file_id_ = H5Fcreate("Probe0D.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//	H5Pclose(plist_id);
 

	bool found=true;
	for(unsigned int p=0;p!=ps_coord[0].size();++p){
        for(unsigned int iDim=0; iDim!=ps_coord.size();++iDim){
            if(smpi_->getDomainLocalMin(iDim)>ps_coord[iDim][p] || smpi_->getDomainLocalMax(iDim)<=ps_coord[iDim][p]) 
				found=false;
		}			
		if (found) {
			Particle *part=new Particle(params->nDim_field);
			for(unsigned int iDim=0;iDim!=ps_coord.size();++iDim){
				part->position(iDim)=ps_coord[iDim][p];
			}
			newParticle.push_back(part);
			ostringstream name_t;
			name_t.str(""); name_t << p;
			DEBUG("created group " << name_t.str());
//			hid_t group_id = H5Gcreate2(probe_file_id_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, plist_id);
//			H5Gclose(group_id);
//			H5Fflush( probe_file_id_, H5F_SCOPE_GLOBAL );
        }
    }


    Eloc_fields.resize(newParticle.size());
    Bloc_fields.resize(newParticle.size());
    
}


void DiagnosticProbe0D::run(int timestep, ElectroMagn* EMfields, Interpolator* interp){
	for (unsigned int count=0; count <newParticle.size(); count++) {
		(*interp)(EMfields,newParticle[count],&Eloc_fields[count],&Bloc_fields[count]);
		DEBUG(smpi_->getRank() << " " << count << " E " << Eloc_fields[count].x << " " << Eloc_fields[count].y << " " << Eloc_fields[count].z);
		DEBUG(smpi_->getRank() << " " << count << " B " << Bloc_fields[count].x << " " << Bloc_fields[count].y << " " << Bloc_fields[count].z);
	}
}
