#include "DiagnosticProbe0D.h"

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "hdf5.h"
#include <iomanip>
#include <string>

using namespace std;


DiagnosticProbe0D::~DiagnosticProbe0D() {
}

void DiagnosticProbe0D::close() {

    for(unsigned int i=0;i!=probeFile_id.size();++i){
        H5Fclose(probeFile_id[i]);
    }
    
}


DiagnosticProbe0D::DiagnosticProbe0D(PicParams* params, SmileiMPI* smpi, vector<vector<double> > ps_coord) :
smpi_(smpi)
{
	// Management of global IO file
    
    dims[0]=7;

     
	bool found=true;
    unsigned int nFound=0;
	for(unsigned int p=0;p!=ps_coord[0].size();++p){
        for(unsigned int iDim=0; iDim!=ps_coord.size();++iDim){
            DEBUG(smpi_->getRank() << " " << iDim << " " <<  p << " " << smpi_->getDomainLocalMin(iDim) << " " << ps_coord[iDim][p]  << " " << smpi_->getDomainLocalMax(iDim));
            if(smpi_->getDomainLocalMin(iDim)>ps_coord[iDim][p] || smpi_->getDomainLocalMax(iDim)<=ps_coord[iDim][p]) {
				found=false;
            }
		}			
		if (found) {
            nFound++;
            ostringstream file_name;
            file_name.str("");file_name<<"probe_"<<p<<".h5";
            hid_t file= H5Fcreate(file_name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT );
            probeFile_id.push_back(file);
            DEBUG("--------------------------------------------------------------");
			Particle *part=new Particle(params->nDim_field);
			for(unsigned int iDim=0;iDim!=ps_coord.size();++iDim){
				part->position(iDim)=ps_coord[iDim][p];
			}
			probeParticles.push_back(part);
        }
    }

//    DEBUG("--------------------------------------------------------------");


    Eloc_fields.resize(probeParticles.size());
    Bloc_fields.resize(probeParticles.size());

       
//    DEBUG("--------------------------------------------------------------");

}


void DiagnosticProbe0D::run(int timestep, ElectroMagn* EMfields, Interpolator* interp){
    for (unsigned int count=0; count <probeParticles.size(); count++) {
		(*interp)(EMfields,probeParticles[count],&Eloc_fields[count],&Bloc_fields[count]);
        
//        DEBUG("--------------------------------------------------------------");
            
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hid_t space_id = H5Screate_simple(1,dims,NULL);
        ostringstream name_t;
        name_t.str(""); name_t << "/T" << setw(8) << timestep << "_" << smpi_->getRank() << "_" << count;

        hid_t dset_id = H5Dcreate(probeFile_id[count], name_t.str().c_str(), H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        
        data[0]=timestep;
        data[1]=Eloc_fields[count].x;
        data[2]=Eloc_fields[count].y;
        data[3]=Eloc_fields[count].z;
        data[4]=Bloc_fields[count].x;
        data[5]=Bloc_fields[count].y;
        data[6]=Bloc_fields[count].z;
        
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, H5S_ALL , H5S_ALL, H5P_DEFAULT, data );
        
        H5Dclose(dset_id);
        
        
//		DEBUG(smpi_->getRank() << " " << count << " E " << Eloc_fields[count].x << " " << Eloc_fields[count].y << " " << Eloc_fields[count].z);
//		DEBUG(smpi_->getRank() << " " << count << " B " << Bloc_fields[count].x << " " << Bloc_fields[count].y << " " << Bloc_fields[count].z);
	}
}
