#include "Probe0D.h"

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"

#include <string>

using namespace std;

Probe0D::Probe0D(PicParams* params, SmileiMPI* smpi, vector<vector<double> > ps_coord)
{
    n_probe=ps_coord[0].size();
    n_probe_loc=0;
    ps_coor=ps_coord;
	here.resize(n_probe);
	params_=params;
	smpi_=smpi;
    set_proc();
    newParticle.resize(n_probe_loc);
    Eloc_fields.resize(n_probe_loc);
    Bloc_fields.resize(n_probe_loc);
    
}



void Probe0D::set_proc(){
    
    for(unsigned int p=0;p!=n_probe;++p){
        bool inside=true;
        for(unsigned int iDim=0; iDim!=ps_coor.size();++iDim){
            if(smpi_->getDomainLocalMin(iDim)>ps_coor[iDim][p] || smpi_->getDomainLocalMax(iDim)<ps_coor[iDim][p]) {
                inside=false;    
            } else { 
                n_probe_loc+=1;
            }
        }
        here[p]=inside;
    }
    
}

void Probe0D::run(int timestep, ElectroMagn* EMfields, Interpolator* interp){
    unsigned int count=0;
    for(unsigned int p=0;p!=n_probe;++p){
        if(here[p]==true){
            newParticle[count]=new Particle(params_->nDim_field);
            for(unsigned int iDim=0;iDim!=params_->nDim_field;++iDim){
                newParticle[count]->position(iDim)=ps_coor[iDim][p];
                (*interp)(EMfields,newParticle[count],&Eloc_fields[count],&Bloc_fields[count]);
            }
            count+=1;
        }
    }
    
    
}
