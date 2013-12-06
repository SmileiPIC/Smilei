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
    DEBUG("-----------------------------------------------------------------------------------------------");
    n_prob=ps_coord[0].size();
    ps_coor=ps_coord;
    weights_primal.resize(n_prob);
    weights_dual.resize(n_prob);
    for (unsigned int i=0;i<n_prob;i++) {
        weights_primal[i].resize(params->interpolation_order*params->nDim_field);
        weights_dual[i].resize(params->interpolation_order*params->nDim_field);
    }

    DEBUG("----------------------------------------------------------------------------------------------- " << n_prob);
	here.resize(n_prob);
	params_=params;
	smpi_=smpi;
    DEBUG("----------------------------------------------------------------------------------------------- " << n_prob);
	set_proc();
    DEBUG("----------------------------------------------------------------------------------------------- " << n_prob);
	
}

void Probe0D::set_proc(){
    
    for(unsigned int p=0;p!=n_prob;++p){
        bool inside=true;
        for(unsigned int iDim=0; iDim!=ps_coor.size();++iDim){
            DEBUG(n_prob << " " << ps_coor[iDim].size() << " " << p << " " << iDim << " " << ps_coor[iDim][p]);

            if(smpi_->getDomainLocalMin(iDim)>ps_coor[iDim][p] || smpi_->getDomainLocalMax(iDim)<ps_coor[iDim][p]) {
                inside=false;    
            }
        }
        here[p]=inside;
    }
    
}

// WARNING!!!Right now it is implemented for 1D!!
void Probe0D::set_weights(){
}


void Probe0D::run(int timestep, ElectroMagn* EMfields, Interpolator* interp){
    for(unsigned int p=0;p!=n_prob;++p){
        DEBUG(smpi_->getRank() <<  " HERE " << timestep << " " << p<< " " << here[p]);
    }
}
