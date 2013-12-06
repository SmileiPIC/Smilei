#include "prob0D.h"

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"

#include <string>

using namespace std;

prob0D::prob0D(PicParams* params, SmileiMPI* smpi,unsigned int num_prob):

n_prob(num_prob),
ps_coor(num_prob,vector<double>(params->nDim_field)),
weights_primal(num_prob,vector<double>(pow(2,params->nDim_field)*params->nDim_field)),
weights_dual(num_prob,vector<double>(pow(2,params->nDim_field)*params->nDim_field))
{
	here.resize(n_prob,false);
	params_=params;
	smpi_=smpi;
	vector<vector<double> > my_p(2,vector<double>(1));
	my_p[0][0]=0;
	my_p[1][0]=12;
	set_p_coor(0,my_p[0]);
	set_p_coor(1,my_p[1]);
	set_proc(); //debugged :)
	set_weights();
	
}

void prob0D::set_proc(){
    
    for(unsigned int p=0;p!=n_prob;++p){
        bool inside=true;
        for(unsigned int iDim=0; iDim!=params_->nDim_field;++iDim){
            if(smpi_->getDomainLocalMin(iDim)>get_ps_coor(p)[iDim]||smpi_->getDomainLocalMax(iDim)<get_ps_coor(p)[iDim]) inside=false;
        }
        here[p]=inside;
    }
    
}

// WARNING!!!Right now it is implemented for 1D!!
void prob0D::set_weights(){
   }

