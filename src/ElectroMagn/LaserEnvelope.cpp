        
#include "LaserEnvelope.h"

#include "Params.h"
#include "cField3D.h"
#include "Field3D.h"
#include "ElectroMagn.h"
using namespace std;

LaserEnvelope::LaserEnvelope( Params& params )
{
}


LaserEnvelope::LaserEnvelope( LaserEnvelope *envelope )
{
}


LaserEnvelope::~LaserEnvelope()
{
}


LaserEnvelope3D::LaserEnvelope3D( Params& params )
    : LaserEnvelope(params)
{
    std::vector<unsigned int>  dimPrim( params.nDim_field );
    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<params.nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = params.n_space[i]+1;
        // + Ghost domain
        dimPrim[i] += 2*params.oversize[i];
    }

    
    A_  = new cField3D( dimPrim );
    A0_ = new cField3D( dimPrim );

    // see Python !!!
    cField3D* A3D = static_cast<cField3D*>(A_);
    //(*A3D)(i,j,k) = 
}


LaserEnvelope3D::LaserEnvelope3D( LaserEnvelope *envelope )
    : LaserEnvelope(envelope)
{
    A_  = new cField3D( envelope->A_->dims_  );
    A0_ = new cField3D( envelope->A0_->dims_ );

    // see Python !!!
    cField3D* A3D = static_cast<cField3D*>(A_);
    //(*A3D)(i,j,k) = 
}


LaserEnvelope3D::~LaserEnvelope3D()
{
    delete A_;
    delete A0_;
}

void LaserEnvelope3D::compute(ElectroMagn* EMfields)
{
    //->rho_e- ???;
    cField3D* A3D = static_cast<cField3D*>(A_);
    cField3D* A03D = static_cast<cField3D*>(A0_);

    // find e_idx in all species  
    int e_idx = 0;
    Field3D* rho_e = static_cast<Field3D*>(EMfields->rho_s[e_idx]);
    
    for (unsigned int i=0 ; i <A_->dims_[0]; i++)
        for (unsigned int j=0 ; j < A_->dims_[1] ; j++)
            for (unsigned int k=0 ; k < A_->dims_[2]; k++)
                (*A3D)(i,j,k) = (*A03D)(i,j,k);
}
