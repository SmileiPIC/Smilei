
#include "PXR_Solver2D_GPSTD.h"

#include "ElectroMagn.h"
#include "Field2D.h"
#include "interface.h"

PXR_Solver2D_GPSTD::PXR_Solver2D_GPSTD( Params &params )
    : Solver2D( params )
{
}

PXR_Solver2D_GPSTD::~PXR_Solver2D_GPSTD()
{
}

void PXR_Solver2D_GPSTD::coupling( Params &params, ElectroMagn *EMfields, bool full_domain )
{
#ifdef _PICSAR
    int cdim=2;
    int n0, n1, n2;
    int ov0, ov1, ov2;
    // unable to convert unsigned int to an iso_c_binding supported type
    
    std::vector<unsigned int> n_space(params.n_space);
    if (full_domain)
        n_space = params.n_space_global;
    else if (params.uncoupled_grids)
        n_space = params.n_space_region;
    
    n0=(int) (0 +  n_space[0]);
    n1=(int) (0 +  n_space[1]);
    
    n2=0;
    if (params.uncoupled_grids) {
        ov0=( int ) params.region_oversize[0];
        ov1=( int ) params.region_oversize[1];
    }
    else {
        ov0=( int ) params.oversize[0];
        ov1=( int ) params.oversize[1];
    }
    ov2=0;
    double dzz = std::numeric_limits<double>::infinity() ;
    params.norderx = params.norder[0];
    params.nordery = params.norder[1];
    params.norderz = 2;
    
    Field2D* Ex2D_pxr = static_cast<Field2D*>( EMfields->Ex_);
    Field2D* Ey2D_pxr = static_cast<Field2D*>( EMfields->Ey_);
    Field2D* Ez2D_pxr = static_cast<Field2D*>( EMfields->Ez_);
    Field2D* Bx2D_pxr = static_cast<Field2D*>( EMfields->Bx_);
    Field2D* By2D_pxr = static_cast<Field2D*>( EMfields->By_);
    Field2D* Bz2D_pxr = static_cast<Field2D*>( EMfields->Bz_);
    Field2D* Jx2D_pxr = static_cast<Field2D*>( EMfields->Jx_);
    Field2D* Jy2D_pxr = static_cast<Field2D*>( EMfields->Jy_);
    Field2D* Jz2D_pxr = static_cast<Field2D*>( EMfields->Jz_);
    Field2D* rho2D_pxr = static_cast<Field2D*>( EMfields->rho_);
    Field2D* rhoold2D_pxr = static_cast<Field2D*>( EMfields->rhoold_);
    
    double pxr_dx = -params.cell_length[0];
    double pxr_dy = -params.cell_length[1];
    
    picsar::init_params_picsar( &n0, &n2, &n1,
                                &pxr_dx, &dzz, &pxr_dy, &params.timestep,
                                &ov0, &ov2, &ov1,
                                &params.norderx, &params.norderz, &params.nordery,
                                &params.is_spectral,
                                &( Ex2D_pxr->data_[0] ),
                                &( Ez2D_pxr->data_[0] ),
                                &( Ey2D_pxr->data_[0] ),
                                &( Bx2D_pxr->data_[0] ),
                                &( Bz2D_pxr->data_[0] ),
                                &( By2D_pxr->data_[0] ),
                                &( Jx2D_pxr->data_[0] ),
                                &( Jz2D_pxr->data_[0] ),
                                &( Jy2D_pxr->data_[0] ),
                                &( rho2D_pxr->data_[0] ),
                                &( rhoold2D_pxr->data_[0] ), &cdim );
#else
    ERROR( "Smilei not linked with picsar, use make config=picsar" );
#endif
                                
}

void PXR_Solver2D_GPSTD::operator()( ElectroMagn *fields )
{
    //duplicate_field_into_pxr( fields );
    
#ifdef _PICSAR
    picsar::push_psatd_ebfield_();
#else
    ERROR( "Smilei not linked with picsar, use make config=picsar" );
#endif
    
    //duplicate_field_into_smilei( fields );
    
}

