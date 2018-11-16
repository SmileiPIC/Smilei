
#include "PXR_Solver3D_GPSTD.h"

#include "ElectroMagn.h"
#include "Field3D.h"
#include "interface.h"

PXR_Solver3D_GPSTD::PXR_Solver3D_GPSTD(Params &params)
: Solver3D(params)
{
}

PXR_Solver3D_GPSTD::~PXR_Solver3D_GPSTD()
{
}

void PXR_Solver3D_GPSTD::coupling( Params &params, ElectroMagn* EMfields )
{
#ifdef _PICSAR
    int cdim=3;

    int n0,n1,n2;
    int ov0,ov1,ov2;
    // unable to convert unsigned int to an iso_c_binding supported type 

    n0=(int) (1 +  params.n_space[0]*params.global_factor[0]);
    n1=(int) (1 +  params.n_space[1]*params.global_factor[1]);
    n2=(int) (1 +  params.n_space[2]*params.global_factor[2]);

    ov0=(int) params.oversize[0];
    ov1=(int) params.oversize[1];
    ov2=(int) params.oversize[2];

    params.norderx = params.norder[0];
    params.nordery = params.norder[1];
    params.norderz = params.norder[2];

    Field3D* Ex3D_pxr = static_cast<Field3D*>( EMfields->Ex_pxr);
    Field3D* Ey3D_pxr = static_cast<Field3D*>( EMfields->Ey_pxr);
    Field3D* Ez3D_pxr = static_cast<Field3D*>( EMfields->Ez_pxr);
    Field3D* Bx3D_pxr = static_cast<Field3D*>( EMfields->Bx_pxr);
    Field3D* By3D_pxr = static_cast<Field3D*>( EMfields->By_pxr);
    Field3D* Bz3D_pxr = static_cast<Field3D*>( EMfields->Bz_pxr);
    Field3D* Jx3D_pxr = static_cast<Field3D*>( EMfields->Jx_pxr);
    Field3D* Jy3D_pxr = static_cast<Field3D*>( EMfields->Jy_pxr);
    Field3D* Jz3D_pxr = static_cast<Field3D*>( EMfields->Jz_pxr);
    Field3D* rho3D_pxr = static_cast<Field3D*>( EMfields->rho_pxr);
    Field3D* rhoold3D_pxr = static_cast<Field3D*>( EMfields->rhoold_pxr);

    //call of extern init routine (defined in picsar)
    picsar::init_params_picsar(&n0,&n1,&n2,
                       &params.cell_length[0],&params.cell_length[1],&params.cell_length[2],&params.timestep,
                       &ov0,&ov1,&ov2,
                       &params.norderx,&params.nordery,&params.norderz,
                       &params.is_spectral,
                       &(Ex3D_pxr->data_[0]),
                       &(Ey3D_pxr->data_[0]),
                       &(Ez3D_pxr->data_[0]),
                       &(Bx3D_pxr->data_[0]),
                       &(By3D_pxr->data_[0]),
                       &(Bz3D_pxr->data_[0]),
                       &(Jx3D_pxr->data_[0]),
                       &(Jy3D_pxr->data_[0]),
                       &(Jz3D_pxr->data_[0]),
                       &(rho3D_pxr->data_[0]),
                       &(rhoold3D_pxr->data_[0]),&cdim);
#else
    ERROR( "Smilei not linked with picsar" );
#endif


}

void PXR_Solver3D_GPSTD::operator() ( ElectroMagn* fields )
{
    duplicate_field_into_pxr( fields );

#ifdef _PICSAR
    picsar::push_psatd_ebfield_();
#else
    ERROR( "Smilei not linked with picsar" );
#endif

    duplicate_field_into_smilei( fields );

}

