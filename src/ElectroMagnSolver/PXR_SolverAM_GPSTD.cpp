
#include "PXR_SolverAM_GPSTD.h"

#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "interface.h"

PXR_SolverAM_GPSTD::PXR_SolverAM_GPSTD( Params &params )
    : SolverAM( params )
{
}

PXR_SolverAM_GPSTD::~PXR_SolverAM_GPSTD()
{
}

void PXR_SolverAM_GPSTD::coupling( Params &params, ElectroMagn *EMfields )
{
    int cdim=3;
    
    int nl, nr;
    int ovl, ovr;
    // unable to convert unsigned int to an iso_c_binding supported type
    
    std::vector<unsigned int> n_space(params.n_space);
    if (params.uncoupled_grids)
        n_space = params.n_space_domain;
    
    nl=( int ) (0 + n_space[0]);
    nr=( int ) (0 + n_space[1]);
    
    ovl=( int ) params.oversize[0];
    ovr=( int ) params.oversize[1];
    
    cField2D* El;
    cField2D* Er;
    cField2D* Et;
    cField2D* Bl;
    cField2D* Br;
    cField2D* Bt;
    cField2D* Jl;
    cField2D* Jr;
    cField2D* Jt;
    cField2D* rho;
    cField2D* rhoold;

    double pxr_dl = params.cell_length[0];
    double pxr_dr = params.cell_length[1];

    for ( int imode=0 ; imode<params.nmodes ; imode++ ) {

        El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
        Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
        Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
        Jl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jl_[imode];
        Jr = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jr_[imode];
        Jt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jt_[imode];
        rho = ( static_cast<ElectroMagnAM *>( EMfields ) )->rho_AM_[imode];
        rhoold = ( static_cast<ElectroMagnAM *>( EMfields ) )->rho_old_AM_[imode];

#ifdef _PICSAR    
        //call of extern init routine (defined in picsar)
        picsar::init_params_picsar_AM( &nr, &nl, &nmode, &imode, 
                                    &pxr_dr, &pxr_dl, &params.timestep,
                                    &ovr, &ovl,
                                    &params.norder[1], &params.norder[0],
                                    &params.is_spectral,
                                    &( El->cdata_[0] ),
                                    &( Er->cdata_[0] ),
                                    &( Et->cdata_[0] ),
                                    &( Bl->cdata_[0] ),
                                    &( Br->cdata_[0] ),
                                    &( Bt->cdata_[0] ),
                                    &( Jl->cdata_[0] ),
                                    &( Jr->cdata_[0] ),
                                    &( Jt->cdata_[0] ),
                                    &( rho->cdata_[0] ),
                                    &( rhoold->cdata_[0] ), &cdim );
#else
    ERROR( "Smilei not linked with picsar" );
#endif
    }
                                
                                
}

void PXR_SolverAM_GPSTD::operator()( ElectroMagn *fields )
{
    //duplicate_field_into_pxr( fields );
    
#ifdef _PICSAR
    picsar::push_psatd_ebfield_();
#else
    ERROR( "Smilei not linked with picsar" );
#endif
    
    //duplicate_field_into_smilei( fields );
    
}

