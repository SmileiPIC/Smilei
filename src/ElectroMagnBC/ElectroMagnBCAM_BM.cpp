#include "ElectroMagnBCAM_BM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "cField2D.h"
#include "Tools.h"
#include "Laser.h"
#include <complex>
#include "dcomplex.h"
using namespace std;

ElectroMagnBCAM_BM::ElectroMagnBCAM_BM( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBCAM( params, patch, i_boundary )
{
    //Number of modes
    Nmode=params.nmodes;
    if( i_boundary_ == 0 && patch->isXmin() ) {
        // BCs at the x-border min
        Bl_val.resize( n_d[1], 0. ); // dual in the y-direction
        Br_val.resize( n_p[1], 0. ); // primal in the y-direction
        Bt_val.resize( n_d[1], 0. ); // dual in the y-direction
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
        // BCs at the x-border max
        Bl_val.resize( n_d[1], 0. );
        Br_val.resize( n_p[1], 0. );
        Bt_val.resize( n_d[1], 0. );
    } else if( i_boundary_ == 2 && patch->isYmin() ) {
        // BCs in the y-border min
        Bl_val.resize( n_p[0], 0. ); // primal in the x-direction
        Br_val.resize( n_d[0], 0. ); // dual in the x-direction
        Bt_val.resize( n_d[0], 0. ); // dual in the x-direction
    } else if( i_boundary_ == 3 && patch->isYmax() ) {
        // BCs in the y-border max
        Bl_val.resize( n_p[0], 0. );
        Br_val.resize( n_d[0], 0. );
        Bt_val.resize( n_d[0], 0. );
    }
    
    // -----------------------------------------------------
    // Parameters for the Buneman boundary conditions
    // -----------------------------------------------------
    
    // Rmax boundary
    double cosphi ;
    double Kx, Kr ;
    double one_ov_rlocal;
    if (!params.uncoupled_grids)
        one_ov_rlocal = 1./( params.grid_length[1]+params.oversize[1]*d[1] ); // BM conditions on rmax are written at the last primal r position.
    else
        one_ov_rlocal = 1./( params.grid_length[1]+params.region_oversize[1]*d[1] ); // BM conditions on rmax are written at the last primal r position.
    //std::cout<< "rlocal " <<params.grid_length[1]+params.oversize[1]*d[1]<< "one_ov"<< one_ov_rlocal*10<< std::endl;
    //std::cout<< "grid length " <<params.grid_length[1]<< "   oversize*d[1]  "<< params.oversize[1]*d[1]<< std::endl;
    Kx =  params.EM_BCs_k[3][0];
    Kr = -params.EM_BCs_k[3][1]; // We're only dealing with the Rmax boundary here. The minus sign is the specular reflexion of the given k on the rmax boundary since users are supposed to provide the injection k.
    cosphi = Kr / sqrt( Kx*Kx + Kr*Kr ) ;
    CB_BM  = cosphi/( 1. + cosphi ); // Theta is always taken equal to zero.
    CE_BM  = 1.0 - CB_BM;
    
    //Coeffs for Bl
    double factor= 1./( 1. + dt_ov_d[1] );
    Alpha_Bl_Rmax    = ( 1. - dt_ov_d[1] )*factor;
    Beta_Bl_Rmax     = CE_BM * dt * one_ov_rlocal * factor;
    Gamma_Bl_Rmax    = CB_BM * dt_ov_d[0]*factor;
    
    //Coeffs for Bt
    factor          = 1. / ( 1 + dt_ov_d[1] + 0.5*CB_BM*dt*one_ov_rlocal );
    Alpha_Bt_Rmax   = ( -1 + dt_ov_d[1] - 0.5*CB_BM*dt*one_ov_rlocal ) * factor;
    Beta_Bt_Rmax    = ( 1 - dt_ov_d[1] - 0.5*CB_BM*dt*one_ov_rlocal ) * factor;
    Gamma_Bt_Rmax   = ( 1 + dt_ov_d[1] - 0.5*CB_BM*dt*one_ov_rlocal ) * factor;
    Epsilon_Bt_Rmax = dt * one_ov_rlocal * factor;
    Delta_Bt_Rmax   = dt_ov_d[0]*factor;
   
}

void ElectroMagnBCAM_BM::save_fields( Field *my_field, Patch *patch )
{
    cField2D *field2D=static_cast<cField2D *>( my_field );
    
    if( i_boundary_ == 0 && patch->isXmin() ) {
        if( field2D->name=="Bl" ) {
            for( unsigned int j=0; j<n_d[1]; j++ ) {
                Bl_val[j]=( *field2D )( 0, j );
            }
        }
        
        if( field2D->name=="Br" ) {
            for( unsigned int j=0; j<n_p[1]; j++ ) {
                Br_val[j]=( *field2D )( 0, j );
            }
        }
        
        if( field2D->name=="Bt" ) {
            for( unsigned int j=0; j<n_d[1]; j++ ) {
                Bt_val[j]=( *field2D )( 0, j );
            }
        }
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
        if( field2D->name=="Bl" ) {
            for( unsigned int j=0; j<n_d[1]; j++ ) {
                Bl_val[j]=( *field2D )( n_p[0]-1, j );
            }
        }
        
        if( field2D->name=="Br" ) {
            for( unsigned int j=0; j<n_p[1]; j++ ) {
                Br_val[j]=( *field2D )( n_d[0]-1, j );
            }
        }
        
        if( field2D->name=="Bt" ) {
            for( unsigned int j=0; j<n_d[1]; j++ ) {
                Bt_val[j]=( *field2D )( n_d[0]-1, j );
            }
        }
    } else if( i_boundary_ == 2 && patch->isYmin() ) {
        if( field2D->name=="Bl" ) {
            for( unsigned int i=0; i<n_p[0]; i++ ) {
                Bl_val[i]=( *field2D )( i, 0 );
            }
        }
        
        if( field2D->name=="Br" ) {
            for( unsigned int i=0; i<n_d[0]; i++ ) {
                Br_val[i]=( *field2D )( i, 0 );
            }
        }
        
        if( field2D->name=="Bt" ) {
            for( unsigned int i=0; i<n_d[0]; i++ ) {
                Bt_val[i]=( *field2D )( i, 0 );
            }
        }
    } else if( i_boundary_ == 3 && patch->isYmax() ) {
        if( field2D->name=="Bl" ) {
            for( unsigned int i=0; i<n_p[0]; i++ ) {
                Bl_val[i]=( *field2D )( i, n_d[1]-1 );
            }
        }
        
        if( field2D->name=="Br" ) {
            for( unsigned int i=0; i<n_d[0]; i++ ) {
                Br_val[i]=( *field2D )( i, n_p[1]-1 );
            }
        }
        
        if( field2D->name=="Bt" ) {
            for( unsigned int i=0; i<n_d[0]; i++ ) {
                Bt_val[i]=( *field2D )( i, n_d[1]-1 );
            }
        }
    }
}

void ElectroMagnBCAM_BM::disableExternalFields()
{
    Bl_val.resize( 0 );
    Br_val.resize( 0 );
    Bt_val.resize( 0 );
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCAM_BM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{

    //This condition can only be applied to Rmax
    
    // Loop on imode
    for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
    
        // Static cast of the fields
        cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
        cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
        cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
        cField2D *Bt_old = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];
        cField2D *Bl_old = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
        cField2D *Br_old = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
        
        if( i_boundary_ == 3 && patch->isYmax() ) {
        
            unsigned int j= n_d[1]-2;
            
            
            //   MESSAGE("JGLOB "<< patch->getCellStartingGlobalIndex(1)+j);
            //std::cout<<"come heree "<<patch->getCellStartingGlobalIndex(1)<<"  "<<j<<" \n " ;
            //std::cout<<"come here "<<n_p[1] <<" nr*d[1] "<<n_p[1]*d[1]<<" \n " ;
            // for Bl^(p,d)
            for( unsigned int i=0 ; i<n_p[0]-1; i++ ) {
                ( *Bl )( i, j+1 ) = ( *Bl_old )( i, j )
                                    -      Alpha_Bl_Rmax * ( ( *Bl )( i, j ) - ( *Bl_old )( i, j+1 ) )
                                    +      Gamma_Bl_Rmax * ( ( *Br )( i+1, j ) + ( *Br_old )( i+1, j ) - ( *Br )( i, j ) - ( *Br_old )( i, j ) )
                                    -      Beta_Bl_Rmax * Icpx * ( double )imode * ( ( *Er )( i, j+1 ) + ( *Er )( i, j ) )
                                    - 2. * Beta_Bl_Rmax * ( *Et )( i, j );
                //if (std::abs((*Bl)(i,j+1))>1.){
                //MESSAGE("BlBM");
                //MESSAGE(i);
                //MESSAGE(j+1);
                //MESSAGE((*Bl)(i,j+1));
                //}
            }//i  ---end Bl
            
            // for Bt^(d,d)
            j = n_d[1]-2;
            for( unsigned int i=1 ; i<n_p[0] ; i++ ) { //Undefined in i=0 and i=n_p[0]
                ( *Bt )( i, j+1 ) =     Alpha_Bt_Rmax * ( *Bt )( i, j )
                                        + Beta_Bt_Rmax  * ( *Bt_old )( i, j+1 )
                                        + Gamma_Bt_Rmax * ( *Bt_old )( i, j )
                                        - Icpx * ( double )imode * CB_BM * Epsilon_Bt_Rmax  * ( ( *Br )( i, j ) - ( *Br_old )( i, j ) )
                                        - CE_BM * Delta_Bt_Rmax * ( ( *Er )( i, j+1 )+( *Er )( i, j )-( *Er )( i-1, j+1 ) -( *Er )( i-1, j ) ) ;
                //if (std::abs((*Bt)(i,j+1))>1.){
                //    MESSAGE("BtMF");
                //    MESSAGE(i);
                //    MESSAGE(j+1);
                //    MESSAGE((*Bt)(i,j+1));
                //}
            }//i  ---end Bt
        }
    }
}

