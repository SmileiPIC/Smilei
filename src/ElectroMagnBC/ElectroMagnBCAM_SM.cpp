#include "ElectroMagnBCAM_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
//#include "Field2D.h"
#include "cField2D.h"
#include "Tools.h"
#include "Laser.h"
#include <complex>
#include "dcomplex.h"

using namespace std;

ElectroMagnBCAM_SM::ElectroMagnBCAM_SM( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBCAM( params, patch, i_boundary )
{
    // conversion factor from degree to radian
    conv_deg2rad = M_PI/180.0;
    
    //Number of modes
    Nmode = params.nmodes;
    
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
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    // Xmin boundary
    //double theta  = 0.0*conv_deg2rad; //0.0;
    double factor = 1.0/( 1.0 + dt_ov_d[0] );
    Alpha_Xmin    = 2.0*factor;
    Beta_Xmin     = - ( 1-dt_ov_d[0] )*factor;
    Gamma_Xmin    = 4.0*factor;
    Delta_Xmin    = - dt_ov_d[1]*factor;
    Epsilon_Xmin  = Icpx*factor*dt ;
    // Xmax boundary
    //theta         = M_PI;
    factor        = 1.0/( 1.0 + dt_ov_d[0] );
    Alpha_Xmax    = 2.0*factor;
    Beta_Xmax     = - ( 1.0 -dt_ov_d[0] )*factor;
    Gamma_Xmax    = 4.0*factor;
    Delta_Xmax    = - dt_ov_d[1]*factor;
    Epsilon_Xmax  = - Icpx*factor*dt;
   
}


void ElectroMagnBCAM_SM::save_fields( Field *my_field, Patch *patch )
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

void ElectroMagnBCAM_SM::disableExternalFields()
{
    Bl_val.resize( 0 );
    Br_val.resize( 0 );
    Bt_val.resize( 0 );
}



// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCAM_SM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    //ERROR( "Test SM" )
    // Loop on imode
    for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
        // Static cast of the fields
        cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
        cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
        cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
        bool isYmin = ( static_cast<ElectroMagnAM *>( EMfields ) )->isYmin;
        int     j_glob = ( static_cast<ElectroMagnAM *>( EMfields ) )->j_glob_;
        
        if( i_boundary_ == 0 && patch->isXmin() ) {
            // for Br^(d,p)
            vector<double> yp( 1 );
            for( unsigned int j=3*isYmin ; j<n_p[1] ; j++ ) {
            
                std::complex<double> byW = 0.;
                yp[0] = patch->getDomainLocalMin( 1 ) +( (double)j - (double)EMfields->oversize[1] )*d[1];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    if (vecLaser[ilaser]->spacetime.size() > 2){
                        byW +=          vecLaser[ilaser]->getAmplitudecomplexN(yp, time_dual, 0, 0, 2*imode);
                    } else {                        
                        if( imode==1 ) {
                            byW +=          vecLaser[ilaser]->getAmplitude0( yp, time_dual, 1+2*j, 0 )
                                            + Icpx * vecLaser[ilaser]->getAmplitude1( yp, time_dual, 1+2*j, 0 );
                        }
                    }
                }
                
                //x= Xmin
                unsigned int i=0;
                ( *Br )( i, j ) = Alpha_Xmin   * ( *Et )( i, j )
                                  +              Beta_Xmin    * ( *Br )( i+1, j )
                                  +              Gamma_Xmin   * byW;
                +              Delta_Xmin   *( ( *Bl )( i, j+1 )- ( *Bl )( i, j ) );
            }//j  ---end compute Br
            
            
            // for Bt^(d,d)
            vector<double> yd( 1 );
            for( unsigned int j=3*isYmin; j<n_d[1] ; j++ ) {
            
                std::complex<double> bzW = 0.;
                yd[0] = patch->getDomainLocalMin( 1 ) + ( (double)j - 0.5 - (double)EMfields->oversize[1] )*d[1];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    if (vecLaser[ilaser]->spacetime.size() > 2){
                        bzW +=          vecLaser[ilaser]->getAmplitudecomplexN(yd, time_dual, 0, 0, 2*imode+1);
                    } else {                        
                        if( imode==1 ) {
                            bzW +=          vecLaser[ilaser]->getAmplitude1( yd, time_dual, 2*j, 0 )
                                           - Icpx * vecLaser[ilaser]->getAmplitude0( yd, time_dual, 2*j, 0 );
                        }
                    }
                }
                //x=Xmin
                unsigned int i=0;
                ( *Bt )( i, j ) = -Alpha_Xmin * ( *Er )( i, j )
                                  +               Beta_Xmin  *( ( *Bt )( i+1, j ) )
                                  +               Gamma_Xmin * bzW
                                  +               Epsilon_Xmin *( double )imode/( ( j_glob+j-0.5 )*d[1] )*( *Bl )( i, j ) ;
                
            }//j  ---end compute Bt
            //Redo condition on axis for Bt because it was modified
            if( isYmin && imode != 1 ) {
                ( *Bt )( 0, 2 ) = -( *Bt )( 0, 3 );
            }
            if( isYmin && imode == 1 ) {
                ( *Bt )( 0, 2 )= -2.*Icpx*( *Br )( 0, 2 )-( *Bt )( 0, 3 );
            }
        } else if( i_boundary_ == 1 && patch->isXmax() ) {
            // for Br^(d,p)
            vector<double> yp( 1 );
            for( unsigned int j=3*isYmin ; j<n_p[1] ; j++ ) {
            
                std::complex<double> byE = 0.;
                // Lasers
                if( imode==1 ) {
                    yp[0] = patch->getDomainLocalMin( 1 ) + ( (double)j - (double)EMfields->oversize[1] )*d[1];
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        byE +=          vecLaser[ilaser]->getAmplitude0( yp, time_dual, 1+2*j, 0 )
                                        + Icpx * vecLaser[ilaser]->getAmplitude1( yp, time_dual, 1+2*j, 0 );
                    }
                }
                unsigned int i= n_p[0];
                ( *Br )( i, j ) = - Alpha_Xmax   * ( *Et )( i-1, j )
                                  +                   Beta_Xmax    * ( *Br )( i-1, j )
                                  +                   Gamma_Xmax   * byE
                                  +                   Delta_Xmax   * ( ( *Bl )( i-1, j+1 )- ( *Bl )( i-1, j ) ); // Check x-index
                
            }//j  ---end compute Br
            
            // for Bt^(d,d)
            vector<double> yd( 1 );
            for( unsigned int j=3*isYmin ; j<n_d[1] ; j++ ) {
            
                std::complex<double> bzE = 0.;
                if( imode==1 ) {
                    yd[0] = patch->getDomainLocalMin( 1 ) + ( (double)j - 0.5  - (double)EMfields->oversize[1] )*d[1];
                    // Lasers
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        bzE +=         vecLaser[ilaser]->getAmplitude1( yd, time_dual, 2*j, 0 )
                                       -Icpx * vecLaser[ilaser]->getAmplitude0( yd, time_dual, 2*j, 0 );
                    }
                }
                unsigned int i= n_p[0];
                ( *Bt )( i, j ) = Alpha_Xmax * ( *Er )( i-1, j )
                                  +                    Beta_Xmax  * ( *Bt )( i-1, j )
                                  +                    Gamma_Xmax * bzE
                                  +					  Epsilon_Xmax * ( double )imode /( ( j_glob+j-0.5 )*d[1] )* ( *Bl )( i-1, j )	;
                
            }//j  ---end compute Bt
            //Redo condition on axis for Bt because it was modified
            if( isYmin && imode != 1 ) {
                ( *Bt )( n_p[0], 2 ) = -( *Bt )( n_p[0], 3 );
            }
            if( isYmin && imode == 1 ) {
                ( *Bt )( n_p[0], 2 )= -2.*Icpx*( *Br )( n_p[0], 2 )-( *Bt )( n_p[0], 3 );
            }
        }
    }
}
