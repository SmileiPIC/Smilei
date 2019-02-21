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

ElectroMagnBCAM_SM::ElectroMagnBCAM_SM( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC( params, patch, _min_max )
{
    // conversion factor from degree to radian
    conv_deg2rad = M_PI/180.0;
    
    // number of nodes of the primal and dual grid in the x-direction
    nl_p = params.n_space[0]+1+2*params.oversize[0];
    nl_d = nl_p+1;
    // number of nodes of the primal and dual grid in the y-direction
    nr_p = params.n_space[1]+1+2*params.oversize[1];
    nr_d = nr_p+1;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dl       = params.cell_length[0];
    dt_ov_dl = dt/dl;
    dl_ov_dt = 1.0/dt_ov_dl;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dr       = params.cell_length[1];
    dt_ov_dr = dt/dr;
    dr_ov_dt = 1.0/dt_ov_dr;
    
    //Number of modes
    Nmode= params.nmodes;
    
    if( min_max == 0 && patch->isXmin() ) {
        // BCs at the x-border min
        Bl_val.resize( nr_d, 0. ); // dual in the y-direction
        Br_val.resize( nr_p, 0. ); // primal in the y-direction
        Bt_val.resize( nr_d, 0. ); // dual in the y-direction
    } else if( min_max == 1 && patch->isXmax() ) {
        // BCs at the x-border max
        Bl_val.resize( nr_d, 0. );
        Br_val.resize( nr_p, 0. );
        Bt_val.resize( nr_d, 0. );
    } else if( min_max == 2 && patch->isYmin() ) {
        // BCs in the y-border min
        Bl_val.resize( nl_p, 0. ); // primal in the x-direction
        Br_val.resize( nl_d, 0. ); // dual in the x-direction
        Bt_val.resize( nl_d, 0. ); // dual in the x-direction
    } else if( min_max == 3 && patch->isYmax() ) {
        // BCs in the y-border max
        Bl_val.resize( nl_p, 0. );
        Br_val.resize( nl_d, 0. );
        Bt_val.resize( nl_d, 0. );
    }
    
    
    
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
#ifdef _TODO_AM
#endif
    // Xmin boundary
    //double theta  = 0.0*conv_deg2rad; //0.0;
    double factor = 1.0/( 1.0 + dt_ov_dl );
    Alpha_SM_Xmin    = 2.0*factor;
    Beta_SM_Xmin     = - ( 1-dt_ov_dl )*factor;
    Gamma_SM_Xmin    = 4.0*factor;
    Delta_SM_Xmin    = - dt_ov_dr*factor;
    Epsilon_SM_Xmin  = Icpx*factor*dt ;
    // Xmax boundary
    //theta         = M_PI;
    factor        = 1.0/( 1.0 + dt_ov_dl );
    Alpha_SM_Xmax    = 2.0*factor;
    Beta_SM_Xmax     = - ( 1.0 -dt_ov_dl )*factor;
    Gamma_SM_Xmax    = 4.0*factor;
    Delta_SM_Xmax    = - dt_ov_dr*factor;
    Epsilon_SM_Xmax  = Icpx*factor*dt;
    
    
    
}


void ElectroMagnBCAM_SM::save_fields( Field *my_field, Patch *patch )
{
    cField2D *field2D=static_cast<cField2D *>( my_field );
    
    if( min_max == 0 && patch->isXmin() ) {
        if( field2D->name=="Bl" ) {
            for( unsigned int j=0; j<nr_d; j++ ) {
                Bl_val[j]=( *field2D )( 0, j );
            }
        }
        
        if( field2D->name=="Br" ) {
            for( unsigned int j=0; j<nr_p; j++ ) {
                Br_val[j]=( *field2D )( 0, j );
            }
        }
        
        if( field2D->name=="Bt" ) {
            for( unsigned int j=0; j<nr_d; j++ ) {
                Bt_val[j]=( *field2D )( 0, j );
            }
        }
    } else if( min_max == 1 && patch->isXmax() ) {
        if( field2D->name=="Bl" ) {
            for( unsigned int j=0; j<nr_d; j++ ) {
                Bl_val[j]=( *field2D )( nl_p-1, j );
            }
        }
        
        if( field2D->name=="Br" ) {
            for( unsigned int j=0; j<nr_p; j++ ) {
                Br_val[j]=( *field2D )( nl_d-1, j );
            }
        }
        
        if( field2D->name=="Bt" ) {
            for( unsigned int j=0; j<nr_d; j++ ) {
                Bt_val[j]=( *field2D )( nl_d-1, j );
            }
        }
    } else if( min_max == 2 && patch->isYmin() ) {
        if( field2D->name=="Bl" ) {
            for( unsigned int i=0; i<nl_p; i++ ) {
                Bl_val[i]=( *field2D )( i, 0 );
            }
        }
        
        if( field2D->name=="Br" ) {
            for( unsigned int i=0; i<nl_d; i++ ) {
                Br_val[i]=( *field2D )( i, 0 );
            }
        }
        
        if( field2D->name=="Bt" ) {
            for( unsigned int i=0; i<nl_d; i++ ) {
                Bt_val[i]=( *field2D )( i, 0 );
            }
        }
    } else if( min_max == 3 && patch->isYmax() ) {
        if( field2D->name=="Bl" ) {
            for( unsigned int i=0; i<nl_p; i++ ) {
                Bl_val[i]=( *field2D )( i, nr_d-1 );
            }
        }
        
        if( field2D->name=="Br" ) {
            for( unsigned int i=0; i<nl_d; i++ ) {
                Br_val[i]=( *field2D )( i, nr_p-1 );
            }
        }
        
        if( field2D->name=="Bt" ) {
            for( unsigned int i=0; i<nl_d; i++ ) {
                Bt_val[i]=( *field2D )( i, nr_d-1 );
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
    // Loop on imode
    for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
        // Static cast of the fields
        //cField2D* El = (static_cast<ElectroMagnAM*>(EMfields))->El_[imode];
        cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
        cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
        cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
        bool isYmin = ( static_cast<ElectroMagnAM *>( EMfields ) )->isYmin;
        //bool isYmax = (static_cast<ElectroMagnAM*>(EMfields))->isYmax;
        int     j_glob = ( static_cast<ElectroMagnAM *>( EMfields ) )->j_glob_;
        
        if( min_max == 0 && patch->isXmin() ) {
            //MESSAGE("Xmin");
            // for Br^(d,p)
            vector<double> yp( 1 );
            //yp[0] = patch->getDomainLocalMin(1) - EMfields->oversize[1]*dr;
            for( unsigned int j=2*isYmin ; j<nr_p-1 ; j++ ) {
            
                std::complex<double> byW = 0.;
                //yp[0] += dr;
                yp[0] = patch->getDomainLocalMin( 1 ) +( j - EMfields->oversize[1] )*dr;
                if( imode==1 ) {
                    // Lasers
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        byW +=          vecLaser[ilaser]->getAmplitude0( yp, time_dual, 1+2*j, 0 )
                                        + Icpx * vecLaser[ilaser]->getAmplitude1( yp, time_dual, 1+2*j, 0 );
                        //if (std::abs(byW)>1.0){
                        //MESSAGE("byW");
                        //MESSAGE(byW);
                        //MESSAGE("j");
                        //MESSAGE(j);}
                    }
                }
                
                //x= Xmin
                unsigned int i=0;
                ( *Br )( i, j ) = Alpha_SM_Xmin   * ( *Et )( i, j )
                                  +              Beta_SM_Xmin    * ( *Br )( i+1, j )
                                  +              Gamma_SM_Xmin   * byW;
                +              Delta_SM_Xmin   *( ( *Bl )( i, j+1 )- ( *Bl )( i, j ) );
                //        if (std::abs((*Br)(i,j))>1.){
                //            MESSAGE("BrAMSM");
                //            MESSAGE(i);
                //            MESSAGE(j);
                //            MESSAGE((*Br)(i,j));
                //        }
            }//j  ---end compute Br

            
            // for Bt^(d,d)
            vector<double> yd( 1 );
            for( unsigned int j=3*isYmin; j<nr_d-1 ; j++ ) {
            
                std::complex<double> bzW = 0.;
                yd[0] = patch->getDomainLocalMin( 1 ) + ( j - 0.5 - EMfields->oversize[1] )*dr;
                if( imode==1 ) {
                    // Lasers
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        bzW +=          vecLaser[ilaser]->getAmplitude1( yd, time_dual, 2*j, 0 )
                                        - Icpx * vecLaser[ilaser]->getAmplitude0( yd, time_dual, 2*j, 0 );
                        //if (std::abs(bzW)>1.0){
                        //	MESSAGE("bzW");
                        //	MESSAGE(bzW);
                        //	MESSAGE("j");
                        //	MESSAGE(j);}
                    }
                }
                //x=Xmin
                unsigned int i=0;
                ( *Bt )( i, j ) = -Alpha_SM_Xmin * ( *Er )( i, j )
                                  +               Beta_SM_Xmin  *( ( *Bt )( i+1, j ) )
                                  +               Gamma_SM_Xmin * bzW
                                  +               Epsilon_SM_Xmin *( double )imode/( ( j_glob+j+0.5 )*dr )*( *Bl )( i, j ) ;
                //if (std::abs((*Bt)(i,j))>1.){
                //MESSAGE("BtSM");
                //MESSAGE(i);
                //MESSAGE(j);
                //MESSAGE((*Bt)(i,j));
                //}
                
                //MESSAGE("EPS")
                //MESSAGE(Epsilon_SM_Xmin *(double)imode/((j_glob+j+0.5)*dr));
                
            }//j  ---end compute Bt
            if( isYmin && imode != 1 ) {
                ( *Bt )( 0, 2 ) = -( *Bt )( 0, 3 );
            }
        } else if( min_max == 1 && patch->isXmax() ) {
            //MESSAGE("Xmax");
            // for Br^(d,p)
            vector<double> yp( 1 );
            for( unsigned int j=2*isYmin ; j<nr_p-1 ; j++ ) {
            
                std::complex<double> byE = 0.;
                yp[0] = patch->getDomainLocalMin( 1 ) + ( j - EMfields->oversize[1] )*dr;
                // Lasers
                if( imode==1 ) {
                    // Lasers
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        byE +=          vecLaser[ilaser]->getAmplitude0( yp, time_dual, 1+2*j, 0 )
                                        + Icpx * vecLaser[ilaser]->getAmplitude1( yp, time_dual, 1+2*j, 0 );
                        //if (std::abs(byE)>1.0){
                        //MESSAGE("byE");
                        //MESSAGE(byE);
                        //MESSAGE("j");
                        //MESSAGE(j);}
                    }
                }
                unsigned int i= nl_d-1;
                ( *Br )( i, j ) = - Alpha_SM_Xmax   * ( *Et )( i-1, j )
                                  +                   Beta_SM_Xmax    * ( *Br )( i-1, j )
                                  +                   Gamma_SM_Xmax   * byE
                                  +                   Delta_SM_Xmax   * ( ( *Bl )( i-1, j+1 )- ( *Bl )( i-1, j ) ); // Check x-index
                //if (std::abs((*Br)(i,j))>1.){
                //MESSAGE("BrSMM");
                //MESSAGE(i);
                //MESSAGE(j);
                //MESSAGE((*Br)(i,j));
                //}
                
            }//j  ---end compute Br
            
            // for Bt^(d,d)
            vector<double> yd( 1 );
            for( unsigned int j=3*isYmin ; j<nr_d-1 ; j++ ) {
            
                std::complex<double> bzE = 0.;
                yd[0] = patch->getDomainLocalMin( 1 ) + ( j - 0.5  - EMfields->oversize[1] )*dr;
                // Lasers
                if( imode==1 ) {
                    // Lasers
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        bzE +=         vecLaser[ilaser]->getAmplitude1( yd, time_dual, 2*j, 0 )
                                       -Icpx * vecLaser[ilaser]->getAmplitude0( yd, time_dual, 2*j, 0 );
                        //if (std::abs(bzE)>1.0){
                        //MESSAGE("bzE");
                        //MESSAGE(bzE);
                        //MESSAGE("j");
                        //MESSAGE(j);}
                    }
                }
                unsigned int i= nl_d-1;
                ( *Bt )( i, j ) = Alpha_SM_Xmax * ( *Er )( i-1, j )
                                  +                    Beta_SM_Xmax  * ( *Bt )( i-1, j )
                                  +                    Gamma_SM_Xmax * bzE
                                  +					  Epsilon_SM_Xmax * ( double )imode /( ( j_glob+j+0.5 )*dr )* ( *Bl )( i-1, j )	;
                //if (std::abs((*Bt)(i,j))>1.){
                //MESSAGE("BtSMM");
                //MESSAGE(i);
                //MESSAGE(j);
                //MESSAGE((*Bt)(i,j));
                //}
                
            }//j  ---end compute Bt
            if( isYmin && imode != 1 ) {
                ( *Bt )( nl_d-1, 2 ) = -( *Bt )( nl_d-1, 3 );
            }
        }
    }
}
