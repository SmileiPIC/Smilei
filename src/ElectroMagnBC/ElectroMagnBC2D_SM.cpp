#include "ElectroMagnBC2D_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field2D.h"
#include "Tools.h"
#include "Laser.h"

using namespace std;

ElectroMagnBC2D_SM::ElectroMagnBC2D_SM( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC2D( params, patch, _min_max )
{
    if( min_max == 0 && patch->isXmin() ) {
        // BCs at the x-border min
        Bx_val.resize( ny_d, 0. ); // dual in the y-direction
        By_val.resize( ny_p, 0. ); // primal in the y-direction
        Bz_val.resize( ny_d, 0. ); // dual in the y-direction
    } else if( min_max == 1 && patch->isXmax() ) {
        // BCs at the x-border max
        Bx_val.resize( ny_d, 0. );
        By_val.resize( ny_p, 0. );
        Bz_val.resize( ny_d, 0. );
    } else if( min_max == 2 && patch->isYmin() ) {
        // BCs in the y-border min
        Bx_val.resize( nx_p, 0. ); // primal in the x-direction
        By_val.resize( nx_d, 0. ); // dual in the x-direction
        Bz_val.resize( nx_d, 0. ); // dual in the x-direction
    } else if( min_max == 3 && patch->isYmax() ) {
        // BCs in the y-border max
        Bx_val.resize( nx_p, 0. );
        By_val.resize( nx_d, 0. );
        Bz_val.resize( nx_d, 0. );
    }
    
    
    
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    double pyKx, pyKy; //, pyKz;
    double kx, ky; //, kz;
    double Knorm;
    double omega = 1. ;
    //! \todo (MG) Check optimal angle for Silver-Muller BCs
    
    // Xmin boundary
    pyKx = params.EM_BCs_k[0][0];
    pyKy = params.EM_BCs_k[0][1];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    
    double factor = 1.0 / ( kx + dt_ov_dx );
    Alpha_Xmin    = 2.0                     * factor;
    Beta_Xmin     = - ( kx - dt_ov_dx ) * factor;
    Gamma_Xmin    = 4.0 * kx        * factor;
    Delta_Xmin    = - ( ky + dt_ov_dy ) * factor;
    Epsilon_Xmin  = - ( ky - dt_ov_dy ) * factor;
    
    // Xmax boundary
    pyKx = params.EM_BCs_k[1][0];
    pyKy = params.EM_BCs_k[1][1];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    
    factor        = 1.0 / ( kx - dt_ov_dx );
    Alpha_Xmax    = 2.0                      * factor;
    Beta_Xmax     = - ( kx + dt_ov_dx )  * factor;
    Gamma_Xmax    = 4.0 * kx         * factor;
    Delta_Xmax    = - ( ky + dt_ov_dy )  * factor;
    Epsilon_Xmax  = - ( ky - dt_ov_dy )  * factor;
    
    // Ymin boundary
    pyKx = params.EM_BCs_k[2][0];
    pyKy = params.EM_BCs_k[2][1];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    
    factor = 1.0 / ( ky + dt_ov_dy );
    Alpha_Ymin    = 2.0                     * factor;
    Beta_Ymin     = - ( ky - dt_ov_dy ) * factor;
    Gamma_Ymin    = 4.0 * ky        * factor;
    Delta_Ymin    = - ( kx + dt_ov_dx ) * factor;
    Epsilon_Ymin  = - ( kx - dt_ov_dx ) * factor;
    
    // Ymax boundary
    pyKx = params.EM_BCs_k[3][0];
    pyKy = params.EM_BCs_k[3][1];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    
    factor = 1.0 / ( ky - dt_ov_dy );
    Alpha_Ymax    = 2.0                     * factor;
    Beta_Ymax     = - ( ky + dt_ov_dy ) * factor;
    Gamma_Ymin    = 4.0 * ky        * factor;
    Delta_Ymax    = - ( kx + dt_ov_dx ) * factor;
    Epsilon_Ymax  = - ( kx - dt_ov_dx ) * factor;
    
}


void ElectroMagnBC2D_SM::save_fields( Field *my_field, Patch *patch )
{
    Field2D *field2D=static_cast<Field2D *>( my_field );
    
    if( min_max == 0 && patch->isXmin() ) {
        if( field2D->name=="Bx" ) {
            for( unsigned int j=0; j<ny_d; j++ ) {
                Bx_val[j]=( *field2D )( 0, j );
            }
        }
        
        if( field2D->name=="By" ) {
            for( unsigned int j=0; j<ny_p; j++ ) {
                By_val[j]=( *field2D )( 0, j );
            }
        }
        
        if( field2D->name=="Bz" ) {
            for( unsigned int j=0; j<ny_d; j++ ) {
                Bz_val[j]=( *field2D )( 0, j );
            }
        }
    } else if( min_max == 1 && patch->isXmax() ) {
        if( field2D->name=="Bx" ) {
            for( unsigned int j=0; j<ny_d; j++ ) {
                Bx_val[j]=( *field2D )( nx_p-1, j );
            }
        }
        
        if( field2D->name=="By" ) {
            for( unsigned int j=0; j<ny_p; j++ ) {
                By_val[j]=( *field2D )( nx_d-1, j );
            }
        }
        
        if( field2D->name=="Bz" ) {
            for( unsigned int j=0; j<ny_d; j++ ) {
                Bz_val[j]=( *field2D )( nx_d-1, j );
            }
        }
    } else if( min_max == 2 && patch->isYmin() ) {
        if( field2D->name=="Bx" ) {
            for( unsigned int i=0; i<nx_p; i++ ) {
                Bx_val[i]=( *field2D )( i, 0 );
            }
        }
        
        if( field2D->name=="By" ) {
            for( unsigned int i=0; i<nx_d; i++ ) {
                By_val[i]=( *field2D )( i, 0 );
            }
        }
        
        if( field2D->name=="Bz" ) {
            for( unsigned int i=0; i<nx_d; i++ ) {
                Bz_val[i]=( *field2D )( i, 0 );
            }
        }
    } else if( min_max == 3 && patch->isYmax() ) {
        if( field2D->name=="Bx" ) {
            for( unsigned int i=0; i<nx_p; i++ ) {
                Bx_val[i]=( *field2D )( i, ny_d-1 );
            }
        }
        
        if( field2D->name=="By" ) {
            for( unsigned int i=0; i<nx_d; i++ ) {
                By_val[i]=( *field2D )( i, ny_p-1 );
            }
        }
        
        if( field2D->name=="Bz" ) {
            for( unsigned int i=0; i<nx_d; i++ ) {
                Bz_val[i]=( *field2D )( i, ny_d-1 );
            }
        }
    }
}

void ElectroMagnBC2D_SM::disableExternalFields()
{
    Bx_val.resize( 0 );
    By_val.resize( 0 );
    Bz_val.resize( 0 );
}



// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_SM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    
    if( min_max == 0 && patch->isXmin() ) {
    
        // Static cast of the fields
        //Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
        Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
        Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
        Field2D *By2D = static_cast<Field2D *>( EMfields->By_ );
        Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // for By^(d,p)
        vector<double> yp( 1 );
        for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
        
            double by = 0.;
            yp[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*dy;
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                by += vecLaser[ilaser]->getAmplitude0( yp, time_dual, j, 0 );
            }
            
            ( *By2D )( 0, j ) = Alpha_Xmin   * ( *Ez2D )( 0, j )
                                           +              Beta_Xmin    *( ( *By2D )( 1, j )-By_val[j] )
                                           +              Gamma_Xmin   * by
                                           +              Delta_Xmin   *( ( *Bx2D )( 0, j+1 )-Bx_val[j+1] )
                                           +              Epsilon_Xmin *( ( *Bx2D )( 0, j )-Bx_val[j] )
                                           +              By_val[j];
                                           
        }//j  ---end compute By
        
        
        // for Bz^(d,d)
        vector<double> yd( 1 );
        for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
        
            double bz = 0.;
            yd[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*dy;
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bz += vecLaser[ilaser]->getAmplitude1( yd, time_dual, j, 0 );
            }
            
            /*(*Bz2D)(0,j) = -Alpha_Xmin * (*Ey2D)(0,j)
             +               Beta_Xmin  * (*Bz2D)(1,j)
             +               Gamma_Xmin * bzW;*/
            ( *Bz2D )( 0, j ) = -Alpha_Xmin * ( *Ey2D )( 0, j )
                                           +               Beta_Xmin  *( ( *Bz2D )( 1, j )- Bz_val[j] )
                                           +               Gamma_Xmin * bz
                                           +               Bz_val[j];
                                           
        }//j  ---end compute Bz
        
    } else if( min_max == 1 && patch->isXmax() ) {
    
        // Static cast of the fields
        //Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
        Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
        Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
        Field2D *By2D = static_cast<Field2D *>( EMfields->By_ );
        Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // for By^(d,p)
        vector<double> yp( 1 );
        for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
        
            double by = 0.;
            yp[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*dy;
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                by += vecLaser[ilaser]->getAmplitude0( yp, time_dual, j, 0 );
            }
            
            /*(*By2D)(nx_d-1,j) = Alpha_Xmax   * (*Ez2D)(nx_p-1,j)
             +                   Beta_Xmax    * (*By2D)(nx_d-2,j)
             +                   Gamma_Xmax   * byE
             +                   Delta_Xmax   * (*Bx2D)(nx_p-1,j+1) // Check x-index
             +                   Epsilon_Xmax * (*Bx2D)(nx_p-1,j);*/
            ( *By2D )( nx_d-1, j ) = Alpha_Xmax   * ( *Ez2D )( nx_p-1, j )
                                                +                   Beta_Xmax    *( ( *By2D )( nx_d-2, j ) -By_val[j] )
                                                +                   Gamma_Xmax   * by
                                                +                   Delta_Xmax   *( ( *Bx2D )( nx_p-1, j+1 ) -Bx_val[j+1] ) // Check x-index
                                                +                   Epsilon_Xmax *( ( *Bx2D )( nx_p-1, j ) -Bx_val[j] )
                                                +                   By_val[j];
                                                
        }//j  ---end compute By
        
        
        // for Bz^(d,d)
        vector<double> yd( 1 );
        for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
        
            double bz = 0.;
            yd[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*dy;
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bz += vecLaser[ilaser]->getAmplitude1( yd, time_dual, j, 0 );
            }
            
            /*(*Bz2D)(nx_d-1,j) = -Alpha_Xmax * (*Ey2D)(nx_p-1,j)
             +                    Beta_Xmax  * (*Bz2D)(nx_d-2,j)
             +                    Gamma_Xmax * bzE;*/
            ( *Bz2D )( nx_d-1, j ) = -Alpha_Xmax * ( *Ey2D )( nx_p-1, j )
                                                +                    Beta_Xmax  *( ( *Bz2D )( nx_d-2, j ) -Bz_val[j] )
                                                +                    Gamma_Xmax * bz
                                                +                    Bz_val[j];
                                                
        }//j  ---end compute Bz
    
    } else if( min_max == 2 && patch->isYmin() ) {
    
        // Static cast of the fields
        Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
        //Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
        Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
        Field2D *By2D = static_cast<Field2D *>( EMfields->By_ );
        Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // for Bx^(p,d)
        vector<double> xp( 1 );
        for( unsigned int j=patch->isXmin() ; j<nx_p-patch->isXmax() ; j++ ) {
            
            double bx = 0.;
            xp[0] = patch->getDomainLocalMin( 0 ) + ( ( int )j - ( int )EMfields->oversize[0] )*dx;
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bx += vecLaser[ilaser]->getAmplitude0( xp, time_dual, j, 0 );
            }
            
            ( *Bx2D )( j, 0 ) = -Alpha_Ymin   * ( *Ez2D )( j, 0 )
                                           +               Beta_Ymin    *( ( *Bx2D )( j, 1 )-Bx_val[j] )
                                           +               Gamma_Ymin   * bx
                                           +               Delta_Ymin   *( ( *By2D )( j+1, 0 )-By_val[j+1] )
                                           +               Epsilon_Ymin *( ( *By2D )( j, 0 )-By_val[j] )
                                           +               Bx_val[j];
        }//j  ---end Bx
        
        
        // for Bz^(d,d)
        vector<double> xd( 1 );
        for( unsigned int j=patch->isXmin() ; j<nx_d-patch->isXmax() ; j++ ) {
            
            double bz = 0.;
            xd[0] = patch->getDomainLocalMin( 0 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[0] )*dx;
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bz += vecLaser[ilaser]->getAmplitude1( xd, time_dual, j, 0 );
            }
            
            ( *Bz2D )( j, 0 ) = Alpha_Ymin * ( *Ex2D )( j, 0 )
                                           +               Beta_Ymin  *( ( *Bz2D )( j, 1 )-Bz_val[j] )
                                           +               Gamma_Ymin * bz
                                           +               Bz_val[j];
        }//j  ---end Bz
        
    } else if( min_max == 3 && patch->isYmax() ) {
    
        // Static cast of the fields
        Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
        //Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
        Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
        Field2D *By2D = static_cast<Field2D *>( EMfields->By_ );
        Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // for Bx^(p,d)
        vector<double> xp( 1 );
        for( unsigned int j=patch->isXmin() ; j<nx_p-patch->isXmax() ; j++ ) {
        
            double bx = 0.;
            xp[0] = patch->getDomainLocalMin( 0 ) + ( ( int )j - ( int )EMfields->oversize[0] )*dx;
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bx += vecLaser[ilaser]->getAmplitude0( xp, time_dual, j, 0 );
            }
            
            ( *Bx2D )( j, ny_d-1 ) = -Alpha_Ymax   * ( *Ez2D )( j, ny_p-1 )
                                                +                   Beta_Ymax    *( ( *Bx2D )( j, ny_d-2 ) -Bx_val[j] )
                                                +                   Gamma_Ymax   * bx
                                                +                   Delta_Ymax   *( ( *By2D )( j+1, ny_p-1 ) -By_val[j+1] )
                                                +                   Epsilon_Ymax *( ( *By2D )( j, ny_p-1 ) -By_val[j] )
                                                +                   Bx_val[j];
        }//j  ---end Bx
        
        // for Bz^(d,d)
        vector<double> xd( 1 );
        for( unsigned int j=patch->isXmin() ; j<nx_d-patch->isXmax() ; j++ ) {
        
            double bz = 0.;
            xd[0] = patch->getDomainLocalMin( 0 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[0] )*dx;
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bz += vecLaser[ilaser]->getAmplitude1( xd, time_dual, j, 0 );
            }
            
            ( *Bz2D )( j, ny_d-1 ) = Alpha_Ymax * ( *Ex2D )( j, ny_p-1 )
                                                +                   Beta_Ymax  *( ( *Bz2D )( j, ny_d-2 )- Bz_val[j] )
                                                +                   Gamma_Ymax * bz
                                                +                   Bz_val[j];
        }//j  ---end Bx
        
        
    }
}
