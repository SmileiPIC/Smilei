#include "ElectroMagnBC3D_BM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field3D.h"
#include "Tools.h"
#include "Laser.h"

using namespace std;

ElectroMagnBC3D_BM::ElectroMagnBC3D_BM( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC3D( params, patch, _min_max )
{

    std::vector<unsigned int> dims( 2, 0 );
    
    Bx_val = By_val = Bz_val = nullptr;
    
    if( min_max==0 && patch->isXmin() ) {
        // BCs at the x-border min
        dims = { ny_d, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims );
        Bx_val->put_to( 0. );
        dims = { ny_p, nz_d }; // By^(d,p,d)
        By_val = new Field2D( dims );
        By_val->put_to( 0. );
        dims = { ny_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims );
        Bz_val->put_to( 0. );
    } else if( min_max==1 && patch->isXmax() ) {
        // BCs at the x-border max
        dims = { ny_d, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims );
        Bx_val->put_to( 0. );
        dims = { ny_p, nz_d }; // By^(d,p,d)
        By_val = new Field2D( dims );
        By_val->put_to( 0. );
        dims = { ny_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims );
        Bz_val->put_to( 0. );
    } else if( min_max==2 && patch->isYmin() ) {
        // BCs in the y-border min
        dims = { nx_p, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims );
        Bx_val->put_to( 0. );
        dims = { nx_d, nz_d }; // By^(d,p,d)
        By_val = new Field2D( dims );
        By_val->put_to( 0. );
        dims = { nx_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims );
        Bz_val->put_to( 0. );
    } else if( min_max==3 && patch->isYmax() ) {
        // BCs in the y-border mix
        dims = { nx_p, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims );
        Bx_val->put_to( 0. );
        dims = { nx_d, nz_d }; // By^(d,p,d)
        By_val = new Field2D( dims );
        By_val->put_to( 0. );
        dims = { nx_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims );
        Bz_val->put_to( 0. );
    } else if( min_max==4 && patch->isZmin() ) {
        // BCs in the z-border min
        dims = { nx_p, ny_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims );
        Bx_val->put_to( 0. );
        dims = { nx_d, ny_p }; // By^(d,p,d)
        By_val = new Field2D( dims );
        By_val->put_to( 0. );
        dims = { nx_d, ny_d }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims );
        Bz_val->put_to( 0. );
    } else if( min_max==5 && patch->isZmax() ) {
        // BCs in the z-border max
        dims = { nx_p, ny_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims );
        Bx_val->put_to( 0. );
        dims = { nx_d, ny_p }; // By^(d,p,d)
        By_val = new Field2D( dims );
        By_val->put_to( 0. );
        dims = { nx_d, ny_d }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims );
        Bz_val->put_to( 0. );
    }
    
    
    // -----------------------------------------------------
    // Parameters for the Buneman boundary conditions
    // -----------------------------------------------------
    
    //! \todo (AB) Check optimal angle for buneman BC
    
    for( unsigned int idim=0; idim <3; idim++ ) {
        for( unsigned int minmax=0; minmax <2; minmax++ ) {
            double kn  = params.EM_BCs_k[idim*2 + minmax][idim];
            cb[idim][minmax] =  kn / ( 1.0 + kn );
            ce[idim][minmax] =  1. - cb[idim][minmax] ;
        }
    }
    
    // X boundary
    Alpha_BM_x    = ( dt_ov_dx - 1. )  / ( dt_ov_dx + 1. ) ;
    Beta_BM_x     =  dt_ov_dy        / ( dt_ov_dx + 1. ) ;
    Gamma_BM_x    =  dt_ov_dz        / ( dt_ov_dx + 1. ) ;
    
    // Y boundary
    Alpha_BM_y    = ( dt_ov_dy - 1. )  / ( dt_ov_dy + 1. ) ;
    Beta_BM_y     =  dt_ov_dx        / ( dt_ov_dy + 1. ) ;
    Gamma_BM_y    =  dt_ov_dz        / ( dt_ov_dy + 1. ) ;
    
    
    // Zmin boundary
    Alpha_BM_z    = ( dt_ov_dz - 1. )  / ( dt_ov_dz + 1. ) ;
    Beta_BM_z     =  dt_ov_dy        / ( dt_ov_dz + 1. ) ;
    Gamma_BM_z    =  dt_ov_dx        / ( dt_ov_dz + 1. ) ;
    
    
}

ElectroMagnBC3D_BM::~ElectroMagnBC3D_BM()
{
    if( Bx_val ) {
        delete Bx_val ;
    }
    if( By_val ) {
        delete By_val ;
    }
    if( Bz_val ) {
        delete Bz_val ;
    }
}

void ElectroMagnBC3D_BM::save_fields( Field *my_field, Patch *patch )
{
    Field3D *field3D=static_cast<Field3D *>( my_field );
    
    if( min_max==0 && patch->isXmin() ) {
    
        if( field3D->name=="Bx" ) {
            field3D->extract_slice_yz( 0,      Bx_val );
        } else if( field3D->name=="By" ) {
            field3D->extract_slice_yz( 0,      By_val );
        } else if( field3D->name=="Bz" ) {
            field3D->extract_slice_yz( 0,      Bz_val );
        }
    } else if( min_max==1 && patch->isXmax() ) {
        if( field3D->name=="Bx" ) {
            field3D->extract_slice_yz( 0,      Bx_val );
        } else if( field3D->name=="By" ) {
            field3D->extract_slice_yz( 0,      By_val );
        } else if( field3D->name=="Bz" ) {
            field3D->extract_slice_yz( 0,      Bz_val );
        }
    } else if( min_max==2 && patch->isYmin() ) {
        if( field3D->name=="Bx" ) {
            field3D->extract_slice_xz( 0,      Bx_val );
        } else if( field3D->name=="By" ) {
            field3D->extract_slice_xz( 0,      By_val );
        } else if( field3D->name=="Bz" ) {
            field3D->extract_slice_xz( 0,      Bz_val );
        }
    } else if( min_max==3 && patch->isYmax() ) {
        if( field3D->name=="Bx" ) {
            field3D->extract_slice_xz( ny_d-1, Bx_val );
        } else if( field3D->name=="By" ) {
            field3D->extract_slice_xz( ny_p-1, By_val );
        } else if( field3D->name=="Bz" ) {
            field3D->extract_slice_xz( ny_d-1, Bz_val );
        }
    } else if( min_max==4 && patch->isZmin() ) {
    
        if( field3D->name=="Bx" ) {
            field3D->extract_slice_xy( 0,      Bx_val );
        } else if( field3D->name=="By" ) {
            field3D->extract_slice_xy( 0,      By_val );
        } else if( field3D->name=="Bz" ) {
            field3D->extract_slice_xy( 0,      Bz_val );
        }
    } else if( min_max==5 && patch->isZmax() ) {
    
        if( field3D->name=="Bx" ) {
            field3D->extract_slice_xy( nz_d-1, Bx_val );
        } else if( field3D->name=="By" ) {
            field3D->extract_slice_xy( nz_d-1, By_val );
        } else if( field3D->name=="Bz" ) {
            field3D->extract_slice_xy( nz_p-1, Bz_val );
        }
    }
}


void ElectroMagnBC3D_BM::disableExternalFields()
{
    delete Bx_val;
    Bx_val = NULL;
    delete By_val;
    By_val = NULL;
    delete Bz_val;
    Bz_val = NULL;
}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_BM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    if( min_max==0 && patch->isXmin() ) {
    
        // Static cast of the fields
        Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        //Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *Bx3D_old = static_cast<Field3D *>( EMfields->Bx_m );
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        Field3D *By3D_old = static_cast<Field3D *>( EMfields->By_m );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        Field3D *Bz3D_old = static_cast<Field3D *>( EMfields->Bz_m );
        
        unsigned const int i = 0 ;
        
        
        // for By^(d,p,d)
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) { //Undefined for extreme k
                ( *By3D )( i, j, k ) = ( *By3D_old )( i+1, j, k )
                                       +                 Alpha_BM_x    * ( ( *By3D )( i+1, j, k ) - ( *By3D_old )( i, j, k ) )
                                       -                 cb[0][0]*Beta_BM_x  * ( ( *Bx3D )( i, j+1, k ) - ( *Bx3D )( i, j, k ) + ( *Bx3D_old )( i, j+1, k ) - ( *Bx3D_old )( i, j, k ) )
                                       -                 ce[0][0]*Gamma_BM_x * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j, k-1 ) + ( *Ex3D )( i+1, j, k ) - ( *Ex3D )( i+1, j, k-1 ) ) ;
            }// k  ---end compute By
        }//j  ---end compute By
        
        // for Bz^(d,d,p)
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) { //Undefined for extreme j
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( i, j, k ) = ( *Bz3D_old )( i+1, j, k )
                                       +                 Alpha_BM_x    * ( ( *Bz3D )( i+1, j, k ) - ( *Bz3D_old )( i, j, k ) )
                                       -                 cb[0][0]*Gamma_BM_x * ( ( *Bx3D )( i, j, k+1 ) - ( *Bx3D )( i, j, k ) + ( *Bx3D_old )( i, j, k+1 ) - ( *Bx3D_old )( i, j, k ) )
                                       +                 ce[0][0]*Beta_BM_x  * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j-1, k ) + ( *Ex3D )( i+1, j, k ) - ( *Ex3D )( i+1, j-1, k ) ) ;
                                       
            }// k  ---end compute Bz
        }//j  ---end compute Bz
    } else if( min_max==1 && patch->isXmax() ) {
    
        // Static cast of the fields
        Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        //Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *Bx3D_old = static_cast<Field3D *>( EMfields->Bx_m );
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        Field3D *By3D_old = static_cast<Field3D *>( EMfields->By_m );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        Field3D *Bz3D_old = static_cast<Field3D *>( EMfields->Bz_m );
        
        unsigned const int i = nx_d -2 ;
        
        // for By^(d,p,d)
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) { //Undefined for extreme k
                ( *By3D )( i+1, j, k ) = ( *By3D_old )( i, j, k )
                                         +                 Alpha_BM_x    * ( ( *By3D )( i, j, k ) - ( *By3D_old )( i+1, j, k ) )
                                         +                 cb[0][1]*Beta_BM_x  * ( ( *Bx3D )( i, j+1, k ) - ( *Bx3D )( i, j, k ) + ( *Bx3D_old )( i, j+1, k ) - ( *Bx3D_old )( i, j, k ) )
                                         -                 ce[0][1]*Gamma_BM_x * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j, k-1 ) + ( *Ex3D )( i+1, j, k ) - ( *Ex3D )( i+1, j, k-1 ) ) ;
                                         
            }//k  ---end compute By
        }//j  ---end compute By
        
        // for Bz^(d,d,p)
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) { //Undefined for extreme j
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( i+1, j, k ) = ( *Bz3D_old )( i, j, k )
                                         +                 Alpha_BM_x    * ( ( *Bz3D )( i, j, k ) - ( *Bz3D_old )( i+1, j, k ) )
                                         +                 cb[0][1]*Gamma_BM_x * ( ( *Bx3D )( i, j, k+1 ) - ( *Bx3D )( i, j, k ) + ( *Bx3D_old )( i, j, k+1 ) - ( *Bx3D_old )( i, j, k ) )
                                         +                 ce[0][1]*Beta_BM_x  * ( ( *Ex3D )( i, j, k )   - ( *Ex3D )( i, j-1, k ) + ( *Ex3D )( i+1, j, k ) - ( *Ex3D )( i+1, j-1, k ) ) ;
            }//k  ---end compute Bz
        }//j  ---end compute Bz
    } else if( min_max==2 && patch->isYmin() ) {
    
        // Static cast of the fields
        //Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
        //Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *Bx3D_old = static_cast<Field3D *>( EMfields->Bx_m );
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        Field3D *By3D_old = static_cast<Field3D *>( EMfields->By_m );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        Field3D *Bz3D_old = static_cast<Field3D *>( EMfields->Bz_m );
        
        unsigned const int j = 0;
        
        // for Bx^(p,d,d)
        for( unsigned int i=0 ; i<nx_p ; i++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) { // undefined for k=0 and k = nz_d-1 !!
                ( *Bx3D )( i, j, k ) = ( *Bx3D_old )( i, j+1, k )
                                       +                 Alpha_BM_y    * ( ( *Bx3D )( i, j+1, k ) - ( *Bx3D_old )( i, j, k ) )
                                       -                 cb[1][0]*Beta_BM_y  * ( ( *By3D )( i+1, j, k ) - ( *By3D )( i, j, k )  + ( *By3D_old )( i+1, j, k ) - ( *By3D_old )( i, j, k ) )
                                       +                 ce[1][0]*Gamma_BM_y * ( ( *Ey3D )( i, j, k ) - ( *Ey3D )( i, j, k-1 )+ ( *Ey3D )( i, j+1, k ) - ( *Ey3D )( i, j+1, k-1 ) ) ;
            }// k  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bz^(d,d,p)
        for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {    // undefined for i=0 and i = nx_d-1 !!
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( i, j, k ) = ( *Bz3D_old )( i, j+1, k )
                                       +                 Alpha_BM_y    * ( ( *Bz3D )( i, j+1, k ) - ( *Bz3D_old )( i, j, k ) )
                                       -                 cb[1][0]*Gamma_BM_y * ( ( *By3D )( i, j, k+1 ) - ( *By3D )( i, j, k ) + ( *By3D_old )( i, j, k+1 ) - ( *By3D_old )( i, j, k ) )
                                       -                 ce[1][0]*Beta_BM_y  * ( ( *Ey3D )( i, j, k ) - ( *Ey3D )( i-1, j, k ) + ( *Ey3D )( i, j+1, k )   - ( *Ey3D )( i-1, j+1, k ) ) ;
            }// k  ---end compute Bz
        }//i  ---end compute Bz
        
    } else if( min_max==3 && patch->isYmax() ) {
    
        // Static cast of the fields
        //Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
        //Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *Bx3D_old = static_cast<Field3D *>( EMfields->Bx_m );
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        Field3D *By3D_old = static_cast<Field3D *>( EMfields->By_m );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        Field3D *Bz3D_old = static_cast<Field3D *>( EMfields->Bz_m );
        
        unsigned const int j = ny_d-2 ;
        
        // for Bx^(p,d,d)
        for( unsigned int i=0 ; i<nx_p ; i++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) { // undefined for k=0 and k = nz_d-1 !!
                ( *Bx3D )( i, j+1, k ) = ( *Bx3D_old )( i, j, k )
                                         +                 Alpha_BM_y    * ( ( *Bx3D )( i, j, k ) - ( *Bx3D_old )( i, j+1, k ) )
                                         +                 cb[1][1]*Beta_BM_y  * ( ( *By3D )( i+1, j, k ) - ( *By3D )( i, j, k ) + ( *By3D_old )( i+1, j, k ) - ( *By3D_old )( i, j, k ) )
                                         +                 ce[1][1]*Gamma_BM_y * ( ( *Ey3D )( i, j, k ) - ( *Ey3D )( i, j, k-1 ) + ( *Ey3D )( i, j+1, k ) - ( *Ey3D )( i, j+1, k-1 ) ) ;
            }// k  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bz^(d,d,p)
        for( unsigned int i=1 ; i<nx_d-1 ; i++ ) { // undefined for i=0 and i = nx_d-1 !!
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( i, j+1, k ) = ( *Bz3D_old )( i, j, k )
                                         +                 Alpha_BM_y    * ( ( *Bz3D )( i, j, k ) - ( *Bz3D_old )( i, j+1, k ) )
                                         +                 cb[1][1]*Gamma_BM_y * ( ( *By3D )( i, j, k+1 ) - ( *By3D )( i, j, k ) + ( *By3D_old )( i, j, k+1 ) - ( *By3D_old )( i, j, k ) )
                                         -                 ce[1][1]*Beta_BM_y  * ( ( *Ey3D )( i, j, k )   - ( *Ey3D )( i-1, j, k ) + ( *Ey3D )( i, j+1, k )   - ( *Ey3D )( i-1, j+1, k ) ) ;
            }// k  ---end compute Bz
        }//i  ---end compute Bz
        
    } else if( min_max==4 && patch->isZmin() ) {
    
        // Static cast of the fields
        //Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *Bx3D_old = static_cast<Field3D *>( EMfields->Bx_m );
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        Field3D *By3D_old = static_cast<Field3D *>( EMfields->By_m );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        Field3D *Bz3D_old = static_cast<Field3D *>( EMfields->Bz_m );
        
        unsigned const int k = 0;
        
        // for By^(d,p,d)
        for( unsigned int i=1 ; i<nx_d-1 ; i++ ) { // undefined for i=0 and i = nx_d-1 !!
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
            
                ( *By3D )( i, j, k ) = ( *By3D_old )( i, j, k+1 )
                                       +                 Alpha_BM_z    * ( ( *By3D )( i, j, k+1 ) - ( *By3D_old )( i, j, k ) )
                                       -                 cb[2][0]*Beta_BM_z  * ( ( *Bz3D )( i, j+1, k )   - ( *Bz3D )( i, j, k ) + ( *Bz3D_old )( i, j+1, k )   - ( *Bz3D_old )( i, j, k ) )
                                       +                 ce[2][0]*Gamma_BM_z * ( ( *Ez3D )( i, j, k )   - ( *Ez3D )( i-1, j, k ) + ( *Ez3D )( i, j, k+1 ) - ( *Ez3D )( i-1, j, k+1 ) ) ;
            }// j  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bx^(p,d,d)
        for( unsigned int i=0 ; i<nx_p ; i++ ) {
            for( unsigned int j=1 ; j<ny_d-1 ; j++ ) { // undefined for j=0 and j = ny_d-1 !!
            
                ( *Bx3D )( i, j, k ) = ( *Bx3D_old )( i, j, k+1 )
                                       +                 Alpha_BM_z    * ( ( *Bx3D )( i, j, k+1 ) - ( *Bx3D_old )( i, j, k ) )
                                       -                 cb[2][0]*Gamma_BM_z * ( ( *Bz3D )( i+1, j, k ) - ( *Bz3D )( i, j, k ) + ( *Bz3D_old )( i+1, j, k )   - ( *Bz3D_old )( i, j, k ) )
                                       -                 ce[2][0]*Beta_BM_z  * ( ( *Ez3D )( i, j, k ) - ( *Ez3D )( i, j-1, k ) + ( *Ez3D )( i, j, k+1 ) - ( *Ez3D )( i, j-1, k+1 ) ) ;
                                       
            }// j  ---end compute By
        }//i  ---end compute By
        
    } else if( min_max==5 && patch->isZmax() ) {
    
        // Static cast of the fields
        //Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *Bx3D_old = static_cast<Field3D *>( EMfields->Bx_m );
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        Field3D *By3D_old = static_cast<Field3D *>( EMfields->By_m );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        Field3D *Bz3D_old = static_cast<Field3D *>( EMfields->Bz_m );
        
        unsigned const int k = nz_d-2 ;
        
        // for By^(d,p,d)
        for( unsigned int i=1 ; i<nx_d-1 ; i++ ) { // undefined for i=0 and i = nx_d-1 !!
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
            
                ( *By3D )( i, j, k+1 ) = ( *By3D_old )( i, j, k )
                                         +                 Alpha_BM_z    * ( ( *By3D )( i, j, k ) - ( *By3D_old )( i, j, k+1 ) )
                                         +                 cb[2][1]*Beta_BM_z  * ( ( *Bz3D )( i, j+1, k )   - ( *Bz3D )( i, j, k ) + ( *Bz3D_old )( i, j+1, k )   - ( *Bz3D_old )( i, j, k ) )
                                         +                 ce[2][1]*Gamma_BM_z * ( ( *Ez3D )( i, j, k )   - ( *Ez3D )( i-1, j, k ) + ( *Ez3D )( i, j, k+1 ) - ( *Ez3D )( i-1, j, k+1 ) ) ;
            }// j  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bx^(p,d,d)
        for( unsigned int i=0 ; i<nx_p ; i++ ) {
            for( unsigned int j=1 ; j<ny_d-1 ; j++ ) { // undefined for j=0 and j = ny_d-1 !!
            
                ( *Bx3D )( i, j, k+1 ) = ( *Bx3D_old )( i, j, k )
                                         +                 Alpha_BM_z    * ( ( *Bx3D )( i, j, k ) - ( *Bx3D_old )( i, j, k+1 ) )
                                         +                 cb[2][1]*Gamma_BM_z * ( ( *Bz3D )( i+1, j, k ) - ( *Bz3D )( i, j, k ) + ( *Bz3D_old )( i+1, j, k )   - ( *Bz3D_old )( i, j, k ) )
                                         -                 ce[2][1]*Beta_BM_z  * ( ( *Ez3D )( i, j, k ) - ( *Ez3D )( i, j-1, k ) + ( *Ez3D )( i, j, k+1 ) - ( *Ez3D )( i, j-1, k+1 ) ) ;
                                         
            }// j  ---end compute By
        }//i  ---end compute By
        
    }
}
