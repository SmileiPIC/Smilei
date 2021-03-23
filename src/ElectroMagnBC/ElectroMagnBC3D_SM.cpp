#include "ElectroMagnBC3D_SM.h"

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

ElectroMagnBC3D_SM::ElectroMagnBC3D_SM( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC3D( params, patch, _min_max )
{

    std::vector<unsigned int> dims( 2, 0 );
    
    Bx_val = By_val = Bz_val = nullptr;
    
    if( min_max==0 && patch->isXmin() ) {
        // BCs at the x-border min
        dims = { ny_d, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims, "Bx_val" );
        Bx_val->put_to( 0. );
        dims = { ny_p, nz_d }; // By^(d,p,d)
        By_val = new Field2D( dims, "By_val" );
        By_val->put_to( 0. );
        dims = { ny_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims, "Bz_val" );
        Bz_val->put_to( 0. );
    } else if( min_max==1 && patch->isXmax() ) {
        // BCs at the x-border max
        dims = { ny_d, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims, "Bx_val" );
        Bx_val->put_to( 0. );
        dims = { ny_p, nz_d }; // By^(d,p,d)
        By_val = new Field2D( dims, "By_val" );
        By_val->put_to( 0. );
        dims = { ny_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims, "Bz_val" );
        Bz_val->put_to( 0. );
    } else if( min_max==2 && patch->isYmin() ) {
        // BCs in the y-border min
        dims = { nx_p, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims, "Bx_val" );
        Bx_val->put_to( 0. );
        dims = { nx_d, nz_d }; // By^(d,p,d)
        By_val = new Field2D( dims, "By_val" );
        By_val->put_to( 0. );
        dims = { nx_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims, "Bz_val" );
        Bz_val->put_to( 0. );
    } else if( min_max==3 && patch->isYmax() ) {
        // BCs in the y-border mix
        dims = { nx_p, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims, "Bx_val" );
        Bx_val->put_to( 0. );
        dims = { nx_d, nz_d }; // By^(d,p,d)
        By_val = new Field2D( dims, "By_val" );
        By_val->put_to( 0. );
        dims = { nx_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims, "Bz_val" );
        Bz_val->put_to( 0. );
    } else if( min_max==4 && patch->isZmin() ) {
        // BCs in the z-border min
        dims = { nx_p, ny_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims, "Bx_val" );
        Bx_val->put_to( 0. );
        dims = { nx_d, ny_p }; // By^(d,p,d)
        By_val = new Field2D( dims, "By_val" );
        By_val->put_to( 0. );
        dims = { nx_d, ny_d }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims, "Bz_val" );
        Bz_val->put_to( 0. );
    } else if( min_max==5 && patch->isZmax() ) {
        // BCs in the z-border max
        dims = { nx_p, ny_d }; // Bx^(p,d,d)
        Bx_val = new Field2D( dims, "Bx_val" );
        Bx_val->put_to( 0. );
        dims = { nx_d, ny_p }; // By^(d,p,d)
        By_val = new Field2D( dims, "By_val" );
        By_val->put_to( 0. );
        dims = { nx_d, ny_d }; // Bz^(d,d,p)
        Bz_val = new Field2D( dims, "Bz_val" );
        Bz_val->put_to( 0. );
    }
    
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    //! \todo (MG) Check optimal angle for Silver-Muller BCs
    double pyKx, pyKy, pyKz;
    double kx, ky, kz;
    double Knorm;
    double omega = 1. ;
    //kx = w cos(theta) cos(phi)
    //ky = w sin(theta)
    //kz = w cos(theta) sin(phi)
    
    // Xmin boundary
    pyKx = params.EM_BCs_k[0][0];
    pyKy = params.EM_BCs_k[0][1];
    pyKz = params.EM_BCs_k[0][2];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    kz = omega*pyKz/Knorm;
    
    double factor = 1.0 / ( kx + dt_ov_dx );
    Alpha_Xmin    = 2.0 * factor;
    Beta_Xmin     = - ( kx-dt_ov_dx ) * factor;
    Gamma_Xmin    = 4.0 * kx * factor;
    Delta_Xmin    = - ( ky + dt_ov_dy ) * factor;
    Epsilon_Xmin  = - ( ky - dt_ov_dy ) * factor;
    Zeta_Xmin     = - ( kz + dt_ov_dz ) * factor;
    Eta_Xmin      = - ( kz - dt_ov_dz ) * factor;
    
    // Xmax boundary
    pyKx = params.EM_BCs_k[1][0];
    pyKy = params.EM_BCs_k[1][1];
    pyKz = params.EM_BCs_k[1][2];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    kz = omega*pyKz/Knorm;
    
    factor        = 1.0 / ( kx - dt_ov_dx );
    Alpha_Xmax    = 2.0 * factor;
    Beta_Xmax     = - ( kx+dt_ov_dx )  * factor;
    Gamma_Xmax    = 4.0 * kx * factor;
    Delta_Xmax    = - ( ky + dt_ov_dy )  * factor;
    Epsilon_Xmax  = - ( ky - dt_ov_dy )  * factor;
    Zeta_Xmax     = - ( kz + dt_ov_dz ) * factor;
    Eta_Xmax      = - ( kz - dt_ov_dz ) * factor;
    
    // Ymin boundary
    pyKx = params.EM_BCs_k[2][0];
    pyKy = params.EM_BCs_k[2][1];
    pyKz = params.EM_BCs_k[2][2];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    kz = omega*pyKz/Knorm;
    
    factor = 1.0 / ( ky + dt_ov_dy );
    Alpha_Ymin    = 2.0 * factor;
    Beta_Ymin     = - ( ky - dt_ov_dy ) * factor;
    Gamma_Ymin    = 4.0 * ky * factor;
    Delta_Ymin    = - ( kz + dt_ov_dz ) * factor;
    Epsilon_Ymin  = - ( kz -dt_ov_dz ) * factor;
    Zeta_Ymin     = - ( kx + dt_ov_dx ) * factor;
    Eta_Ymin      = - ( kx - dt_ov_dx ) * factor;
    
    // Ymax boundary
    pyKx = params.EM_BCs_k[3][0];
    pyKy = params.EM_BCs_k[3][1];
    pyKz = params.EM_BCs_k[3][2];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    kz = omega*pyKz/Knorm;
    
    factor = 1.0 / ( ky - dt_ov_dy );
    Alpha_Ymax    = 2.0                     * factor;
    Beta_Ymax     = - ( ky + dt_ov_dy ) * factor;
    Gamma_Ymax    = 4.0 * ky * factor;
    Delta_Ymax    = - ( kz + dt_ov_dz ) * factor;
    Epsilon_Ymax  = - ( kz - dt_ov_dz ) * factor;
    Zeta_Ymax     = - ( kx + dt_ov_dx ) * factor;
    Eta_Ymax      = - ( kx - dt_ov_dx ) * factor;
    
    // Zmin boundary
    pyKx = params.EM_BCs_k[4][0];
    pyKy = params.EM_BCs_k[4][1];
    pyKz = params.EM_BCs_k[4][2];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    kz = omega*pyKz/Knorm;
    
    factor = 1.0 / ( kz + dt_ov_dz );
    Alpha_Zmin    = 2.0                     * factor;
    Beta_Zmin     = - ( kz - dt_ov_dz ) * factor;
    Delta_Zmin    = - ( kx + dt_ov_dx ) * factor;
    Epsilon_Zmin  = - ( kx - dt_ov_dx ) * factor;
    Zeta_Zmin     = - ( ky + dt_ov_dy ) * factor;
    Eta_Zmin      = - ( ky - dt_ov_dy ) * factor;
    
    // Zmax boundary
    pyKx = params.EM_BCs_k[5][0];
    pyKy = params.EM_BCs_k[5][1];
    pyKz = params.EM_BCs_k[5][2];
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
    kx = omega*pyKx/Knorm;
    ky = omega*pyKy/Knorm;
    kz = omega*pyKz/Knorm;
    
    factor        = 1.0 / ( kz - dt_ov_dz );
    Alpha_Zmax    = 2.0                      * factor;
    Beta_Zmax     = - ( kz + dt_ov_dz )  * factor;
    Delta_Zmax    = - ( kx + dt_ov_dx )  * factor;
    Epsilon_Zmax  = - ( kx - dt_ov_dx )  * factor;
    Zeta_Zmax     = - ( ky + dt_ov_dy ) * factor;
    Eta_Zmax      = - ( ky - dt_ov_dy ) * factor;
    
}

ElectroMagnBC3D_SM::~ElectroMagnBC3D_SM()
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


// Magnetic field Bx^(p,d,d)
// Magnetic field By^(d,p,d)
// Magnetic field Bz^(d,d,p)

void ElectroMagnBC3D_SM::save_fields( Field *my_field, Patch *patch )
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


void ElectroMagnBC3D_SM::disableExternalFields()
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
void ElectroMagnBC3D_SM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{

    // Static cast of the fields
    double *Ex3D = &(EMfields->Ex_->data_[0]);
    double *Ey3D = &(EMfields->Ey_->data_[0]);
    double *Ez3D = &(EMfields->Ez_->data_[0]);
    double *Bx3D = &(EMfields->Bx_->data_[0]);
    double *By3D = &(EMfields->By_->data_[0]);
    double *Bz3D = &(EMfields->Bz_->data_[0]);
    vector<double> pos( 2 );
    
    double* Bx_ext = NULL;
    if( Bx_val!=nullptr ) {
        Bx_ext = &(Bx_val->data_[0]);
    }
    double* By_ext = NULL;
    if( By_val!=nullptr ) {
        By_ext = &(By_val->data_[0]);
    }
    double* Bz_ext = NULL;
    if( Bz_val!=nullptr ) {
        Bz_ext = &(Bz_val->data_[0]);
    }
    
    if( min_max==0 && patch->isXmin() ) {
    
        // for By^(d,p,d)
        vector<double> by( ny_p*nz_d, 0. );
        for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*dy;
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k -0.5 - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    by[ j*nz_d+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, k );
                }
            }
        }
        for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                
                By3D[ 0*(ny_p*nz_d) + j*nz_d + k ] = Alpha_Xmin   * Ez3D[ 0*(ny_p*nz_d) + j*nz_d + k ]
                                       +              Beta_Xmin    *( By3D[ 1*(ny_p*nz_d) + j*nz_d + k ]-By_ext[ j*nz_d + k ] )
                                       +              Gamma_Xmin   * by[ j*nz_d+k ]
                                       +              Delta_Xmin   *( Bx3D[ 0*(ny_d*nz_d) + (j+1)*nz_d + k ]-Bx_ext[ (j+1)*nz_d + k ] )
                                       +              Epsilon_Xmin *( Bx3D[ 0*(ny_d*nz_d) +  j   *nz_d + k ]-Bx_ext[  j   *nz_d + k ] )
                                       + By_ext[ j*nz_d + k ];
            }// k  ---end compute By
        }//j  ---end compute By
        
        // for Bz^(d,d,p)
        vector<double> bz( ny_d*nz_p, 0. );
        for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*dy;
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bz[ j*nz_p+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, k );
                }
            }
        }
        for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                Bz3D[ 0*(ny_d*nz_p) + j*nz_p + k ] = - Alpha_Xmin   * Ey3D[ 0*(ny_d*nz_p) + j*nz_p + k ]
                                       +              Beta_Xmin    *( Bz3D[ 1*(ny_d*nz_p) + j*nz_p + k ]-Bz_ext[ j*nz_p + k ] )
                                       +              Gamma_Xmin   * bz[ j*nz_p+k ]
                                       +              Zeta_Xmin    *( Bx3D[ 0*(ny_d*nz_d) + j*nz_d + k+1 ]-Bx_ext[ j*nz_d + (k+1) ] )
                                       +              Eta_Xmin     *( Bx3D[ 0*(ny_d*nz_d) + j*nz_d + k   ]-Bx_ext[ j*nz_d +  k    ] )
                                       + Bz_ext[ j*nz_p + k ];
                                       
            }// k  ---end compute Bz
        }//j  ---end compute Bz
    } else if( min_max==1 && patch->isXmax() ) {
    
        // for By^(d,p,d)
        vector<double> by( ny_p*nz_d, 0. );
        for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*dy;
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - 0.5 - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    by[ j*nz_d+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, k );
                }
            }
        }
        for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                By3D[ (nx_d-1)*(ny_p*nz_d) + j*nz_d + k ] = Alpha_Xmax   * Ez3D[ (nx_p-1)*(ny_p*nz_d) + j*nz_d + k ]
                                            +                   Beta_Xmax    *( By3D[ (nx_d-2)*(ny_p*nz_d) + j*nz_d + k ] -By_ext[ j*nz_d + k ] )
                                            +                   Gamma_Xmax   * by[ j*nz_d+k ]
                                            +                   Delta_Xmax   *( Bx3D[ (nx_p-1)*(ny_d*nz_d) + (j+1)*nz_d + k ] -Bx_ext[ (j+1)*nz_d + k ] ) // Check x-index
                                            +                   Epsilon_Xmax *( Bx3D[ (nx_p-1)*(ny_d*nz_d) +  j   *nz_d + k ] -Bx_ext[  j   *nz_d + k ] )
                                            + By_ext[ j*nz_d + k ];
                                            
            }//k  ---end compute By
        }//j  ---end compute By
        
        // for Bz^(d,d,p)
        vector<double> bz( ny_d*nz_p, 0. );
        for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax(); j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*dy;
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bz[ j*nz_p+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, k );
                }
            }
        }
        for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax(); j++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                Bz3D[ (nx_d-1)*(ny_d*nz_p) + j*nz_p + k ] = -Alpha_Xmax * Ey3D[ (nx_p-1)*(ny_d*nz_p) + j*nz_p + k ]
                                            +                    Beta_Xmax  *( Bz3D[ (nx_d-2)*(ny_d*nz_p) + j*nz_p + k ] -Bz_ext[ j*nz_p + k ] )
                                            +                    Gamma_Xmax * bz[ j*nz_p+k ]
                                            +                    Zeta_Xmax  *( Bx3D[ (nx_p-1)*(ny_d*nz_d) + j*nz_d + k+1 ]-Bx_ext[ j*nz_d + (k+1) ] )
                                            +                    Eta_Xmax   *( Bx3D[ (nx_p-1)*(ny_d*nz_d) + j*nz_d + k   ]-Bx_ext[ j*nz_d +  k    ] )
                                            + Bz_ext[ j*nz_p + k ];
            }//k  ---end compute Bz
        }//j  ---end compute Bz
    
    } else if( min_max==2 && patch->isYmin() ) {
    
        // for Bx^(p,d,d)
        vector<double> bx( nx_p*nz_d, 0. );
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - ( int )EMfields->oversize[0] )*dx;
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k -0.5 - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bx[ i*nz_d+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, i, k );
                }
            }
        }
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                Bx3D[ i*(ny_d*nz_d) + 0*nz_d + k ] = - Alpha_Ymin   * Ez3D[ i*(ny_p*nz_d) + 0*nz_d + k ]
                                       +              Beta_Ymin     *( Bx3D[  i   *(ny_d*nz_d) + 1*nz_d + k ]-Bx_ext[  i   *nz_d + k ] )
                                       +              Gamma_Ymin   * bx[ i*nz_d+k ]
                                       +              Zeta_Ymin     *( By3D[ (i+1)*(ny_p*nz_d) + 0*nz_d + k ]-By_ext[ (i+1)*nz_d + k ] )
                                       +              Eta_Ymin      *( By3D[  i   *(ny_p*nz_d) + 0*nz_d + k ]-By_ext[  i   *nz_d + k ] )
                                       + Bx_ext[ i*nz_d + k ];
            }// k  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bz^(d,d,p)
        vector<double> bz( nx_d*nz_p, 0. );
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i -0.5 - ( int )EMfields->oversize[0] )*dx;
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bz[ i*nz_p+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, i, k );
                }
            }
        }
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                Bz3D[ i*(ny_d*nz_p) + 0*nz_p + k ] = Alpha_Ymin   * Ex3D[ i*(ny_p*nz_p) + 0*nz_p + k ]
                                       +              Beta_Ymin    *( Bz3D[ i*(ny_d*nz_p) + 1*nz_p + k   ]-Bz_ext[ i*nz_p +  k    ] )
                                       +              Gamma_Ymin   * bz[ i*nz_p+k ]
                                       +              Delta_Ymin   *( By3D[ i*(ny_p*nz_d) + 0*nz_d + k+1 ]-By_ext[ i*nz_d + (k+1) ] )
                                       +              Epsilon_Ymin *( By3D[ i*(ny_p*nz_d) + 0*nz_d + k   ]-By_ext[ i*nz_d +  k    ] )
                                       + Bz_ext[ i*nz_p + k ];
            }// k  ---end compute Bz
        }//i  ---end compute Bz
    
    } else if( min_max==3 && patch->isYmax() ) {
    
        // for Bx^(p,d,d)
        vector<double> bx( nx_p*nz_d, 0. );
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - ( int )EMfields->oversize[0] )*dx;
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k -0.5 - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bx[ i*nz_d+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, i, k );
                }
            }
        }
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
            
                Bx3D[ i*(ny_d*nz_d) + (ny_d-1)*nz_d + k ] = -Alpha_Ymax * Ez3D[ i*(ny_p*nz_d) + (ny_p-1)*nz_d + k ]
                                            +                    Beta_Ymax  *( Bx3D[  i   *(ny_d*nz_d) + (ny_d-2)*nz_d + k ]-Bx_ext[  i   *nz_d + k ] )
                                            +                    Gamma_Ymax   * bx[ i*nz_d+k ]
                                            +                    Zeta_Ymax  *( By3D[ (i+1)*(ny_p*nz_d) + (ny_p-1)*nz_d + k ]-By_ext[ (i+1)*nz_d + k ] )
                                            +                    Eta_Ymax   *( By3D[  i   *(ny_p*nz_d) + (ny_p-1)*nz_d + k ]-By_ext[  i   *nz_d + k ] )
                                            + Bx_ext[ i*nz_d + k ];
                                            
            }//k  ---end compute Bz
        }//j  ---end compute Bz
        
        // for Bz^(d,d,p)
        vector<double> bz( nx_d*nz_p, 0. );
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i -0.5 - ( int )EMfields->oversize[0] )*dx;
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bz[ i*nz_p+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, i, k );
                }
            }
        }
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
            
                Bz3D[ i*(ny_d*nz_p) + (ny_d-1)*nz_p + k ] = Alpha_Ymax   * Ex3D[ i*(ny_p*nz_p) + (ny_p-1)*nz_p + k ]
                                            +                   Beta_Ymax    *( Bz3D[ i*(ny_d*nz_p) + (ny_d-2)*nz_p + k   ] -Bz_ext[ i*nz_p +  k    ] )
                                            +                   Gamma_Ymax * bz[ i*nz_p+k ]
                                            +                   Delta_Ymax   *( By3D[ i*(ny_p*nz_d) + (ny_p-1)*nz_d + k+1 ] -By_ext[ i*nz_d + (k+1) ] )
                                            +                   Epsilon_Ymax *( By3D[ i*(ny_p*nz_d) + (ny_p-1)*nz_d + k   ] -By_ext[ i*nz_d +  k    ] )
                                            + Bz_ext[ i*nz_p + k ];
                                            
            }//k  ---end compute Bz
        }//j  ---end compute Bz
    
    } else if( min_max==4 && patch->isZmin() ) {
    
        // for Bx^(p,d,d)
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
            
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + 0 ] = Alpha_Zmin   * Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + 0 ]
                                       +              Beta_Zmin    *( Bx3D[  i   *(ny_d*nz_d) + j*(nz_d) + 1 ]-Bx_ext[  i   *ny_d + j ] )
                                       +              Delta_Zmin   *( Bz3D[ (i+1)*(ny_d*nz_p) + j*(nz_p) + 0 ]-Bz_ext[ (i+1)*ny_d + j ] )
                                       +              Epsilon_Zmin *( Bz3D[  i   *(ny_d*nz_p) + j*(nz_p) + 0 ]-Bz_ext[  i   *ny_d + j ] )
                                       + Bx_ext[ i*ny_d + j ];
            }// j  ---end compute Bx
        }//i  ---end compute Bx
        
        // for By^(d,p,d)
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            
                By3D[ i*(ny_p*nz_d)+ j*(nz_d) + 0 ] = - Alpha_Zmin   * Ex3D[ i*(ny_p*nz_p)+ j*(nz_p) + 0 ]
                                       +              Beta_Zmin   *( By3D[ i*(ny_p*nz_d)+  j   *(nz_d) + 1 ]-By_ext[ i*ny_p +  j   ] )
                                       +              Zeta_Zmin   *( Bz3D[ i*(ny_d*nz_p)+ (j+1)*(nz_p) + 0 ]-Bz_ext[ i*ny_d + (j+1) ] )
                                       +              Eta_Zmin    *( Bz3D[ i*(ny_d*nz_p)+  j   *(nz_p) + 0 ]-Bz_ext[ i*ny_d +  j   ] )
                                       + By_ext[ i*ny_p + j ];
                                       
            }// j  ---end compute By
        }//i  ---end compute By
        
    } else if( min_max==5 && patch->isZmax() ) {
    
        // for Bx^(p,d,d)
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
            
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + (nz_d-1) ] = Alpha_Zmax   * Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + (nz_p-1) ]
                                            +                   Beta_Zmax    *( Bx3D[  i   *(ny_d*nz_d) + j*(nz_d) + (nz_d-2) ] -Bx_ext[  i   *ny_d + j ] )
                                            +                   Delta_Zmax   *( Bz3D[ (i+1)*(ny_d*nz_p) + j*(nz_p) + (nz_p-1) ] -Bz_ext[ (i+1)*ny_d + j ] )
                                            +                   Epsilon_Zmax *( Bz3D[  i   *(ny_d*nz_p) + j*(nz_p) + (nz_p-1) ] -Bz_ext[  i   *ny_d + j ] )
                                            + Bx_ext[ i*ny_d + j ];
                                            
            }//j  ---end compute Bx
        }//i  ---end compute Bx
        
        
        // for By^(d,p,d)
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + (nz_d-1) ] = -Alpha_Zmax * Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + (nz_p-1) ]
                                            +                    Beta_Zmax  *( By3D[ i*(ny_p*nz_d) +  j   *(nz_d) + (nz_d-2) ]-By_ext[ i*ny_p +  j ] )
                                            +                    Zeta_Zmax  *( Bz3D[ i*(ny_d*nz_p) + (j+1)*(nz_p) + (nz_p-1) ]-Bz_ext[ i*ny_d + (j+1) ] )
                                            +                    Eta_Zmax   *( Bz3D[ i*(ny_d*nz_p) +  j   *(nz_p) + (nz_p-1) ]-Bz_ext[ i*ny_d +  j ] )
                                            + By_ext[ i*ny_p + j ];
                                            
            }//j  ---end compute By
        }//i  ---end compute By
        
    }
}
