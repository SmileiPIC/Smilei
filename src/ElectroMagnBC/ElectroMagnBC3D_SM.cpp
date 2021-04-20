#include "ElectroMagnBC3D_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>
#ifdef _GPU
#include <openacc.h>
#endif

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field3D.h"
#include "Tools.h"
#include "Laser.h"

using namespace std;

ElectroMagnBC3D_SM::ElectroMagnBC3D_SM( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC3D( params, patch, i_boundary )
{
    // Calculate axes indices
    axis0_ = i_boundary_ / 2; // axis normal to border
    axis1_ = axis0_ == 0 ? 1 : 0; // other axis (tangent to border)
    axis2_ = axis0_ == 2 ? 1 : 2; // other axis (tangent to border)
    sign_ = (double) (i_boundary_ % 2) *2 - 1.; // -1 or 1 for min or max
    
    // Index where to set the field along axis 0
    iB_.resize( 3 );
    if( sign_ < 0 ) {
        iB_[0] = 0;
        iB_[1] = 0;
        iB_[2] = 0;
    } else {
        iB_[axis0_] = n_p[axis0_] - 1;
        iB_[axis1_] = n_d[axis0_] - 1;
        iB_[axis2_] = n_d[axis0_] - 1;
    }
    
    B_val.resize( 3, nullptr );
    if( patch->isBoundary( i_boundary_ ) ) {
        std::vector<unsigned int> dims0 = { n_d[axis1_], n_d[axis2_] };
        std::vector<unsigned int> dims1 = { n_p[axis1_], n_d[axis2_] };
        std::vector<unsigned int> dims2 = { n_d[axis1_], n_p[axis2_] };
        B_val[axis0_] = new Field2D( dims0, "B_val" );
        B_val[axis1_] = new Field2D( dims1, "B_val" );
        B_val[axis2_] = new Field2D( dims2, "B_val" );
        B_val[0]->put_to( 0. );
        B_val[1]->put_to( 0. );
        B_val[2]->put_to( 0. );
    }
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    vector<double> K = params.EM_BCs_k[i_boundary_];
    double Knorm = sqrt( K[0]*K[0] + K[1]*K[1] + K[2]*K[2] ) ;
    double omega = 1.;
    double k0 = omega*K[axis0_] / Knorm;
    double k1 = omega*K[axis1_] / Knorm;
    double k2 = omega*K[axis2_] / Knorm;
    
    double factor = 1.0 / ( k0 - sign_ * dt_ov_d[axis0_] );
    Alpha_   = 2.0 * factor;
    Beta_    = - ( k0 + sign_ * dt_ov_d[axis0_] ) * factor;
    Gamma_   = 4.0 * k0 * factor;
    Delta_   = - ( k1 + dt_ov_d[axis1_] ) * factor;
    Epsilon_ = - ( k1 - dt_ov_d[axis1_] ) * factor;
    Zeta_    = - ( k2 + dt_ov_d[axis2_] ) * factor;
    Eta_     = - ( k2 - dt_ov_d[axis2_] ) * factor;
}

ElectroMagnBC3D_SM::~ElectroMagnBC3D_SM()
{
    if( B_val[0] ) {
        delete B_val[0] ;
    }
    if( B_val[1] ) {
        delete B_val[1] ;
    }
    if( B_val[2] ) {
        delete B_val[2] ;
    }
}


// Magnetic field Bx^(p,d,d)
// Magnetic field By^(d,p,d)
// Magnetic field Bz^(d,d,p)

void ElectroMagnBC3D_SM::save_fields( Field *my_field, Patch *patch )
{
    Field3D *field3D=static_cast<Field3D *>( my_field );
    
    if( patch->isBoundary( i_boundary_ ) ) {
        
        unsigned int xyz = 0;
        if( field3D->name=="Bx" ) {
            xyz = 0;
        } else if( field3D->name=="By" ) {
            xyz = 1;
        } else if( field3D->name=="Bz" ) {
            xyz = 2;
        }
        
        if( axis0_ == 0 ) {
            field3D->extract_slice_yz( iB_[xyz], B_val[xyz] );
        } else if( axis0_ == 1 ) {
            field3D->extract_slice_xz( iB_[xyz], B_val[xyz] );
        } else {
            field3D->extract_slice_xy( iB_[xyz], B_val[xyz] );
        }
        
    }
}


void ElectroMagnBC3D_SM::disableExternalFields()
{
    delete B_val[0];
    B_val[0] = NULL;
    delete B_val[1];
    B_val[1] = NULL;
    delete B_val[2];
    B_val[2] = NULL;
}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_SM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
<<<<<<< HEAD
    // Static cast of the fields
    double *Ex3D = &(EMfields->Ex_->data_[0]);
    double *Ey3D = &(EMfields->Ey_->data_[0]);
    double *Ez3D = &(EMfields->Ez_->data_[0]);
    double *Bx3D = &(EMfields->Bx_->data_[0]);
    double *By3D = &(EMfields->By_->data_[0]);
    double *Bz3D = &(EMfields->Bz_->data_[0]);
    int sizeofEx = EMfields->Ex_->globalDims_;
    int sizeofEy = EMfields->Ey_->globalDims_;
    int sizeofEz = EMfields->Ez_->globalDims_;
    int sizeofBx = EMfields->Bx_->globalDims_;
    int sizeofBy = EMfields->By_->globalDims_;
    int sizeofBz = EMfields->Bz_->globalDims_;
    int byw_size = ny_p*nz_d;
    int bzw_size = ny_d*nz_p;

    vector<double> pos( 2 );

    vector<double> byW( ny_p*nz_d, 0. );
    vector<double> byE( ny_p*nz_d, 0. );
    vector<double> bzW( ny_d*nz_p, 0. );
    vector<double> bzE( ny_d*nz_p, 0. );
    double* dbyW = &(byW[0]);
    double* dbyE = &(byE[0]);
    double* dbzW = &(bzW[0]);
    double* dbzE = &(bzE[0]);
    int isymin = patch->isYmin();
    int isymax = patch->isYmax();
    int iszmin = patch->isZmin();
    int iszmax = patch->isZmax();

    int Bx_ext_size, By_ext_size, Bz_ext_size;
    double* Bx_ext = NULL;
    if (Bx_val!=nullptr) {
        Bx_ext = &(Bx_val->data_[0]);
        Bx_ext_size = Bx_val->globalDims_;
    }
    double* By_ext = NULL;
    if (By_val!=nullptr) {
        By_ext = &(By_val->data_[0]);
        By_ext_size = By_val->globalDims_;
    }
    double* Bz_ext = NULL;
    if (Bz_val!=nullptr) {
        Bz_ext = &(Bz_val->data_[0]);
        Bz_ext_size = Bz_val->globalDims_;
    }
#ifdef _GPU
    if ( (Bx_ext!=NULL) && (!acc_deviceptr( Bx_ext )) ) {
        #pragma acc enter data copyin(Bx_ext[0:Bx_ext_size])
    }
    if ( (By_ext!=NULL) && (!acc_deviceptr( By_ext )) ) {
        #pragma acc enter data copyin(By_ext[0:By_ext_size])
    }
    if ( (Bz_ext!=NULL) && !acc_deviceptr( Bz_ext ) ) {
        #pragma acc enter data copyin(Bz_ext[0:Bz_ext_size])
    }
#endif

    if( min_max==0 && patch->isXmin() ) {
    
        // for By^(d,p,d)
        for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*dy;
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k -0.5 - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                byW[ j*nz_d+k ] = 0.;
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    byW[ j*nz_d+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, k );
                }
            }
        }
#ifdef _GPU
        #pragma acc parallel present(Ez3D[0:sizeofEz],Bx3D[0:sizeofBx],By3D[0:sizeofBy],By_ext[0:By_ext_size],Bx_ext[0:Bx_ext_size]) copyin(dbyW[0:byw_size])
        #pragma acc loop gang
#endif
        for( unsigned int j=isymin ; j<ny_p-isymax ; j++ ) {
#ifdef _GPU
            #pragma acc loop worker vector
#endif
            for( unsigned int k=iszmin ; k<nz_d-iszmax ; k++ ) {
                By3D[ 0*(ny_p*nz_d) + j*nz_d + k ] = Alpha_SM_W   * Ez3D[ 0*(ny_p*nz_d) + j*nz_d + k ]
                                       +              Beta_SM_W    *( By3D[ 1*(ny_p*nz_d) + j*nz_d + k ]-By_ext[ j*nz_d + k ] )
                                       +              Gamma_SM_W   * dbyW[ j*nz_d+k ]
                                       +              Delta_SM_W   *( Bx3D[ 0*(ny_d*nz_d) + (j+1)*nz_d + k ]-Bx_ext[ (j+1)*nz_d + k ] )
                                       +              Epsilon_SM_W *( Bx3D[ 0*(ny_d*nz_d) +  j   *nz_d + k ]-Bx_ext[  j   *nz_d + k ] )
                                       + By_ext[ j*nz_d + k ];
            }// k  ---end compute By
        }//j  ---end compute By

        // for Bz^(d,d,p)
        for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*dy;
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                bzW[ j*nz_p+k ] = 0.;
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bzW[ j*nz_p+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, k );
                }
            }
        }
#ifdef _GPU
        #pragma acc parallel present(Ey3D[0:sizeofEy],Bx3D[0:sizeofBx],Bz3D[0:sizeofBz],Bz_ext[0:Bz_ext_size],Bx_ext[0:Bx_ext_size]) copyin(dbzW[0:bzw_size])
        #pragma acc loop gang
#endif
        for( unsigned int j=isymin ; j<ny_d-isymax ; j++ ) {
#ifdef _GPU
            #pragma acc loop worker vector
#endif
            for( unsigned int k=iszmin ; k<nz_p-iszmax ; k++ ) {
                Bz3D[ 0*(ny_d*nz_p) + j*nz_p + k ] = - Alpha_SM_W   * Ey3D[ 0*(ny_d*nz_p) + j*nz_p + k ]
                                       +              Beta_SM_W    *( Bz3D[ 1*(ny_d*nz_p) + j*nz_p + k ]-Bz_ext[ j*nz_p + k ] )
                                       +              Gamma_SM_W   * dbzW[ j*nz_p+k ]
                                       +              Zeta_SM_W    *( Bx3D[ 0*(ny_d*nz_d) + j*nz_d + k+1 ]-Bx_ext[ j*nz_d + (k+1) ] )
                                       +              Eta_SM_W     *( Bx3D[ 0*(ny_d*nz_d) + j*nz_d + k   ]-Bx_ext[ j*nz_d +  k    ] )
                                       + Bz_ext[ j*nz_p + k ];
                                       
            }// k  ---end compute Bz
        }//j  ---end compute Bz
    } else if( min_max==1 && patch->isXmax() ) {
    
        // for By^(d,p,d)
        for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*dy;
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - 0.5 - ( int )EMfields->oversize[2] )*dz;
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    byE[ j*nz_d+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, k );
                }
            }
        }
#ifdef _GPU
        #pragma acc parallel present(Ez3D[0:sizeofEz],Bx3D[0:sizeofBx],By3D[0:sizeofBy],By_ext[0:By_ext_size],Bx_ext[0:Bx_ext_size]) copyin(dbyE[0:byw_size])
        #pragma acc loop gang
#endif
        for( unsigned int j=isymin ; j<ny_p-isymax ; j++ ) {
#ifdef _GPU
            #pragma acc loop worker vector
#endif
            for( unsigned int k=iszmin ; k<nz_d-iszmax ; k++ ) {
                By3D[ (nx_d-1)*(ny_p*nz_d) + j*nz_d + k ] = Alpha_SM_E   * Ez3D[ (nx_p-1)*(ny_p*nz_d) + j*nz_d + k ]
                                            +                   Beta_SM_E    *( By3D[ (nx_d-2)*(ny_p*nz_d) + j*nz_d + k ] -By_ext[ j*nz_d + k ] )
                                            +                   Gamma_SM_E   * dbyE[ j*nz_d+k ]
                                            +                   Delta_SM_E   *( Bx3D[ (nx_p-1)*(ny_d*nz_d) + (j+1)*nz_d + k ] -Bx_ext[ (j+1)*nz_d + k ] ) // Check x-index
                                            +                   Epsilon_SM_E *( Bx3D[ (nx_p-1)*(ny_d*nz_d) +  j   *nz_d + k ] -Bx_ext[  j   *nz_d + k ] )
                                            + By_ext[ j*nz_d + k ];
                                            
            }//k  ---end compute By
        }//j  ---end compute By
=======

    if( patch->isBoundary( i_boundary_ ) ) {
        
        // Static cast of the fields
        vector<double*> E( 3 );
        E[0] = &( EMfields->Ex_->data_[0] );
        E[1] = &( EMfields->Ey_->data_[0] );
        E[2] = &( EMfields->Ez_->data_[0] );
        vector<double*> B( 3 );
        B[0] = &( EMfields->Bx_->data_[0] );
        B[1] = &( EMfields->By_->data_[0] );
        B[2] = &( EMfields->Bz_->data_[0] );
        
        vector<double*> B_ext = { NULL, NULL, NULL };
        if( B_val[0] ) { B_ext[0] = &(B_val[0]->data_[0]); }
        if( B_val[1] ) { B_ext[1] = &(B_val[1]->data_[0]); }
        if( B_val[2] ) { B_ext[2] = &(B_val[2]->data_[0]); }
        
        vector<double> pos( 2 );
        
        unsigned int nz_p = n_p[2];
        unsigned int nz_d = n_d[2];
        unsigned int nyz_pp = n_p[1]*n_p[2];
        unsigned int nyz_pd = n_p[1]*n_d[2];
        unsigned int nyz_dp = n_d[1]*n_p[2];
        unsigned int nyz_dd = n_d[1]*n_d[2];
        unsigned int n1 = n_p[axis1_];
        unsigned int n2 = n_d[axis2_];
        unsigned int p0 = iB_[axis0_];
        unsigned int p1 = iB_[axis1_] - sign_;
        
        // Component along axis 1
        // Lasers
        vector<double> b1( n1*n2, 0. );
        if( ! vecLaser.empty() ) {
            for( unsigned int j=patch->isBoundary(axis1_,0); j<n1-patch->isBoundary(axis1_,1) ; j++ ) {
                pos[0] = patch->getDomainLocalMin( axis1_ ) + ( ( int )j - ( int )EMfields->oversize[axis1_] )*d[axis1_];
                for( unsigned int k=patch->isBoundary(axis2_,0) ; k<n2-patch->isBoundary(axis2_,1) ; k++ ) {
                    pos[1] = patch->getDomainLocalMin( axis2_ ) + ( ( int )k -0.5 - ( int )EMfields->oversize[axis2_] )*d[axis2_];
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        b1[ j*n2+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, k );
                    }
                }
            }
        }
        // B1
        if( axis0_ == 0 ) {
            for( unsigned int j=patch->isBoundary(axis1_,0); j<n1-patch->isBoundary(axis1_,1) ; j++ ) {
                for( unsigned int k=patch->isBoundary(axis2_,0) ; k<n2-patch->isBoundary(axis2_,1) ; k++ ) {
                    B[axis1_][ iB_[axis1_]*nyz_pd + j*nz_d + k ]
                        = Alpha_   *  E[axis2_][ p0*nyz_pd + j*nz_d + k ]
                        + Beta_    *( B[axis1_][ p1*nyz_pd + j*nz_d + k ]-B_ext[axis1_][ j*n2 + k ] )
                        + Gamma_   * b1[ j*n2 + k ]
                        + Delta_   *( B[axis0_][ p0*nyz_dd + (j+1)*nz_d + k ]-B_ext[axis0_][ (j+1)*nz_d + k ] )
                        + Epsilon_ *( B[axis0_][ p0*nyz_dd +  j   *nz_d + k ]-B_ext[axis0_][  j   *nz_d + k ] )
                        + B_ext[axis1_][ j*n2 + k ];
                }
            }
        } else if( axis0_ == 1 ) {
            for( unsigned int i=patch->isBoundary(axis1_,0); i<n1-patch->isBoundary(axis1_,1) ; i++ ) {
                for( unsigned int k=patch->isBoundary(axis2_,0) ; k<n2-patch->isBoundary(axis2_,1) ; k++ ) {
                    B[axis1_][ i*nyz_dd + iB_[axis1_]*nz_d + k ]
                        =-Alpha_   *  E[axis2_][ i*nyz_pd + p0*nz_d + k ]
                        + Beta_    *( B[axis1_][ i*nyz_dd + p1*nz_d + k ]-B_ext[axis1_][ i*n2 + k ] )
                        + Gamma_   * b1[ i*n2 + k ]
                        + Delta_   *( B[axis0_][ (i+1)*nyz_pd + p0*nz_d + k ]-B_ext[axis0_][ (i+1)*nz_d + k ] )
                        + Epsilon_ *( B[axis0_][  i   *nyz_pd + p0*nz_d + k ]-B_ext[axis0_][  i   *nz_d + k ] )
                        + B_ext[axis1_][ i*n2 + k ];
                }
            }
        } else {
            for( unsigned int i=patch->isBoundary(axis1_,0); i<n1-patch->isBoundary(axis1_,1) ; i++ ) {
                for( unsigned int j=patch->isBoundary(axis2_,0) ; j<n2-patch->isBoundary(axis2_,1) ; j++ ) {
                    B[axis1_][ i*nyz_dd + j*nz_d + iB_[axis1_] ]
                        = Alpha_   *  E[axis2_][ i*nyz_dp + j*nz_p + p0 ]
                        + Beta_    *( B[axis1_][ i*nyz_dd + j*nz_d + p1 ]-B_ext[axis1_][ i*n2 + j ] )
                        + Gamma_   * b1[ i*n2 + j ]
                        + Delta_   *( B[axis0_][ (i+1)*nyz_dp + j*nz_p + p0 ]-B_ext[axis0_][ (i+1)*n2 + j ] )
                        + Epsilon_ *( B[axis0_][  i   *nyz_dp + j*nz_p + p0 ]-B_ext[axis0_][  i   *n2 + j ] )
                        + B_ext[axis1_][ i*n2 + j ];
                }
            }
        }
>>>>>>> develop
        
        // Component along axis 2
        // Lasers
        unsigned int n1d = n_d[axis1_];
        unsigned int n2p = n_p[axis2_];
        vector<double> b2( n1d*n2p, 0. );
        if( ! vecLaser.empty() ) {
            for( unsigned int j=patch->isBoundary(axis1_,0); j<n1d-patch->isBoundary(axis1_,1) ; j++ ) {
                pos[0] = patch->getDomainLocalMin( axis1_ ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[axis1_] )*d[axis1_];
                for( unsigned int k=patch->isBoundary(axis2_,0) ; k<n2p-patch->isBoundary(axis2_,1) ; k++ ) {
                    pos[1] = patch->getDomainLocalMin( axis2_ ) + ( ( int )k - ( int )EMfields->oversize[axis2_] )*d[axis2_];
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        b2[ j*n2p+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, k );
                    }
                }
            }
        }
        // B2
        if( axis0_ == 0 ) {
            for( unsigned int j=patch->isBoundary(axis1_,0); j<n1d-patch->isBoundary(axis1_,1) ; j++ ) {
                for( unsigned int k=patch->isBoundary(axis2_,0) ; k<n2p-patch->isBoundary(axis2_,1) ; k++ ) {
                    B[axis2_][ iB_[axis1_]*nyz_dp + j*nz_p + k ]
                        = -Alpha_ *  E[axis1_][ p0*nyz_dp + j*nz_p + k ]
                        +  Beta_  *( B[axis2_][ p1*nyz_dp + j*nz_p + k ]-B_ext[axis2_][ j*n2p + k ] )
                        +  Gamma_ * b2[ j*n2p + k ]
                        +  Zeta_  *( B[axis0_][ p0*nyz_dd + j*nz_d + k+1 ]-B_ext[axis0_][ j*n2 + (k+1) ] )
                        +  Eta_   *( B[axis0_][ p0*nyz_dd + j*nz_d + k   ]-B_ext[axis0_][ j*n2 +  k    ] )
                        +  B_ext[axis2_][ j*n2p + k ];
                }
            }
        } else if( axis0_ == 1 ) {
            for( unsigned int i=patch->isBoundary(axis1_,0); i<n1d-patch->isBoundary(axis1_,1) ; i++ ) {
                for( unsigned int k=patch->isBoundary(axis2_,0) ; k<n2p-patch->isBoundary(axis2_,1) ; k++ ) {
                    B[axis2_][ i*nyz_dp + iB_[axis1_]*nz_p + k ]
                        =  Alpha_ *  E[axis1_][ i*nyz_pp + p0*nz_p + k ]
                        +  Beta_  *( B[axis2_][ i*nyz_dp + p1*nz_p + k ]-B_ext[axis2_][ i*n2p + k ] )
                        +  Gamma_ * b2[ i*n2p + k ]
                        +  Zeta_  *( B[axis0_][ i*nyz_pd + p0*nz_d + k+1 ]-B_ext[axis0_][ i*n2 + (k+1) ] )
                        +  Eta_   *( B[axis0_][ i*nyz_pd + p0*nz_d + k   ]-B_ext[axis0_][ i*n2 +  k    ] )
                        +  B_ext[axis2_][ i*n2p + k ];
                }
            }
        } else {
            for( unsigned int i=patch->isBoundary(axis1_,0); i<n1d-patch->isBoundary(axis1_,1) ; i++ ) {
                for( unsigned int j=patch->isBoundary(axis2_,0) ; j<n2p-patch->isBoundary(axis2_,1) ; j++ ) {
                    B[axis2_][ i*nyz_pd + j*nz_d + iB_[axis1_] ]
                        = -Alpha_ *  E[axis1_][ i*nyz_pp + j*nz_p + p0 ]
                        +  Beta_  *( B[axis2_][ i*nyz_pd + j*nz_d  + p1 ]-B_ext[axis2_][ i*n2p + j ] )
                        +  Gamma_ * b2[ i*n2p + j ]
                        +  Zeta_  *( B[axis0_][ i*nyz_dp + (j+1)*nz_p + p0 ]-B_ext[axis0_][ i*n2 + (j+1) ] )
                        +  Eta_   *( B[axis0_][ i*nyz_dp +  j   *nz_p + p0 ]-B_ext[axis0_][ i*n2 +  j    ] )
                        +  B_ext[axis2_][ i*n2p + j ];
                }
            }
        }
<<<<<<< HEAD
#ifdef _GPU
        #pragma acc parallel present(Ey3D[0:sizeofEy],Bx3D[0:sizeofBx],Bz3D[0:sizeofBz],Bz_ext[0:Bz_ext_size],Bx_ext[0:Bx_ext_size]) copyin(dbzE[0:bzw_size])
        #pragma acc loop gang
#endif
        for( unsigned int j=isymin ; j<ny_d-isymax; j++ ) {
#ifdef _GPU
            #pragma acc loop worker vector
#endif
            for( unsigned int k=iszmin ; k<nz_p-iszmax ; k++ ) {
                Bz3D[ (nx_d-1)*(ny_d*nz_p) + j*nz_p + k ] = -Alpha_SM_E * Ey3D[ (nx_p-1)*(ny_d*nz_p) + j*nz_p + k ]
                                            +                    Beta_SM_E  *( Bz3D[ (nx_d-2)*(ny_d*nz_p) + j*nz_p + k ] -Bz_ext[ j*nz_p + k ] )
                                            +                    Gamma_SM_E * dbzE[ j*nz_p+k ]
                                            +                    Zeta_SM_E  *( Bx3D[ (nx_p-1)*(ny_d*nz_d) + j*nz_d + k+1 ]-Bx_ext[ j*nz_d + (k+1) ] )
                                            +                    Eta_SM_E   *( Bx3D[ (nx_p-1)*(ny_d*nz_d) + j*nz_d + k   ]-Bx_ext[ j*nz_d +  k    ] )
                                            + Bz_ext[ j*nz_p + k ];
            }//k  ---end compute Bz
        }//j  ---end compute Bz
    } else if( min_max==2 && patch->isYmin() ) {
    
        // for Bx^(p,d,d)
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                Bx3D[ i*(ny_d*nz_d) + 0*nz_d + k ] = - Alpha_SM_S   * Ez3D[ i*(ny_p*nz_d) + 0*nz_d + k ]
                                       +              Beta_SM_S     *( Bx3D[  i   *(ny_d*nz_d) + 1*nz_d + k ]-Bx_ext[  i   *nz_d + k ] )
                                       +              Zeta_SM_S     *( By3D[ (i+1)*(ny_p*nz_d) + 0*nz_d + k ]-By_ext[ (i+1)*nz_d + k ] )
                                       +              Eta_SM_S      *( By3D[  i   *(ny_p*nz_d) + 0*nz_d + k ]-By_ext[  i   *nz_d + k ] )
                                       + Bx_ext[ i*nz_d + k ];
            }// k  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bz^(d,d,p)
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
                Bz3D[ i*(ny_d*nz_p) + 0*nz_p + k ] = Alpha_SM_S   * Ex3D[ i*(ny_p*nz_p) + 0*nz_p + k ]
                                       +              Beta_SM_S    *( Bz3D[ i*(ny_d*nz_p) + 1*nz_p + k   ]-Bz_ext[ i*nz_p +  k    ] )
                                       +              Delta_SM_S   *( By3D[ i*(ny_p*nz_d) + 0*nz_d + k+1 ]-By_ext[ i*nz_d + (k+1) ] )
                                       +              Epsilon_SM_S *( By3D[ i*(ny_p*nz_d) + 0*nz_d + k   ]-By_ext[ i*nz_d +  k    ] )
                                       + Bz_ext[ i*nz_p + k ];
            }// k  ---end compute Bz
        }//i  ---end compute Bz       }
    } else if( min_max==3 && patch->isYmax() ) {
    
        // for Bx^(p,d,d)
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
            
                Bx3D[ i*(ny_d*nz_d) + (ny_d-1)*nz_d + k ] = -Alpha_SM_N * Ez3D[ i*(ny_p*nz_d) + (ny_p-1)*nz_d + k ]
                                            +                    Beta_SM_N  *( Bx3D[  i   *(ny_d*nz_d) + (ny_d-2)*nz_d + k ]-Bx_ext[  i   *nz_d + k ] )
                                            +                    Zeta_SM_N  *( By3D[ (i+1)*(ny_p*nz_d) + (ny_p-1)*nz_d + k ]-By_ext[ (i+1)*nz_d + k ] )
                                            +                    Eta_SM_N   *( By3D[  i   *(ny_p*nz_d) + (ny_p-1)*nz_d + k ]-By_ext[  i   *nz_d + k ] )
                                            + Bx_ext[ i*nz_d + k ];
                                            
            }//k  ---end compute Bz
        }//j  ---end compute Bz
        
        // for Bz^(d,d,p)
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            for( unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++ ) {
            
                Bz3D[ i*(ny_d*nz_p) + (ny_d-1)*nz_p + k ] = Alpha_SM_N   * Ex3D[ i*(ny_p*nz_p) + (ny_p-1)*nz_p + k ]
                                            +                   Beta_SM_N    *( Bz3D[ i*(ny_d*nz_p) + (ny_d-2)*nz_p + k   ] -Bz_ext[ i*nz_p +  k    ] )
                                            +                   Delta_SM_N   *( By3D[ i*(ny_p*nz_d) + (ny_p-1)*nz_d + k+1 ] -By_ext[ i*nz_d + (k+1) ] )
                                            +                   Epsilon_SM_N *( By3D[ i*(ny_p*nz_d) + (ny_p-1)*nz_d + k   ] -By_ext[ i*nz_d +  k    ] )
                                            + Bz_ext[ i*nz_p + k ];
                                            
            }//k  ---end compute Bz
        }//j  ---end compute Bz
    } else if( min_max==4 && patch->isZmin() ) {
    
        // for Bx^(p,d,d)
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
            
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + 0 ] = Alpha_SM_B   * Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + 0 ]
                                       +              Beta_SM_B    *( Bx3D[  i   *(ny_d*nz_d) + j*(nz_d) + 1 ]-Bx_ext[  i   *ny_d + j ] )
                                       +              Delta_SM_B   *( Bz3D[ (i+1)*(ny_d*nz_p) + j*(nz_p) + 0 ]-Bz_ext[ (i+1)*ny_d + j ] )
                                       +              Epsilon_SM_B *( Bz3D[  i   *(ny_d*nz_p) + j*(nz_p) + 0 ]-Bz_ext[  i   *ny_d + j ] )
                                       + Bx_ext[ i*ny_d + j ];
            }// j  ---end compute Bx
        }//i  ---end compute Bx
        
        // for By^(d,p,d)
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            
                By3D[ i*(ny_p*nz_d)+ j*(nz_d) + 0 ] = - Alpha_SM_B   * Ex3D[ i*(ny_p*nz_p)+ j*(nz_p) + 0 ]
                                       +              Beta_SM_B   *( By3D[ i*(ny_p*nz_d)+  j   *(nz_d) + 1 ]-By_ext[ i*ny_p +  j   ] )
                                       +              Zeta_SM_B   *( Bz3D[ i*(ny_d*nz_p)+ (j+1)*(nz_p) + 0 ]-Bz_ext[ i*ny_d + (j+1) ] )
                                       +              Eta_SM_B    *( Bz3D[ i*(ny_d*nz_p)+  j   *(nz_p) + 0 ]-Bz_ext[ i*ny_d +  j   ] )
                                       + By_ext[ i*ny_p + j ];
                                       
            }// j  ---end compute By
        }//i  ---end compute By
        
    } else if( min_max==5 && patch->isZmax() ) {
    
        // for Bx^(p,d,d)
        for( unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++ ) {
            for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
            
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + (nz_d-1) ] = Alpha_SM_T   * Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + (nz_p-1) ]
                                            +                   Beta_SM_T    *( Bx3D[  i   *(ny_d*nz_d) + j*(nz_d) + (nz_d-2) ] -Bx_ext[  i   *ny_d + j ] )
                                            +                   Delta_SM_T   *( Bz3D[ (i+1)*(ny_d*nz_p) + j*(nz_p) + (nz_p-1) ] -Bz_ext[ (i+1)*ny_d + j ] )
                                            +                   Epsilon_SM_T *( Bz3D[  i   *(ny_d*nz_p) + j*(nz_p) + (nz_p-1) ] -Bz_ext[  i   *ny_d + j ] )
                                            + Bx_ext[ i*ny_d + j ];
                                            
            }//j  ---end compute Bx
        }//i  ---end compute Bx
        
        
        // for By^(d,p,d)
        for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
            for( unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++ ) {
            
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + (nz_d-1) ] = -Alpha_SM_T * Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + (nz_p-1) ]
                                            +                    Beta_SM_T  *( By3D[ i*(ny_p*nz_d) +  j   *(nz_d) + (nz_d-2) ]-By_ext[ i*ny_p +  j ] )
                                            +                    Zeta_SM_T  *( Bz3D[ i*(ny_d*nz_p) + (j+1)*(nz_p) + (nz_p-1) ]-Bz_ext[ i*ny_d + (j+1) ] )
                                            +                    Eta_SM_T   *( Bz3D[ i*(ny_d*nz_p) +  j   *(nz_p) + (nz_p-1) ]-Bz_ext[ i*ny_d +  j ] )
                                            + By_ext[ i*ny_p + j ];
                                            
            }//j  ---end compute By
        }//i  ---end compute By
        
=======
>>>>>>> develop
    }
}
