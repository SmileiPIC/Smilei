#include "ElectroMagnBC3D_SM.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include "ElectroMagn.h"
#include "Field3D.h"
#include "Laser.h"
#include "Params.h"
#include "Patch.h"
#include "Tools.h"
#include "gpu.h"

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
    
    std::vector<double> K = params.EM_BCs_k[i_boundary_];
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
    if( patch->isBoundary( i_boundary_ ) ) {

        // Static cast of the fields
        const Field  *E[3]{ EMfields->Ex_, EMfields->Ey_, EMfields->Ez_ };
        const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };

        // double *E0 = E[axis0_]->data_;
        const double *const __restrict__ E1 = E[axis1_]->data_;
        const double *const __restrict__ E2 = E[axis2_]->data_;
        const double *const __restrict__ B0 = B[axis0_]->data_;
        double *const __restrict__ B1       = B[axis1_]->data_;
        double *const __restrict__ B2       = B[axis2_]->data_;

        const double *const __restrict__ B_ext0 = B_val[axis0_]->data_;
        const double *const __restrict__ B_ext1 = B_val[axis1_]->data_;
        const double *const __restrict__ B_ext2 = B_val[axis2_]->data_;

        const unsigned int nz_p   = n_p[2];
        const unsigned int nz_d   = n_d[2];
        const unsigned int nyz_pp = n_p[1] * n_p[2];
        const unsigned int nyz_pd = n_p[1] * n_d[2];
        const unsigned int nyz_dp = n_d[1] * n_p[2];
        const unsigned int nyz_dd = n_d[1] * n_d[2];
        const unsigned int n1p    = n_p[axis1_];
        const unsigned int n1d    = n_d[axis1_];
        const unsigned int n2p    = n_p[axis2_];
        const unsigned int n2d    = n_d[axis2_];
        const unsigned int p0     = iB_[axis0_];
        const unsigned int p1     = iB_[axis1_] - sign_;
        const unsigned int iB1    = iB_[axis1_];

        const int b1_size = n1p * n2d;
        const int b2_size = n1d * n2p;

        std::vector<double> b1( b1_size, 0. );
        std::vector<double> b2( b2_size, 0. );
        std::vector<double> pos( 2 );

        double *const __restrict__ db1 = b1.data();
        double *const __restrict__ db2 = b2.data();

        const int isBoundary1min = patch->isBoundary( axis1_, 0 );
        const int isBoundary1max = patch->isBoundary( axis1_, 1 );
        const int isBoundary2min = patch->isBoundary( axis2_, 0 );
        const int isBoundary2max = patch->isBoundary( axis2_, 1 );

#ifdef _GPU
        const int sizeofE0 = E[axis0_]->globalDims_;
        const int sizeofE1 = E[axis1_]->globalDims_;
        const int sizeofE2 = E[axis2_]->globalDims_;
        const int sizeofB0 = B[axis0_]->globalDims_;
        const int sizeofB1 = B[axis1_]->globalDims_;
        const int sizeofB2 = B[axis2_]->globalDims_;
#endif

        const int B_ext_size0 = B_val[axis0_]->globalDims_;
        const int B_ext_size1 = B_val[axis1_]->globalDims_;
        const int B_ext_size2 = B_val[axis2_]->globalDims_;

        smilei::tools::gpu::HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( B_ext0, B_ext_size0 );
        smilei::tools::gpu::HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( B_ext1, B_ext_size1 );
        smilei::tools::gpu::HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( B_ext2, B_ext_size2 );

        if( vecLaser.empty() ) {
            smilei::tools::gpu::HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( db1, b1_size );
            smilei::tools::gpu::HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( db2, b2_size );
        }

        // Component along axis 1
        // Lasers
        if( !vecLaser.empty() ) {
            for( unsigned int j=isBoundary1min; j<n1p-isBoundary1max ; j++ ) {
                pos[0] = patch->getDomainLocalMin( axis1_ ) + ( ( int )j - ( int )EMfields->oversize[axis1_] )*d[axis1_];
                for( unsigned int k=isBoundary2min; k<n2d-isBoundary2max; k++ ) {
                    pos[1] = patch->getDomainLocalMin( axis2_ ) + ( ( int )k -0.5 - ( int )EMfields->oversize[axis2_] )*d[axis2_];
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        db1[ j*n2d+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, k );
                    }
                }
            }

            smilei::tools::gpu::HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( db1, b1_size );
        }
        // B1
        if( axis0_ == 0 ) {
#ifdef _GPU
            #pragma acc parallel present(E2[0:sizeofE2],B0[0:sizeofB0],B1[0:sizeofB1],B_ext1[0:B_ext_size1],B_ext0[0:B_ext_size0]) copyin(db1[0:b1_size])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams
            #pragma omp distribute parallel for collapse( 2 )
#endif
            for( unsigned int j=isBoundary1min; j<n1p-isBoundary1max ; j++ ) {
#ifdef _GPU
            #pragma acc loop worker vector
#endif
                for( unsigned int k=isBoundary2min ; k<n2d-isBoundary2max ; k++ ) {
                    B1[ iB1*nyz_pd + j*nz_d + k ]
                        = Alpha_   *  E2[ p0*nyz_pd + j*nz_d + k ]
                        + Beta_    *( B1[ p1*nyz_pd + j*nz_d + k ]-B_ext1[ j*n2d + k ] )
                        + Gamma_   * db1[ j*n2d + k ]
                        + Delta_   *( B0[ p0*nyz_dd + (j+1)*nz_d + k ]-B_ext0[ (j+1)*nz_d + k ] )
                        + Epsilon_ *( B0[ p0*nyz_dd +  j   *nz_d + k ]-B_ext0[  j   *nz_d + k ] )
                        + B_ext1[ j*n2d + k ];
                }
            }
        } else if( axis0_ == 1 ) {
#ifdef _GPU
            #pragma acc parallel present(E2[0:sizeofE2],B0[0:sizeofB0],B1[0:sizeofB1],B_ext1[0:B_ext_size1],B_ext0[0:B_ext_size0]) copyin(db1[0:b1_size])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams
            #pragma omp distribute parallel for collapse( 2 )
#endif
            for( unsigned int i=isBoundary1min; i<n1p-isBoundary1max ; i++ ) {
#ifdef _GPU
            #pragma acc loop worker vector
#endif
                for( unsigned int k=isBoundary2min ; k<n2d-isBoundary2max ; k++ ) {
                    B1[ i*nyz_dd + iB1*nz_d + k ]
                        =-Alpha_   *  E2[ i*nyz_pd + p0*nz_d + k ]
                        + Beta_    *( B1[ i*nyz_dd + p1*nz_d + k ]-B_ext1[ i*n2d + k ] )
                        + Gamma_   * db1[ i*n2d + k ]
                        + Delta_   *( B0[ (i+1)*nyz_pd + p0*nz_d + k ]-B_ext0[ (i+1)*nz_d + k ] )
                        + Epsilon_ *( B0[  i   *nyz_pd + p0*nz_d + k ]-B_ext0[  i   *nz_d + k ] )
                        + B_ext1[ i*n2d + k ];
                }
            }
        } else {
#ifdef _GPU
            #pragma acc parallel present(E2[0:sizeofE2],B0[0:sizeofB0],B1[0:sizeofB1],B_ext1[0:B_ext_size1],B_ext0[0:B_ext_size0]) copyin(db1[0:b1_size])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams
            #pragma omp distribute parallel for collapse( 2 )
#endif
            for( unsigned int i=isBoundary1min; i<n1p-isBoundary1max ; i++ ) {
#ifdef _GPU
            #pragma acc loop worker vector
#endif
                for( unsigned int j=isBoundary2min ; j<n2d-isBoundary2max ; j++ ) {
                    B1[ i*nyz_dd + j*nz_d + iB1 ]
                        = Alpha_   *  E2[ i*nyz_dp + j*nz_p + p0 ]
                        + Beta_    *( B1[ i*nyz_dd + j*nz_d + p1 ]-B_ext1[ i*n2d + j ] )
                        + Gamma_   * db1[ i*n2d + j ]
                        + Delta_   *( B0[ (i+1)*nyz_dp + j*nz_p + p0 ]-B_ext0[ (i+1)*n2d + j ] )
                        + Epsilon_ *( B0[  i   *nyz_dp + j*nz_p + p0 ]-B_ext0[  i   *n2d + j ] )
                        + B_ext1[ i*n2d + j ];
                }
            }
        }
        
        // Component along axis 2
        // Lasers
        if( !vecLaser.empty() ) {
            for( unsigned int j=isBoundary1min; j<n1d-isBoundary1max; j++ ) {
                pos[0] = patch->getDomainLocalMin( axis1_ ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[axis1_] )*d[axis1_];
                for( unsigned int k=isBoundary2min; k<n2p-isBoundary2max; k++ ) {
                    pos[1] = patch->getDomainLocalMin( axis2_ ) + ( ( int )k - ( int )EMfields->oversize[axis2_] )*d[axis2_];
                    for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                        db2[ j*n2p+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, k );
                    }
                }
            }

            smilei::tools::gpu::HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( db2, b2_size );
        }
        // B2
        if( axis0_ == 0 ) {
#ifdef _GPU
            #pragma acc parallel present(E1[0:sizeofE1],B0[0:sizeofB0],B2[0:sizeofB2],B_ext2[0:B_ext_size2],B_ext0[0:B_ext_size0]) copyin(db2[0:b2_size])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams
            #pragma omp distribute parallel for collapse( 2 )
#endif
            for( unsigned int j=isBoundary1min; j<n1d-isBoundary1max ; j++ ) {
#ifdef _GPU
                #pragma acc loop worker vector
#endif
                for( unsigned int k=isBoundary2min; k<n2p-isBoundary2max ; k++ ) {
                    B2[ iB1*nyz_dp + j*nz_p + k ]
                        = -Alpha_ *  E1[ p0*nyz_dp + j*nz_p + k ]
                        +  Beta_  *( B2[ p1*nyz_dp + j*nz_p + k ]-B_ext2[ j*n2p + k ] )
                        +  Gamma_ * db2[ j*n2p + k ]
                        +  Zeta_  *( B0[ p0*nyz_dd + j*nz_d + k+1 ]-B_ext0[ j*n2d + (k+1) ] )
                        +  Eta_   *( B0[ p0*nyz_dd + j*nz_d + k   ]-B_ext0[ j*n2d +  k    ] )
                        +  B_ext2[ j*n2p + k ];
                }
            }
        } else if( axis0_ == 1 ) {
#ifdef _GPU
            #pragma acc parallel present(E1[0:sizeofE1],B0[0:sizeofB0],B2[0:sizeofB2],B_ext2[0:B_ext_size2],B_ext0[0:B_ext_size0]) copyin(db2[0:b2_size])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams
            #pragma omp distribute parallel for collapse( 2 )
#endif
            for( unsigned int i=isBoundary1min; i<n1d-isBoundary1max ; i++ ) {
#ifdef _GPU
                #pragma acc loop worker vector
#endif
                for( unsigned int k=isBoundary2min; k<n2p-isBoundary2max ; k++ ) {
                    B2[ i*nyz_dp + iB1*nz_p + k ]
                        =  Alpha_ *  E1[ i*nyz_pp + p0*nz_p + k ]
                        +  Beta_  *( B2[ i*nyz_dp + p1*nz_p + k ]-B_ext2[ i*n2p + k ] )
                        +  Gamma_ * db2[ i*n2p + k ]
                        +  Zeta_  *( B0[ i*nyz_pd + p0*nz_d + k+1 ]-B_ext0[ i*n2d + (k+1) ] )
                        +  Eta_   *( B0[ i*nyz_pd + p0*nz_d + k   ]-B_ext0[ i*n2d +  k    ] )
                        +  B_ext2[ i*n2p + k ];
                }
            }
        } else {
#ifdef _GPU
            #pragma acc parallel present(E1[0:sizeofE1],B0[0:sizeofB0],B2[0:sizeofB2],B_ext2[0:B_ext_size2],B_ext0[0:B_ext_size0]) copyin(db2[0:b2_size])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams
            #pragma omp distribute parallel for collapse( 2 )
#endif
            for( unsigned int i=isBoundary1min; i<n1d-isBoundary1max ; i++ ) {
#ifdef _GPU
                #pragma acc loop worker vector
#endif
                for( unsigned int j=isBoundary2min; j<n2p-isBoundary2max ; j++ ) {
                    B2[ i*nyz_pd + j*nz_d + iB1 ]
                        = -Alpha_ *  E1[ i*nyz_pp + j*nz_p + p0 ]
                        +  Beta_  *( B2[ i*nyz_pd + j*nz_d  + p1 ]-B_ext2[ i*n2p + j ] )
                        +  Gamma_ * db2[ i*n2p + j ]
                        +  Zeta_  *( B0[ i*nyz_dp + (j+1)*nz_p + p0 ]-B_ext0[ i*n2d + (j+1) ] )
                        +  Eta_   *( B0[ i*nyz_dp +  j   *nz_p + p0 ]-B_ext0[ i*n2d +  j    ] )
                        +  B_ext2[ i*n2p + j ];
                }
            }
        }

        if( !vecLaser.empty() ) {
            smilei::tools::gpu::HostDeviceMemoryManagment::DeviceFree( db1, b1_size );
            smilei::tools::gpu::HostDeviceMemoryManagment::DeviceFree( db2, b2_size );
        }
    }
}
