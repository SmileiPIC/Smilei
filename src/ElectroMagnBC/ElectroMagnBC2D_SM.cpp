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
#include "gpu.h"

//using namespace std;

ElectroMagnBC2D_SM::ElectroMagnBC2D_SM( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC2D( params, patch, i_boundary )
{
    
    // Calculate axes indices
    axis0_ = i_boundary_ / 2; // axis normal to border
    axis1_ = axis0_ == 0 ? 1 : 0; // other axis (tangent to border)
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
        iB_[2     ] = n_d[axis0_] - 1;
    }
    
    // Buffers to save B field
    B_val.resize( 3 );
    if( patch->isBoundary( i_boundary_ ) ) {
        B_val[axis0_].resize( n_d[axis1_], 0. ); // dual in the first direction
        B_val[axis1_].resize( n_p[axis1_], 0. ); // primal in the other direction
        B_val[2     ].resize( n_d[axis1_], 0. ); // dual in the other direction

        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( B_val[0].data(), B_val[0].size() );
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( B_val[1].data(), B_val[1].size() );
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( B_val[2].data(), B_val[2].size() );
    }
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    std::vector<double> K = params.EM_BCs_k[i_boundary_];
    double Knorm = sqrt( K[0]*K[0] + K[1]*K[1] ) ;
    double omega = 1.;
    double k0 = omega*K[axis0_] / Knorm;
    double k1 = omega*K[axis1_] / Knorm;
    
    double factor = 1.0 / ( k0 - sign_ * dt_ov_d[axis0_] );
    Alpha_   = 2.0 * factor;
    Beta_    = - ( k0 + sign_ * dt_ov_d[axis0_] ) * factor;
    Gamma_   = 4.0 * k0 * factor;
    Delta_   = - ( k1 + dt_ov_d[axis1_] ) * factor;
    Epsilon_ = - ( k1 - dt_ov_d[axis1_] ) * factor;
}

ElectroMagnBC2D_SM::~ElectroMagnBC2D_SM()
{
    for( auto B: B_val ){
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( B.data(), B.size() );
        //delete[] B;
    }
}

void ElectroMagnBC2D_SM::save_fields( Field *my_field, Patch *patch )
{
    Field2D *field2D=static_cast<Field2D *>( my_field );
    
    if( patch->isBoundary( i_boundary_ ) ) {
        
        unsigned int xyz = 0;
        if( field2D->name=="Bx" ) {
            xyz = 0;
        } else if( field2D->name=="By" ) {
            xyz = 1;
        } else if( field2D->name=="Bz" ) {
            xyz = 2;
        } else {
            return;
        }
        
        if( axis0_ == 0 ) {
            for( unsigned int j=0; j<B_val[xyz].size(); j++ ) {
                B_val[xyz][j]=( *field2D )( iB_[xyz], j );
            }
        } else {
            for( unsigned int i=0; i<B_val[xyz].size(); i++ ) {
                B_val[xyz][i]=( *field2D )( i, iB_[xyz] );
            }
        }
        
    }
}

void ElectroMagnBC2D_SM::disableExternalFields()
{
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( B_val[0].data(), B_val[0].size() );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( B_val[1].data(), B_val[1].size() );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( B_val[2].data(), B_val[2].size() );
    
    B_val[0].resize( 0 );
    B_val[1].resize( 0 );
    B_val[2].resize( 0 );
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_SM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    if( patch->isBoundary( i_boundary_ ) ) {
        
        const Field  *E[3]{ EMfields->Ex_, EMfields->Ey_, EMfields->Ez_ };
        const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };

        const double *const __restrict__ E0 = E[0]->data_;
        const double *const __restrict__ E1 = E[1]->data_;
        const double *const __restrict__ E2 = E[2]->data_;
        double *const __restrict__ B0       = B[0]->data_;
        double *const __restrict__ B1       = B[1]->data_;
        double *const __restrict__ B2       = B[2]->data_;

        const double *const __restrict__ B_ext0 = B_val[0].data();
        const double *const __restrict__ B_ext1 = B_val[1].data();
        const double *const __restrict__ B_ext2 = B_val[2].data();

#ifdef SMILEI_ACCELERATOR_GPU_OACC
        const int sizeofE0 = E[0]->number_of_points_;
        const int sizeofE1 = E[1]->number_of_points_;
        const int sizeofE2 = E[2]->number_of_points_;
        const int sizeofB0 = B[0]->number_of_points_;
        const int sizeofB1 = B[1]->number_of_points_;
        const int sizeofB2 = B[2]->number_of_points_;

        const int B_ext_size0 = B_val[0].size();
        const int B_ext_size1 = B_val[1].size();
        const int B_ext_size2 = B_val[2].size();
#endif

        const int isBoundary1min = patch->isBoundary( axis1_, 0 );
        const int isBoundary1max = patch->isBoundary( axis1_, 1 );

        // Lasers polarized along axis 1
        std::vector<double> b1( n_p[axis1_], 0. );
        double *const __restrict__ db1 = b1.data();
        const unsigned int n1p    = n_p[axis1_];
        const unsigned int n1d    = n_d[axis1_];

        const unsigned int nyp   = n_p[1];
        const unsigned int nyd   = n_d[1]; 
        const unsigned int iB0    = iB_[0];
        const unsigned int p0     = iB_[0] - sign_;
        const unsigned int p1     = iB_[1] - sign_;
        const unsigned int iB1    = iB_[1];
        const unsigned int iB2    = iB_[2];
        const unsigned int p2     = iB_[2] - sign_;

        const int b1_size = n1p ;
        const int b2_size = n1d ;
        std::vector<double> pos( 1 );

        if( ! vecLaser.empty() ) {
            for( unsigned int j=isBoundary1min ; j<n1p-isBoundary1max ; j++ ) {
                pos[0] = patch->getDomainLocalMin( axis1_ ) + ( ( int )j - ( int )EMfields->oversize[axis1_] )*d[axis1_];
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    db1[j] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, 0 );
                }
            }
        }
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocateAndCopyHostToDevice( db1, b1_size );

        if( axis0_ == 0 ) { // for By^(d,p)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(E2[0:sizeofE2],B0[0:sizeofB0],B1[0:sizeofB1],B_ext1[0:B_ext_size1],B_ext0[0:B_ext_size0],db1[0:b1_size])
            #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#endif
            for( unsigned int j=isBoundary1min; j<n1p-isBoundary1max ; j++ ) {
                B1[ iB1*nyp + j ]
                    = Alpha_  * E2[ iB0*nyp + j ]                                                       
                    + Beta_   * ( B1[p1*nyp + j] - B_ext1[j] )                                          
                    + Gamma_  * db1[j]
                    + Delta_  * ( B0[iB0*nyd + j+1 ] - B_ext0[j+1] )                                    
                    + Epsilon_* ( B0[iB0*nyd + j   ] - B_ext0[j  ] )
                    + B_ext1[j];
            }
        } else { // for Bx^(p,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(E2[0:sizeofE2],B0[0:sizeofB0],B1[0:sizeofB1],B_ext1[0:B_ext_size1],B_ext0[0:B_ext_size0],db1[0:b1_size])
            #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#endif
            for( unsigned int j=isBoundary1min; j<n1p-isBoundary1max ; j++ ) {
                B0[ iB0 + nyd * j ]
                    = -Alpha_  *   E2[iB1+nyp * j]
                    + Beta_   * ( B0[p0+nyd * j ] - B_ext0[j] )
                    + Gamma_  * db1[j]
                    + Delta_  * ( B1[iB1+nyp * (j+1) ] - B_ext1[j+1] )
                    + Epsilon_* ( B1[iB1+nyp * j   ] - B_ext1[j  ] )
                    + B_ext0[j];
            }
        }
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( db1, b1_size );

        // Lasers polarized along axis 2
        std::vector<double> b2( n_d[axis1_], 0. );
        double *const __restrict__ db2 = b2.data();

        if( ! vecLaser.empty() ) {
            for( unsigned int j=isBoundary1min; j<n1d-isBoundary1max ; j++ ) {
                pos[0] = patch->getDomainLocalMin( axis1_ ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[axis1_] )*d[axis1_];
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    db2[j] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, 0 );
                }
            }
        }
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocateAndCopyHostToDevice( db2, b2_size );

        // for Bz^(d,d)
        if( axis0_ == 0 ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(E1[0:sizeofE1],B2[0:sizeofB2],B_ext2[0:B_ext_size2],db2[0:b2_size])
            #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#endif
            for( unsigned int j=isBoundary1min; j<n1d-isBoundary1max ; j++ ) {
                B2[ iB2*nyd + j ] = -Alpha_ * E1[iB0*nyd + j] + Beta_ * (B2[ p2*nyd + j] - B_ext2[j]) 
                    + Gamma_ * db2[j] + B_ext2[j];

            }
        } else {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(E0[0:sizeofE0],B2[0:sizeofB2],B_ext2[0:B_ext_size2],db2[0:b2_size])
            #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#endif
            for( unsigned int j=isBoundary1min; j<n1d-isBoundary1max ; j++ ) {
                B2[ j*nyd + iB2 ] = Alpha_ * E0[j*nyp + iB1] + Beta_ * (B2[j*nyd + p2] - B_ext2[j])
                    + Gamma_ * db2[j] + B_ext2[j];
            }
        }
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( db2, b2_size );
    }
}
