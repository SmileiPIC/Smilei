#include "ElectroMagnBC3D_refl.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field3D.h"
#include "Tools.h"

using namespace std;


ElectroMagnBC3D_refl::ElectroMagnBC3D_refl( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC3D( params, patch, i_boundary )
{
    // oversize
    if (!params.multiple_decomposition) {
        oversize_x = params.oversize[0];
        oversize_y = params.oversize[1];
        oversize_z = params.oversize[2];
    }
    else {
        oversize_x = params.region_oversize[0];
        oversize_y = params.region_oversize[1];
        oversize_z = params.region_oversize[2];
    }
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_refl::apply( ElectroMagn *EMfields, double, Patch *patch )
{
    const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };
    double *const __restrict__ Bx3D = B[0]->data_;
    double *const __restrict__ By3D = B[1]->data_;
    double *const __restrict__ Bz3D = B[2]->data_;
    const unsigned int nzp   = n_p[2];
    const unsigned int nzd   = n_d[2];
    const unsigned int nxp   = n_p[0];
    const unsigned int nxd   = n_d[0];
    const unsigned int nyp   = n_p[1];
    const unsigned int nyd   = n_d[1];

    const unsigned int nyz_pd = n_p[1] * n_d[2];
    const unsigned int nyz_dp = n_d[1] * n_p[2];
    const unsigned int nyz_dd = n_d[1] * n_d[2];
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    const int sizeofB0 = B[0]->number_of_points_;
    const int sizeofB1 = B[1]->number_of_points_;
    const int sizeofB2 = B[2]->number_of_points_;
#endif
    if( i_boundary_ == 0 && patch->isXmin() ) {
    
        // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
        
        // Static cast of the fields
        //Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        //Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        // for By^(d,p,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(By3D[0:sizeofB1])
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=oversize_x; i>0; i-- ) {
            for( unsigned int j=0 ; j<nyp ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=0 ; k<nzd ; k++ ) {
                    By3D[(i-1)*nyz_pd + j*nzd + k] = By3D[i*nyz_pd + j*nzd + k];
                    //( *By3D )( i-1, j, k ) = ( *By3D )( i, j, k );
                }
            }
        }
        
        // for Bz^(d,d,p)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(Bz3D[0:sizeofB2])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=oversize_x; i>0; i-- ) {
            for( unsigned int j=0 ; j<nyd ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=0 ; k<nzp ; k++ ) {
                    Bz3D[(i-1)*nyz_dp + j*nzp + k] = Bz3D[i*nyz_dp + j*nzp + k];
                    //( *Bz3D )( i-1, j, k ) = ( *Bz3D )( i, j, k );
                }
            }
        }
        
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
    
        // Static cast of the fields
        //Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        //Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        // for By^(d,p,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(By3D[0:sizeofB1])
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=nxd-oversize_x; i<nxd; i++ ) {
            for( unsigned int j=0 ; j<nyp ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=0 ; k<nzd ; k++ ) {
                    By3D[i*nyz_pd + j*nzd + k] = By3D[(i-1)*nyz_pd + j*nzd + k-1];
                    //( *By3D )( i, j, k ) = ( *By3D )( i-1, j, k );
                }
            }
        }
        
        // for Bz^(d,d,p)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(Bz3D[0:sizeofB2])
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=nxd-oversize_x; i<nxd; i++ ) {
            for( unsigned int j=0 ; j<nyd ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=0 ; k<nzp ; k++ ) {
                    Bz3D[i*nyz_dp + j*nzp + k] = Bz3D[(i-1)*nyz_dp + j*nzp + k];
                    //( *Bz3D )( i, j, k ) = ( *Bz3D )( i-1, j, k );
                }
            }
        }
        
    } else if( i_boundary_ == 2 && patch->isYmin() ) {
    
        // Static cast of the fields
        //Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        //Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(Bx3D[0:sizeofB0])
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxp; i++ ) {
            for( unsigned int j=oversize_y ; j>0 ; j-- ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=0; k<nzd ; k++ ) {
                    Bx3D[i*nyz_dd + (j-1)*nzd + k] = Bx3D[i*nyz_pd + j*nzd + k-1]; 
                    //( *Bx3D )( i, j-1, k ) = ( *Bx3D )( i, j, k );
                }
            }
        }
        
        // for Bz^(d,d,p)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(Bz3D[0:sizeofB2])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxd; i++ ) {
            for( unsigned int j=oversize_y ; j>0 ; j-- ) {
                for( unsigned int k=0; k<nzp ; k++ ) {
                    Bz3D[i*nyz_dp + (j-1)*nzp + k] = Bz3D[i*nyz_dp + j*nzp + k]; 
		    //( *Bz3D )( i, j-1, k ) = ( *Bz3D )( i, j, k );
                }
            }
        }
        
    } else if( i_boundary_ == 3 && patch->isYmax() ) {
    
        // Static cast of the fields
        //Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        //Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(Bx3D[0:sizeofB0])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxp; i++ ) {
            for( unsigned int j=nyd-oversize_y; j<nyd ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=0; k<nzd ; k++ ) {
                    Bx3D[i*nyz_dd + j*nzd + k] = Bx3D[i*nyz_dd + (j-1)*nzd + k-1];
                    //( *Bx3D )( i, j, k ) = ( *Bx3D )( i, j-1, k );
                }
            }
        }
        
        // for Bz^(d,d,p)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(Bz3D[0:sizeofB2])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxd; i++ ) {
            for( unsigned int j=nyd-oversize_y; j<nyd ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=0; k<nzp ; k++ ) {
                    Bz3D[i*nyz_dp + j*nzp + k] = Bz3D[i*nyz_dp + (j-1)*nzp + k];
                    //( *Bz3D )( i, j, k ) = ( *Bz3D )( i, j-1, k );
                }
            }
        }
        
    } else if( i_boundary_==4 && patch->isZmin() ) {
    
        // Static cast of the fields
        //Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        //Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(Bx3D[0:sizeofB0])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxp; i++ ) {
            for( unsigned int j=0 ; j<nyd ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=oversize_z ; k>0 ; k-- ) {
                    //( *Bx3D )( i, j, k-1 ) = ( *Bx3D )( i, j, k );
                    Bx3D[i*nyz_dd + j*nzd + k-1] = Bx3D[i*nyz_dd + j*nzd + k];
                }
            }
        }
        
        // for By^(d,p,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(By3D[0:sizeofB1])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxd; i++ ) {
            for( unsigned int j=0 ; j<nyp ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=oversize_z ; k>0 ; k-- ) {
                    //( *By3D )( i, j, k-1 ) = ( *By3D )( i, j, k );
                    By3D[i*nyz_pd + j*nzd + k-1] = By3D[i*nyz_pd + j*nzd + k];
                }
            }
        }
        
    } else if( i_boundary_==5 && patch->isZmax() ) {
    
        // Static cast of the fields
        //Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        //Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(Bx3D[0:sizeofB0])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxp; i++ ) {
            for( unsigned int j=0 ; j<nyd ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=nzd-oversize_z; k<nzd ; k++ ) {
                    //( *Bx3D )( i, j, k ) = ( *Bx3D )( i, j, k-1 );
                    Bx3D[i*nyz_dd + j*nzd + k] = Bx3D[i*nyz_dd + j*nzd + k-1];
                }
            }
        }
        
        // for By^(d,p,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(By3D[0:sizeofB1])
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxd; i++ ) {
            for( unsigned int j=0 ; j<nyp ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                #pragma acc loop worker vector
#endif
                for( unsigned int k=nzd-oversize_z; k<nzd ; k++ ) {
                    //( *By3D )( i, j, k ) = ( *By3D )( i, j, k-1 );
                    By3D[i*nyz_pd + j*nzd + k] = By3D[i*nyz_pd + j*nzd + k-1];
                }
            }
        }
        
    }
}

