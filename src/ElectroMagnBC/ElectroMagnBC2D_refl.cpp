#include "ElectroMagnBC2D_refl.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field2D.h"
#include "Tools.h"

using namespace std;

ElectroMagnBC2D_refl::ElectroMagnBC2D_refl( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC2D( params, patch, i_boundary )
{
    // oversize
    if (!params.multiple_decomposition)
        oversize_ = params.oversize[0];
    else
        oversize_ = params.region_oversize[0];
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_refl::apply( ElectroMagn *EMfields, double, Patch *patch )
{
    if( i_boundary_ == 0 && patch->isXmin() ) {
    
        const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };
        double *const __restrict__ By2D = B[1]->data_;
        double *const __restrict__ Bz2D = B[2]->data_;
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        const int sizeofB1 = B[1]->number_of_points_;
        const int sizeofB2 = B[2]->number_of_points_;
#endif
        const unsigned int nyp   = n_p[1];
        const unsigned int nyd   = n_d[1]; 
        // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
        // Static cast of the fields
        //Field2D *By2D = static_cast<Field2D *>( EMfields->By_ );
        //Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        
        // for By^(d,p)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(By2D[0:sizeofB1])
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=oversize_; i>0; i-- ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop worker vector
#endif
            for( unsigned int j=0 ; j<nyp ; j++ ) {
                //( *By2D )( i-1, j ) = ( *By2D )( i, j );
                By2D[(i-1)*nyp + j] = By2D[i*nyp + j];
            }
        }
        
        // for Bz^(d,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(Bz2D[0:sizeofB2],)
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=oversize_; i>0; i-- ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop worker vector
#endif
            for( unsigned int j=0 ; j<nyd ; j++ ) {
                //( *Bz2D )( i-1, j ) = ( *Bz2D )( i, j );
                Bz2D[(i-1)*nyd + j] = Bz2D[i*nyd + j];
            }
        }
        
        //        // FORCE ZERO ELECTRIC FIELDS
        //
        //        // for Ex^(d,p)
        //        for (unsigned int i=0; i<n_d[0]; i++) {
        //            for (unsigned int j=0 ; j<oversize_ ; j++) {
        //                (*Ex2D)(i,j) = 0.0;
        //            }//j
        //        }//i
        //
        //        // for Ey^(p,d)
        //        for (unsigned int i=0; i<n_p[0]; i++) {
        //            for (unsigned int j=0 ; j<oversize_+1 ; j++) {
        //                (*Ey2D)(i,j) = 0.0;
        //            }//j
        //        }//i
        //
        //        // for Ez^(p,p)
        //        for (unsigned int i=0; i<n_p[0]; i++) {
        //            for (unsigned int j=0 ; j<oversize_ ; j++) {
        //                (*Ez2D)(i,j) = 0.0;
        //            }
        //        }
        
        /* DEFINITION BY NICO
         // Static cast of the fields
         // The other fields are already defined, so no need for boundary conditions.
         Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
         Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);
        
         // perfect conducting wall (no flux through it)
         // normal derivative of tangential B = 0 <--> dBt/dn
         // normal component of B is zero <--> Bn=0
         // By and Bz just outside equal By and Bz just inside.
         // for By^(d,p)
         for (unsigned int j=0 ; j<n_p[1] ; j++) {
         (*By2D)(0,j) = (*By2D)(1,j);
         }
        
        
         // for Bz^(d,d)
         for (unsigned int j=0 ; j<n_d[1] ; j++) {
         (*Bz2D)(0,j) = (*Bz2D)(1,j);
         }
         */
        
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
    
        const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };
        double *const __restrict__ By2D = B[1]->data_;
        double *const __restrict__ Bz2D = B[2]->data_;
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        const int sizeofB1 = B[1]->number_of_points_;
        const int sizeofB2 = B[2]->number_of_points_;
#endif
        const unsigned int nxp   = n_p[0];
        const unsigned int nxd   = n_d[0]; 
        const unsigned int nyp   = n_p[1];
        const unsigned int nyd   = n_d[1]; 
        // Static cast of the fields
        //Field2D *By2D = static_cast<Field2D *>( EMfields->By_ );
        //Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for By^(d,p)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(By2D[0:sizeofB1],)
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for
#endif
        for( unsigned int i = nxd - oversize_; i < nxd; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop worker vector
#endif
            for( unsigned int j = 0 ; j < nyp ; j++ ) {
                By2D[i*nyp + j] = By2D[ (i-1)*nyp + j];
                //( *By2D )( i, j ) = ( *By2D )( i-1, j );
            }
        }
        
        // for Bz^(d,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(Bz2D[0:sizeofB2],)
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i = nxd - oversize_; i < nxd; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop worker vector
#endif
            for( unsigned int j = 0 ; j < nyd ; j++ ) {
                Bz2D[i*nyd + j] = Bz2D[(i-1)*nyd + j];
                //( *Bz2D )( i, j ) = ( *Bz2D )( i-1, j );
            }
        }
        
        //        // FORCE ZERO ELECTRIC FIELDS
        //
        //        // for Ex^(d,p)
        //        for (unsigned int i=n_d[0]-oversize_; i<n_d[0]; i++) {
        //            for (unsigned int j=0 ; j<n_p[1] ; j++) {
        //                (*Ex2D)(i,j) = 0.0;
        //            }//j
        //        }//i
        //
        //        // for Ey^(p,d)
        //        for (unsigned int i=n_p[0]-oversize_; i<n_p[0]; i++) {
        //            for (unsigned int j=0 ; j<n_d[1] ; j++) {
        //                (*Ey2D)(i,j) = 0.0;
        //            }//j
        //        }//i
        //
        //        // for Ez^(p,p)
        //        for (unsigned int i=n_p[0]-oversize_; i<n_p[0]; i++) {
        //            for (unsigned int j=0 ; j<n_p[1] ; j++) {
        //                (*Ez2D)(i,j) = 0.0;
        //            }
        //        }
        
        /* DEFINITION BY NICO
         // Static cast of the fields
         // The other fields are already defined, so no need for boundary conditions.
         Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
         Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);
        
         // perfect conducting wall (no flux through it)
         // normal derivative of tangential B = 0 <--> dBt/dn
         // normal component of B is zero <--> Bn=0
         // By and Bz just outside equal By and Bz just inside.
        
         // for By^(d,p)
         for (unsigned int j=0 ; j<n_p[1] ; j++) {
         (*By2D)(n_d[0]-1,j) = (*By2D)(n_d[0]-2,j);
         }//j
        
         // for Bz^(d,d)
         for (unsigned int j=0 ; j<n_d[1] ; j++) {
         (*Bz2D)(n_d[0]-1,j) = (*Bz2D)(n_d[0]-2,j);
         }//j
         */
        
        
        
    } else if( i_boundary_ == 2 && patch->isYmin() ) {
    
        const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };
        double *const __restrict__ Bx2D = B[0]->data_;
        double *const __restrict__ Bz2D = B[2]->data_;
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        const int sizeofB0 = B[0]->number_of_points_;
        const int sizeofB2 = B[2]->number_of_points_;
#endif
        const unsigned int nxp   = n_p[0];
        const unsigned int nxd   = n_d[0]; 
        const unsigned int nyp   = n_p[1];
        const unsigned int nyd   = n_d[1]; 
        // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
        // Static cast of the fields
        //Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
        //Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(Bx2D[0:sizeofB0],)
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxp; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop worker vector
#endif
            for( unsigned int j=oversize_ ; j>0 ; j-- ) {
                //( *Bx2D )( i, j-1 ) = ( *Bx2D )( i, j );
                Bx2D[i*nyd + j-1] = Bx2D[i*nyd + j];
            }
        }
        
        // for Bz^(d,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(Bz2D[0:sizeofB2],)
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxd; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop worker vector
#endif
            for( unsigned int j = oversize_ ; j>0 ; j-- ) {
                //( *Bz2D )( i, j-1 ) = ( *Bz2D )( i, j );
                Bz2D[i*nyd + j-1] = Bz2D[i*nyd + j];
            }
        }
        
        //        // FORCE ZERO ELECTRIC FIELDS
        //
        //        // for Ex^(d,p)
        //        for (unsigned int i=0; i<n_d[0]; i++){
        //            for (unsigned int j=0 ; j<oversize_ ; j++) {
        //                (*Ex2D)(i,j) = 0.0;
        //            }//j
        //        }//i
        //
        //        // for Ey^(p,d)
        //        for (unsigned int i=0; i<n_p[0]; i++) {
        //            for (unsigned int j=0 ; j<oversize_ ; j++) {//j<oversize_+1???
        //                (*Ey2D)(i,j) = 0.0;
        //            }//j
        //        }//i
        //
        //        // for Ez^(p,p)
        //        for (unsigned int i=0; i<n_p[0]; i++) {
        //            for (unsigned int j=0 ; j<oversize_ ; j++) {
        //                (*Ez2D)(i,j) = 0.0;
        //            }
        //        }
        /* NICO BCs
         // Static cast of the fields
         // The other fields are already defined, so no need for boundary conditions.
         Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
         Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);
        
         // perfect conducting wall (no flux through it)
         // normal derivative of tangential B = 0 <--> dBt/dn
         // normal component of B is zero <--> Bn=0
         // By and Bz just outside equal By and Bz just inside.
        
         // for Bx^(p,d)
         for (unsigned int i=0 ; i<n_p[0]; i++) {
         (*Bx2D)(i,0) = (*Bx2D)(i,1);
         }
        
         // for Bz^(d,d)
         for (unsigned int i=0 ; i<n_d[0] ; i++) {
         (*Bz2D)(i,0) = (*Bz2D)(i,1);
         }
         */
        
    } else if( i_boundary_ == 3 && patch->isYmax() ) {
    
        const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };
        double *const __restrict__ Bx2D = B[0]->data_;
        double *const __restrict__ Bz2D = B[2]->data_;
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        const int sizeofB0 = B[0]->number_of_points_;
        const int sizeofB2 = B[2]->number_of_points_;
#endif
        const unsigned int nxp   = n_p[0];
        const unsigned int nxd   = n_d[0]; 
        const unsigned int nyp   = n_p[1];
        const unsigned int nyd   = n_d[1]; 
        // Static cast of the fields
        //Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
        //Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(Bx2D[0:sizeofB0],)
        #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
        #pragma omp teams distribute parallel for
#endif
        for( unsigned int i=0; i<nxp; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop worker vector
#endif
            for( unsigned int j=nyd-oversize_; j<nyd ; j++ ) {
                //( *Bx2D )( i, j ) = ( *Bx2D )( i, j-1 );
                Bx2D[i*nyd + j] = Bx2D[i*nyd + j-1];
            }
        }
        
        // for Bz^(d,d)
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(Bz2D[0:sizeofB2],)
            #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for 
#endif
        for( unsigned int i=0; i<nxd; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop worker vector
#endif
            for( unsigned int j=nyd-oversize_; j<nyd ; j++ ) {
                //( *Bz2D )( i, j ) = ( *Bz2D )( i, j-1 );
                Bz2D[i*nyd + j] = Bz2D[i*nyd + j-1];
            }
        }
        
        //        // FORCE ZERO ELECTRIC FIELDS
        //
        //        // for Ex^(d,p)
        //        for (unsigned int i=0; i<n_d[0]; i++) {
        //            for (unsigned int j=n_p[1]-oversize_ ; j<n_p[1] ; j++) {
        //                (*Ex2D)(i,j) = 0.0;
        //            }//j
        //        }//i
        //
        //        // for Ey^(p,d)
        //        for (unsigned int i=0; i<n_p[0]; i++) {
        //            for (unsigned int j=n_d[1]-oversize_ ; j<n_d[1] ; j++) {
        //                (*Ey2D)(i,j) = 0.0;
        //            }//j
        //        }//i
        //
        //        // for Ez^(p,p)
        //        for (unsigned int i=0; i<n_p[0]; i++) {
        //            for (unsigned int j=n_p[1]-oversize_ ; j<n_p[1] ; j++) {
        //                (*Ez2D)(i,j) = 0.0;
        //            }
        //        }
        
        
        /* DEFINITION BY NICO
         // Static cast of the fields
         // The other fields are already defined, so no need for boundary conditions.
         Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
         Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);
        
         // perfect conducting wall (no flux through it)
         // normal derivative of tangential B = 0 <--> dBt/dn
         // normal component of B is zero <--> Bn=0
         // By and Bz just outside equal By and Bz just inside.
        
         // for Bx^(p,d)
         for (unsigned int i=0 ; i<n_p[0] ; i++) {
         (*Bx2D)(i,n_d[1]-1) = (*Bx2D)(i,n_d[1]-2);
         }//i
        
         // for Bz^(d,d)
         for (unsigned int i=0 ; i<n_d[0] ; i++) {
         (*Bz2D)(i,n_d[1]-1) = (*Bz2D)(i,n_d[1]-2);
         }//i
         */
        
    }
}
