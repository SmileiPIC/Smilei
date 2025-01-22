#include "ElectroMagnBC1D_refl.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Tools.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Reflective boundary conditions
// The ghost cells behave as a Perfect Condutor:
//      - electric fields are put to zero in the ghost cells
//      - magnetic fields are constant in the ghost cells
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagnBC1D_refl::ElectroMagnBC1D_refl( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC1D( params, patch, i_boundary )
{

    // oversize
    oversize_ = params.oversize[0];
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Reflective Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_refl::apply( ElectroMagn *EMfields, double, Patch *patch )
{
    if( i_boundary_ == 0 ) {
        if( patch->isXmin() ) {
            const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };
            double *const __restrict__ By1D = B[1]->data_;
            double *const __restrict__ Bz1D = B[2]->data_;
        
            // Application over the full-ghost cell
            //Field1D *By1D   = static_cast<Field1D *>( EMfields->By_ );
            //Field1D *Bz1D   = static_cast<Field1D *>( EMfields->Bz_ );
            
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        const int sizeofB1 = B[1]->number_of_points_;
        const int sizeofB2 = B[2]->number_of_points_;
#endif
            // force constant magnetic fields in the ghost cells
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(By1D[0:sizeofB1],Bz1D[0:sizeofB2])
            #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#endif
            for( unsigned int i=oversize_; i>0; i-- ) {
                //( *By1D )( i-1 ) = ( *By1D )( i );
                //( *Bz1D )( i-1 ) = ( *Bz1D )( i );
                By1D[i-1] = By1D[i];
                Bz1D[i-1] = Bz1D[i];
            }
            
            //        // force 0 electric fields in the ghost cells
            //        for (unsigned int i=0; i<oversize_+1; i++) {
            //            (*Ex1D)(i) = 0.0;
            //        }
            //        for (unsigned int i=0; i<oversize_; i++) {
            //            (*Ey1D)(i) = 0.0;
            //            (*Ez1D)(i) = 0.0;
            //        }
            
            /* DEFINITION BY NICO
            
             // The other fields are already defined everywhere, so no need for boundary conditions.
             Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
             Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
            
             // normal derivative of tangential B is zero.
             // By and Bz just outside equal By and Bz just inside.
             (*By1D)(0) = (*By1D)(1);
             (*Bz1D)(0) = (*Bz1D)(1);
             */
            
        }//if Xmin
    } else {
        if( patch->isXmax() ) {
            const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };
            double *const __restrict__ By1D = B[1]->data_;
            double *const __restrict__ Bz1D = B[2]->data_;
        
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            const int sizeofB1 = B[1]->number_of_points_;
            const int sizeofB2 = B[2]->number_of_points_;
#endif
            const unsigned int nxd   = n_d[0]; 
            // application of Bcs over the full ghost cells
            //Field1D *By1D   = static_cast<Field1D *>( EMfields->By_ );
            //Field1D *Bz1D   = static_cast<Field1D *>( EMfields->Bz_ );
            
            // force constant magnetic fields in the ghost cells
            //        for (unsigned int i=n_p[0]-oversize_; i<n_p[0]; i++)
            //            (*Bx1D)(i) = (*Bx1D)(i-1);
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc parallel present(By1D[0:sizeofB1],Bz1D[0:sizeofB2])
            #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#endif
            for( unsigned int i=nxd-oversize_; i<nxd; i++ ) {
                //( *By1D )( i ) = ( *By1D )( i-1 );
                //( *Bz1D )( i ) = ( *Bz1D )( i-1 );
                By1D[i] = By1D[i-1];
                Bz1D[i] = Bz1D[i-1];
            }
            
            //        // force 0 electric fields in the ghost cells
            //        for (unsigned int i=n_d[0]-oversize_; i<n_d[0]; i++)
            //            (*Ex1D)(i) = 0.0;
            //        for (unsigned int i=n_p[0]-oversize_; i<n_p[0]; i++) {
            //            (*Ey1D)(i) = 0.0;
            //            (*Ez1D)(i) = 0.0;
            //        }
            
            /* DEFINITION BY NICO
             Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
             Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
            
             // normal derivative of tangential B is zero.
             // By and Bz just outside equal By and Bz just inside.
             (*By1D)(n_d[0]-1) = (*By1D)(n_d[0]-2);
             (*Bz1D)(n_d[0]-1) = (*Bz1D)(n_d[0]-2);
             */
            
        }//if Xmax
        
    }
}

