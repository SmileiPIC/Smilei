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


ElectroMagnBC3D_refl::ElectroMagnBC3D_refl( Params &params, Patch* patch, unsigned int _min_max )
: ElectroMagnBC( params, patch, _min_max )
{
    // oversize
    oversize_x = params.oversize[0];
    oversize_y = params.oversize[0];
    oversize_z = params.oversize[0];
    
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params.n_space[0]*params.global_factor[0]+1+2*params.oversize[0];
    nx_d = nx_p+1;
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params.n_space[1]*params.global_factor[1]+1+2*params.oversize[1];
    ny_d = ny_p+1;
    // number of nodes of the primal and dual grid in the z-direction
    nz_p = params.n_space[2]*params.global_factor[2]+1+2*params.oversize[2];
    nz_d = nz_p+1;
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_refl::apply(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if (min_max == 0 && patch->isXmin() ) {
        
        // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
        
        // Static cast of the fields
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        // FORCE CONSTANT MAGNETIC FIELDS
        // for By^(d,p,d)
        for (unsigned int i=oversize_x; i>0; i--) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                for (unsigned int k=0 ; k<nz_d ; k++) {
                    (*By3D)(i-1,j,k) = (*By3D)(i,j,k);
                }
            }
        }
        
        // for Bz^(d,d,p)
        for (unsigned int i=oversize_x; i>0; i--) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                for (unsigned int k=0 ; k<nz_p ; k++) {
                    (*Bz3D)(i-1,j,k) = (*Bz3D)(i,j,k);
                }
            }
        }
        
    }
    else if (min_max == 1 && patch->isXmax() ) {
        
        // Static cast of the fields
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        // FORCE CONSTANT MAGNETIC FIELDS
        // for By^(d,p,d)
        for (unsigned int i=nx_d-oversize_x; i<nx_d; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                for (unsigned int k=0 ; k<nz_d ; k++) {
                    (*By3D)(i,j,k) = (*By3D)(i-1,j,k);
                }
            }
        }
        
        // for Bz^(d,d,p)
        for (unsigned int i=nx_d-oversize_x; i<nx_d; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                for (unsigned int k=0 ; k<nz_p ; k++) {
                    (*Bz3D)(i,j,k) = (*Bz3D)(i-1,j,k);
                }
            }
        }
        
    }
    else if (min_max == 2 && patch->isYmin() ) {
        
        // Static cast of the fields
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
        for (unsigned int i=0; i<nx_p; i++) {
            for (unsigned int j=oversize_y ; j>0 ; j--) {
                for (unsigned int k=0; k<nz_d ; k++) {
                    (*Bx3D)(i,j-1,k) = (*Bx3D)(i,j,k);
                }
            }
        }
        
        // for Bz^(d,d,p)
        for (unsigned int i=0; i<nx_d; i++) {
            for (unsigned int j=oversize_y ; j>0 ; j--) {
                for (unsigned int k=0; k<nz_p ; k++) {
                    (*Bz3D)(i,j-1,k) = (*Bz3D)(i,j,k);
                }
            }
        }
        
    }
    else if (min_max == 3 && patch->isYmax() ) {
        
        // Static cast of the fields
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
        for (unsigned int i=0; i<nx_p; i++) {
            for (unsigned int j=ny_d-oversize_y; j<ny_d ; j++) {
                for (unsigned int k=0; k<nz_d ; k++) {
                    (*Bx3D)(i,j,k) = (*Bx3D)(i,j-1,k);
                }
            }
        }
        
        // for Bz^(d,d,p)
        for (unsigned int i=0; i<nx_d; i++) {
            for (unsigned int j=ny_d-oversize_y; j<ny_d ; j++) {
                for (unsigned int k=0; k<nz_p ; k++) {
                    (*Bz3D)(i,j,k) = (*Bz3D)(i,j-1,k);
                }
            }
        }
        
    }
    else if (min_max==4 && patch->isZmin() ) {
        
        // Static cast of the fields
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
        for (unsigned int i=0; i<nx_p; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                for (unsigned int k=oversize_z ; k>0 ; k--) {
                    (*Bx3D)(i,j,k-1) = (*Bx3D)(i,j,k);
                }
            }
        }
        
        // for By^(d,p,d)
        for (unsigned int i=0; i<nx_d; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                for (unsigned int k=oversize_z ; k>0 ; k--) {
                    (*By3D)(i,j,k-1) = (*By3D)(i,j,k);
                }
            }
        }
        
    }
    else if (min_max==5 && patch->isZmax() ) {
        
        // Static cast of the fields
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
        for (unsigned int i=0; i<nx_p; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                for (unsigned int k=nz_d-oversize_z; k<nz_d ; k++) {
                    (*Bx3D)(i,j,k) = (*Bx3D)(i,j,k-1);
                }
            }
        }
        
        // for By^(d,p,d)
        for (unsigned int i=0; i<nx_d; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                for (unsigned int k=nz_d-oversize_z; k<nz_d ; k++) {
                    (*By3D)(i,j,k) = (*By3D)(i,j,k-1);
                }
            }
        }
        
    }
}

