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

ElectroMagnBC2D_refl::ElectroMagnBC2D_refl( Params &params, LaserParams &laser_params )
: ElectroMagnBC( params, laser_params )
{
    // oversize
    oversize_ = params.oversize[0];
    
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params.n_space[0]+1+2*params.oversize[0];
    nx_d = params.n_space[0]+2+2*params.oversize[0];
 
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params.n_space[1]+1+2*params.oversize[1];
    ny_d = params.n_space[1]+2+2*params.oversize[1];
    
}

ElectroMagnBC2D_refl::~ElectroMagnBC2D_refl()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_refl::apply_xmin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isWestern() ) {
        
        // Static cast of the fields
        Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
        Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
        Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
        Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);

        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d)
        for (unsigned int i=oversize_; i>0; i--) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Bx2D)(i-1,j) = (*Bx2D)(i,j);
            }//j
        }//i
        
        // for By^(d,p)
        for (unsigned int i=oversize_; i>0; i--) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*By2D)(i-1,j) = (*By2D)(i,j);
            }//j
        }//i
        
        // for Bz^(d,d)
        for (unsigned int i=oversize_; i>0; i--) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Bz2D)(i-1,j) = (*Bz2D)(i,j);
            }
        }
        
        // FORCE ZERO ELECTRIC FIELDS
        
        // for Ex^(d,p)
        for (unsigned int i=0; i<oversize_+1; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Ex2D)(i,j) = 0.0;
            }//j
        }//i
        
        // for Ey^(p,d)
        for (unsigned int i=0; i<oversize_; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Ey2D)(i,j) = 0.0;
            }//j
        }//i
        
        // for Ez^(p,p)
        for (unsigned int i=0; i<oversize_; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Ez2D)(i,j) = 0.0;
            }
        }
        
    }//if Western
    
}
// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_refl::apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isEastern() ) {
        
        // Static cast of the fields
        Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
        Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
        Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
        Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d)
        for (unsigned int i=nx_p-oversize_; i<nx_p; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Bx2D)(i,j) = (*Bx2D)(i-1,j);
            }//j
        }//i
        
        // for By^(d,p)
        for (unsigned int i=nx_d-oversize_; i<nx_d; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*By2D)(i,j) = (*By2D)(i-1,j);
            }//j
        }//i
        
        // for Bz^(d,d)
        for (unsigned int i=nx_d-oversize_; i<nx_d; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Bz2D)(i,j) = (*Bz2D)(i-1,j);
            }
        }
        
        // FORCE ZERO ELECTRIC FIELDS
        
        // for Ex^(d,p)
        for (unsigned int i=nx_d-oversize_; i<nx_d; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Ex2D)(i,j) = 0.0;
            }//j
        }//i
        
        // for Ey^(p,d)
        for (unsigned int i=nx_p-oversize_; i<nx_p; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Ey2D)(i,j) = 0.0;
            }//j
        }//i
        
        // for Ez^(p,p)
        for (unsigned int i=nx_p-oversize_; i<nx_p; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Ez2D)(i,j) = 0.0;
            }
        }

        
    }//if Eastern
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_refl::apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isSouthern() ) {
        
        // Static cast of the fields
        //Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        //Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        //Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
        //Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
        //Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
        //Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);
        
        ERROR("NOT YET DEFINED");
        
    }//if Southern
    
}
// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_refl::apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isNorthern() ) {
        
        // Static cast of the fields
        //Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        //Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        //Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
        //Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
        //Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
        //Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);
        
        ERROR("NOT YET DEFINED");
        
    }//if Northern
    
}


