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


/*

 Perfect Conducting Boundary Conditions for electromagnetic fields
 dBt/dn = 0 (tangential magnetic components have zero derivative)
 Bn = 0 (normal flux is zero)

 we only need to fix the magnetic field. The electric field is
 calculated in Maxwell-Ampere on the whole grid already
 (knowing B and J on the whole grid).

*/




ElectroMagnBC2D_refl::ElectroMagnBC2D_refl( Params &params, Patch* patch )
  : ElectroMagnBC( params, patch )
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
        // The other fields are already defined, so no need for boundary conditions.
        Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
        Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);

        // perfect conducting wall (no flux through it)
        // normal derivative of tangential B = 0 <--> dBt/dn
        // normal component of B is zero <--> Bn=0
        // By and Bz just outside equal By and Bz just inside.
        // for By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D)(0,j) = (*By2D)(1,j);
        }


        // for Bz^(d,d)
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D)(0,j) = (*Bz2D)(1,j);
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
        // The other fields are already defined, so no need for boundary conditions.
        Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
        Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);

        // perfect conducting wall (no flux through it)
        // normal derivative of tangential B = 0 <--> dBt/dn
        // normal component of B is zero <--> Bn=0
        // By and Bz just outside equal By and Bz just inside.

        // for By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D)(nx_d-1,j) = (*By2D)(nx_d-2,j);
        }//j

        // for Bz^(d,d)
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D)(nx_d-1,j) = (*Bz2D)(nx_d-2,j);
        }//j



    }//if Eastern

}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_refl::apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isSouthern() ) {

        // Static cast of the fields
        // The other fields are already defined, so no need for boundary conditions.
        Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
        Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);

        // perfect conducting wall (no flux through it)
        // normal derivative of tangential B = 0 <--> dBt/dn
        // normal component of B is zero <--> Bn=0
        // By and Bz just outside equal By and Bz just inside.

        // for Bx^(p,d)
        for (unsigned int i=0 ; i<nx_p; i++) {
            (*Bx2D)(i,0) = (*Bx2D)(i,1);
        }


        // for Bz^(d,d)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            (*Bz2D)(i,0) = (*Bz2D)(i,1);
        }


    }//if Southern

}
// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_refl::apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isNorthern() ) {

        // Static cast of the fields
        // The other fields are already defined, so no need for boundary conditions.
        Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
        Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);

        // perfect conducting wall (no flux through it)
        // normal derivative of tangential B = 0 <--> dBt/dn
        // normal component of B is zero <--> Bn=0
        // By and Bz just outside equal By and Bz just inside.

        // for Bx^(p,d)
        for (unsigned int i=0 ; i<nx_p ; i++) {
            (*Bx2D)(i,ny_d-1) = (*Bx2D)(i,ny_d-2);
        }//i

        // for Bz^(d,d)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            (*Bz2D)(i,ny_d-1) = (*Bz2D)(i,ny_d-2);
        }//i



    }//if Northern

}


