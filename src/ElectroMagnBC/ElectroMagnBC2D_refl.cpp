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

 MG to NICO:
 THIS IS UNCLEAR FOR ME, WE WANT THE ALL GHOST CELL TO BEHAVE AS A PERFECT CONDUCTOR
 SO THAT WE HAVE REFLECTION AT THE EXACT SIMULATION BORDER - NOT INSIDE THE GHOST CELL

 Perfect Conducting Boundary Conditions for electromagnetic fields
 dBt/dn = 0 (tangential magnetic components have zero derivative)
 Bn = 0 (normal flux is zero)

 we only need to fix the magnetic field. The electric field is
 calculated in Maxwell-Ampere on the whole grid already
 (knowing B and J on the whole grid).

 */




ElectroMagnBC2D_refl::ElectroMagnBC2D_refl( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC2D( params, patch, i_boundary )
{
    // oversize
    if (!params.uncoupled_grids)
        oversize_ = params.oversize[0];
    else
        oversize_ = params.region_oversize[0];
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_refl::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    if( i_boundary_ == 0 && patch->isXmin() ) {
    
        // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
        // Static cast of the fields
        //Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        //Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        //Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
        //Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
        Field2D *By2D = static_cast<Field2D *>( EMfields->By_ );
        Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        //        // for Bx^(p,d)
        //        for (unsigned int i=oversize_; i>0; i--) {
        //            for (unsigned int j=0 ; j<n_d[1] ; j++) {
        //                (*Bx2D)(i-1,j) = (*Bx2D)(i,j);
        //            }//j
        //        }//i
        
        // for By^(d,p)
        for( unsigned int i=oversize_; i>0; i-- ) {
            for( unsigned int j=0 ; j<n_p[1] ; j++ ) {
                ( *By2D )( i-1, j ) = ( *By2D )( i, j );
            }//j
        }//i
        
        // for Bz^(d,d)
        for( unsigned int i=oversize_; i>0; i-- ) {
            for( unsigned int j=0 ; j<n_d[1] ; j++ ) {
                ( *Bz2D )( i-1, j ) = ( *Bz2D )( i, j );
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
    
        // Static cast of the fields
        //Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        //Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        //Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
        //Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
        Field2D *By2D = static_cast<Field2D *>( EMfields->By_ );
        Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        //        // for Bx^(p,d)
        //        for (unsigned int i=n_p[0]-oversize_; i<n_p[0]; i++) {
        //            for (unsigned int j=0 ; j<n_d[1] ; j++) {
        //                (*Bx2D)(i,j) = (*Bx2D)(i-1,j);
        //            }//j
        //        }//i
        
        // for By^(d,p)
        for( unsigned int i=n_d[0]-oversize_; i<n_d[0]; i++ ) {
            for( unsigned int j=0 ; j<n_p[1] ; j++ ) {
                ( *By2D )( i, j ) = ( *By2D )( i-1, j );
            }//j
        }//i
        
        // for Bz^(d,d)
        for( unsigned int i=n_d[0]-oversize_; i<n_d[0]; i++ ) {
            for( unsigned int j=0 ; j<n_d[1] ; j++ ) {
                ( *Bz2D )( i, j ) = ( *Bz2D )( i-1, j );
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
    
        // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
        // Static cast of the fields
        //Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        //Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        //Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
        Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
        //Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
        Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d)
        for( unsigned int i=0; i<n_p[0]; i++ ) {
            for( unsigned int j=oversize_ ; j>0 ; j-- ) {
                ( *Bx2D )( i, j-1 ) = ( *Bx2D )( i, j );
            }//j
        }//i
        
        //        // for By^(d,p)
        //        for (unsigned int i=0; i<n_d[0]; i++) {
        //            for (unsigned int j=oversize_ ; j>0 ; j--) {
        //                (*By2D)(i,j-1) = (*By2D)(i,j);
        //            }//j
        //        }//i
        
        // for Bz^(d,d)
        for( unsigned int i=0; i<n_d[0]; i++ ) {
            for( unsigned int j=oversize_ ; j>0 ; j-- ) {
                ( *Bz2D )( i, j-1 ) = ( *Bz2D )( i, j );
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
    
        // Static cast of the fields
        //Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
        //Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
        //Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
        Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
        //Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
        Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d)
        for( unsigned int i=0; i<n_p[0]; i++ ) {
            for( unsigned int j=n_d[1]-oversize_; j<n_d[1] ; j++ ) {
                ( *Bx2D )( i, j ) = ( *Bx2D )( i, j-1 );
            }//j
        }//i
        
        //        // for By^(d,p)
        //        for (unsigned int i=0; i<n_d[0]; i++) {
        //            for (unsigned int j=n_p[1]-oversize_ ; j<n_p[1] ; j++) {
        //                (*By2D)(i,j) = (*By2D)(i,j-1);
        //            }//j
        //        }//i
        
        // for Bz^(d,d)
        for( unsigned int i=0; i<n_d[0]; i++ ) {
            for( unsigned int j=n_d[1]-oversize_; j<n_d[1] ; j++ ) {
                ( *Bz2D )( i, j ) = ( *Bz2D )( i, j-1 );
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

