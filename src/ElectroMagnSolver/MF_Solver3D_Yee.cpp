
#include "MF_Solver3D_Yee.h"

#include "ElectroMagn.h"
#include "Field3D.h"

MF_Solver3D_Yee::MF_Solver3D_Yee(Params &params)
: Solver3D(params)
{
}

MF_Solver3D_Yee::~MF_Solver3D_Yee()
{
}

void MF_Solver3D_Yee::operator() ( ElectroMagn* fields )
{
    // Static-cast of the fields
    Field3D* Ex3D = static_cast<Field3D*>(fields->Ex_);
    Field3D* Ey3D = static_cast<Field3D*>(fields->Ey_);
    Field3D* Ez3D = static_cast<Field3D*>(fields->Ez_);
    Field3D* Bx3D = static_cast<Field3D*>(fields->Bx_);
    Field3D* By3D = static_cast<Field3D*>(fields->By_);
    Field3D* Bz3D = static_cast<Field3D*>(fields->Bz_);
    
    // Magnetic field Bx^(p,d,d)
    for (unsigned int i=0 ; i<nx_p;  i++) {
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            for (unsigned int k=1 ; k<nz_d-1 ; k++) {
#ifdef _PATCH3D_TODO
                (*Bx3D)(i,j,k) -= ;
#endif
            }
        }
    }
        
    // Magnetic field By^(d,p,d)
    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            for (unsigned int k=1 ; k<nz_d-1 ; k++) {
#ifdef _PATCH3D_TODO
                (*By3D)(i,j,k) += ;
#endif
            }
        }
    }
        
    // Magnetic field Bz^(d,d,p)
    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            for (unsigned int k=0 ; k<nz_p ; k++) {
#ifdef _PATCH3D_TODO
                (*Bz3D)(i,j) += ;
#endif
            }
        }
    }

}

