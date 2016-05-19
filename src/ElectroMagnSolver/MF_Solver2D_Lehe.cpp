
#include "MF_Solver2D_Lehe.h"

#include "ElectroMagn.h"
#include "Field2D.h"

#include <algorithm>

MF_Solver2D_Lehe::MF_Solver2D_Lehe(Params &params)
    : Solver2D(params)
{
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    delta_x = (1./4.)*(1.-pow( sin(M_PI*dt_ov_dx/2.)/dt_ov_dx,2 ) );
    beta_x =   1./8.;
    beta_y =   pow(dx/dy,2)/8.;
    Beta_x =  1.-2.*beta_x;
    Beta_y =  1.-2.*beta_y-3.*delta_x;
}

MF_Solver2D_Lehe::~MF_Solver2D_Lehe()
{
}

void MF_Solver2D_Lehe::operator() ( ElectroMagn* fields )
{
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(fields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(fields->Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(fields->Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(fields->Bx_);
    Field2D* By2D = static_cast<Field2D*>(fields->By_);
    Field2D* Bz2D = static_cast<Field2D*>(fields->Bz_);



//    for (unsigned int i=0 ; i<nx_p;  i++) {
//    for (unsigned int i=1 ; i<nx_d-1;  i++) {
    for (unsigned int i=1 ; i<nx_d-2;  i++) {
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*Bx2D)(i,j) -= dt_ov_dy * ( Beta_x*((*Ez2D)(i,j) - (*Ez2D)(i,j-1)) +beta_x*( (*Ez2D)(i+1,j)- (*Ez2D)(i+1,j-1) +(*Ez2D)(i-1,j)- (*Ez2D)(i-1,j-1)));
	}
    }
    
    // Magnetic field By^(d,p)
    for (unsigned int i=2 ; i<nx_d-2 ; i++) {
        for (unsigned int j=1 ; j<ny_p-1 ; j++) {
            (*By2D)(i,j) += dt_ov_dx * ( Beta_y*((*Ez2D)(i,j) - (*Ez2D)(i-1,j)) +beta_y*((*Ez2D)(i,j+1)- (*Ez2D)(i-1,j+1) +(*Ez2D)(i,j-1)- (*Ez2D)(i-1,j-1) ) +delta_x*( (*Ez2D)(i+1,j) - (*Ez2D)(i-2,j)));
        }
	//(*By2D)(i,0) += dt_ov_dx *( (1-2*beta_y)*((*Ez2D)(i,0) - (*Ez2D)(i-1,0)) +beta_y*((*Ez2D)(i,1)- (*Ez2D)(i-1,1) +(*Ez2D)(i,ny_p-1)- (*Ez2D)(i-1,ny_p-1) ) );
	//(*By2D)(i,ny_p-1) += dt_ov_dx * ( (1-2*beta_y)*((*Ez2D)(i,ny_p-1) - (*Ez2D)(i-1,ny_p-1)) +beta_y*((*Ez2D)(i,0)- (*Ez2D)(i-1,0) +(*Ez2D)(i,ny_p-2)- (*Ez2D)(i-1,ny_p-2) ) );
    
    
    // Magnetic field Bz^(d,d)
    //for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*Bz2D)(i,j) += dt_ov_dy * (Beta_x*( (*Ex2D)(i,j) - (*Ex2D)(i,j-1) ) +beta_x*((*Ex2D)(i+1,j)- (*Ex2D)(i+1,j-1)+ (*Ex2D)(i-1,j)- (*Ex2D)(i-1,j-1)))
            -               dt_ov_dx * (Beta_y*( (*Ey2D)(i,j) - (*Ey2D)(i-1,j) ) +beta_y*((*Ey2D)(i,j+1)- (*Ey2D)(i-1,j+1)+ (*Ey2D)(i,j-1)- (*Ey2D)(i-1,j-1)) +delta_x*( (*Ey2D)(i+1,j) - (*Ey2D)(i-2,j)));
	}
    }
//}// end parallel
}//END solveMaxwellFaraday



