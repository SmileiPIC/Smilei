#include "MF_Solver2D_M4.h"

#include "ElectroMagn.h"
#include "Field2D.h"

#include <algorithm>

MF_Solver2D_M4::MF_Solver2D_M4(Params &params)
    : Solver2D(params)
{
    
    beta_yx = dt_ov_dx * dt_ov_dx / 12.0;
    beta_xy = dt_ov_dy * dt_ov_dy / 12.0;
    delta_x = beta_yx - 1.0 / 12.0;
    delta_y = beta_xy - 1.0 / 12.0;
   
    alpha_x =  1.0 - 2.0 * beta_xy - 3.0 * delta_x;
    alpha_y =  1.0 - 2.0 * beta_yx - 3.0 * delta_y;
   
    Ax    = alpha_x*dt_ov_dx;
    Ay    = alpha_y*dt_ov_dy;
    Bx    = beta_xy*dt_ov_dx;
    By    = beta_yx*dt_ov_dy;
    Dx    = delta_x*dt_ov_dx;
    Dy    = delta_y*dt_ov_dy;

    isEFilterApplied = params.Friedman_filter;
}

MF_Solver2D_M4::~MF_Solver2D_M4()
{
}

void MF_Solver2D_M4::operator() ( ElectroMagn* fields )
{
    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    // Static-cast of the fields
    // Field2D* Ex2D = static_cast<Field2D*>(fields->Ex_);
    // Field2D* Ey2D = static_cast<Field2D*>(fields->Ey_);
    // Field2D* Ez2D = static_cast<Field2D*>(fields->Ez_);
    // Field2D* Bx2D = static_cast<Field2D*>(fields->Bx_);
    // Field2D* By2D = static_cast<Field2D*>(fields->By_);
    // Field2D* Bz2D = static_cast<Field2D*>(fields->Bz_);

    // Static-cast of the fields
    Field2D* Ex2D;
    Field2D* Ey2D;
    Field2D* Ez2D;
    if (isEFilterApplied) {
        Ex2D = static_cast<Field2D*>(fields->filter_->Ex_[0]);
        Ey2D = static_cast<Field2D*>(fields->filter_->Ey_[0]);
        Ez2D = static_cast<Field2D*>(fields->filter_->Ez_[0]);
    } else {
        Ex2D = static_cast<Field2D*>(fields->Ex_);
        Ey2D = static_cast<Field2D*>(fields->Ey_);
        Ez2D = static_cast<Field2D*>(fields->Ez_);
    }
    Field2D* Bx2D = static_cast<Field2D*>(fields->Bx_);
    Field2D* By2D = static_cast<Field2D*>(fields->By_);
    Field2D* Bz2D = static_cast<Field2D*>(fields->Bz_);


    // Magnetic field Bx^(p,d)
    {
        for (unsigned int j=1 ; j<ny_p ; j++) {
            (*Bx2D)(0,j) += dt_ov_dy * ((*Ez2D)(0,j-1)-(*Ez2D)(0,j));
        }
    }
    
    // Magnetic field Bx^(p,d)
    for (unsigned int i=1 ; i<nx_p-1;  i++) {
        for (unsigned int j=2 ; j<ny_d-2 ; j++) {
            (*Bx2D)(i,j) += Ay * ((*Ez2D)(i,j-1)-(*Ez2D)(i,j))
                          + By * ((*Ez2D)(i+1,j-1)-(*Ez2D)(i+1,j) + (*Ez2D)(i-1,j-1)-(*Ez2D)(i-1,j))
                          + Dy * ((*Ez2D)(i,j-2)-(*Ez2D)(i,j+1));
        }
    }

    // Magnetic field Bx^(p,d)
    {
        for (unsigned int j=1 ; j<ny_p ; j++) {
            (*Bx2D)(nx_p-1,j) += dt_ov_dy * ((*Ez2D)(nx_p-1,j-1)-(*Ez2D)(nx_p-1,j));
        }
    }
    
    // Magnetic field By^(d,p)
    for (unsigned int i=2 ; i<nx_d-2 ; i++) {
        for (unsigned int j=1 ; j<ny_p-1 ; j++) {
            (*By2D)(i,j) += Ax * ((*Ez2D)(i,j) - (*Ez2D)(i-1,j))
                          + Bx * ((*Ez2D)(i,j+1)-(*Ez2D)(i-1,j+1) + (*Ez2D)(i,j-1)-(*Ez2D)(i-1,j-1))
                          + Dx * ((*Ez2D)(i+1,j) - (*Ez2D)(i-2,j));
        }
    }
       
    // Magnetic field Bz^(d,d)
    for (unsigned int i=2 ; i<nx_d-2 ; i++) {
        for (unsigned int j=2 ; j<ny_d-2 ; j++) {
            (*Bz2D)(i,j) += Ay * ((*Ex2D)(i,j)-(*Ex2D)(i,j-1))
                          + By * ((*Ex2D)(i+1,j)-(*Ex2D)(i+1,j-1) + (*Ex2D)(i-1,j)-(*Ex2D)(i-1,j-1))
                          + Dy * ((*Ex2D)(i,j+1)-(*Ex2D)(i,j-2))
                          + Ax * ((*Ey2D)(i-1,j)-(*Ey2D)(i,j))
                          + Bx * ((*Ey2D)(i-1,j+1)-(*Ey2D)(i,j+1) + (*Ey2D)(i-1,j-1)-(*Ey2D)(i,j-1))
                          + Dx * ((*Ey2D)(i-2,j)-(*Ey2D)(i+1,j));
        }
    }
    
    // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
    for (unsigned int j=0 ; j<ny_p ; j++) {
        (*By2D)(1,j) += dt_ov_dx * ( (*Ez2D)(1,j) - (*Ez2D)(0,j) );
    }
    // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
    for (unsigned int j=0 ; j<ny_p ; j++) {
        (*By2D)(nx_d-2,j) += dt_ov_dx * ( (*Ez2D)(nx_d-2,j) - (*Ez2D)(nx_d-3,j) );
    }
    
    // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
    for (unsigned int j=1 ; j<ny_p ; j++) {
        (*Bz2D)(1,j) += dt_ov_dx * ( (*Ey2D)(0,j) - (*Ey2D)(1,j)   )
        +               dt_ov_dy * ( (*Ex2D)(1,j) - (*Ex2D)(1,j-1) );
    }
    // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
    for (unsigned int j=1 ; j<ny_p ; j++) {
        (*Bz2D)(nx_d-2,j) += dt_ov_dx * ( (*Ey2D)(nx_d-3,j) - (*Ey2D)(nx_d-2,j)   )
        +                    dt_ov_dy * ( (*Ex2D)(nx_d-2,j) - (*Ex2D)(nx_d-2,j-1) );
    }
    
//}// end parallel
}//END solveMaxwellFaraday



