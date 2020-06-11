#include "MF_Solver2D_Bouchard.h"

#include "ElectroMagn.h"
#include "Field2D.h"

#include <algorithm>

MF_Solver2D_Bouchard::MF_Solver2D_Bouchard(Params &params)
    : Solver2D(params)
{
    //ERROR("Under development, not yet working");
    double dt = params.timestep;
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    double dx_ov_dt  = dx/dt;
    double dy_ov_dt  = dy/dt;
    double dt_ov_dx  = dt/dx;
    double dt_ov_dy  = dt/dy;
    if( dx!=dy ) {
        ERROR( "Bouchard solver requires the same cell-length in x and y directions" );
    }
    if( dx_ov_dt!=2 ) {
        WARNING( "Bouchard solver requires dx/dt = 2 (Magical Timestep)" );
    }
    
    // On the axes v_phi^max = 1.01c and is below c @ 0.54 kxdx/pi
    // So there could existe a numerical cherenkov emission at this point
    // On the diagonal v_phi^max = 1.01c and is below c @ 0.85 sqrt((kxdx)^2+(kydy)^2)
    double delta = 0.1222*(1-pow(2.,2))/4. ;
    double beta = -0.1727*(1-0.5*pow(2.,2)-4.*delta)/4. ;
    double alpha = 1-2.*beta-3.*delta;
 
    beta_xy = beta;
    beta_yx = beta;
    delta_y = delta;
    delta_x = delta;
   
    alpha_y =  alpha;
    alpha_x =  alpha;
   
    Ax    = alpha_x*dt/dx;
    Ay    = alpha_y*dt/dy;
    Bx    = beta_xy*dt/dx;
    By    = beta_yx*dt/dy;
    Dx    = delta_x*dt/dy;
    Dy    = delta_y*dt/dy;

    isEFilterApplied = false;
    if (params.Friedman_filter)
        isEFilterApplied = true;
}

MF_Solver2D_Bouchard::~MF_Solver2D_Bouchard()
{
}

void MF_Solver2D_Bouchard::operator() ( ElectroMagn* fields )
{
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
        Ex2D = static_cast<Field2D*>(fields->Exfilter[0]);
        Ey2D = static_cast<Field2D*>(fields->Eyfilter[0]);
        Ez2D = static_cast<Field2D*>(fields->Ezfilter[0]);
    } else {
        Ex2D = static_cast<Field2D*>(fields->Ex_);
        Ey2D = static_cast<Field2D*>(fields->Ey_);
        Ez2D = static_cast<Field2D*>(fields->Ez_);
    }
    Field2D* Bx2D = static_cast<Field2D*>(fields->Bx_);
    Field2D* By2D = static_cast<Field2D*>(fields->By_);
    Field2D* Bz2D = static_cast<Field2D*>(fields->Bz_);


    // Magnetic field Bx^(p,d)
    for (unsigned int i=1 ; i<nx_p-1;  i++) {
        for (unsigned int j=2 ; j<ny_d-2 ; j++) {
            (*Bx2D)(i,j) += Ay * ((*Ez2D)(i,j-1)-(*Ez2D)(i,j))
                          + By * ((*Ez2D)(i+1,j-1)-(*Ez2D)(i+1,j) + (*Ez2D)(i-1,j-1)-(*Ez2D)(i-1,j))
                          + Dy * ((*Ez2D)(i,j-2)-(*Ez2D)(i,j+1));
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
    for (unsigned int j=2 ; j<ny_d-2 ; j++) {
        (*Bz2D)(1,j) += dt_ov_dx * ( (*Ey2D)(0,j) - (*Ey2D)(1,j)   )
        +               dt_ov_dy * ( (*Ex2D)(1,j) - (*Ex2D)(1,j-1) );
    }
    // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
    for (unsigned int j=2 ; j<ny_d-2 ; j++) {
        (*Bz2D)(nx_d-2,j) += dt_ov_dx * ( (*Ey2D)(nx_d-3,j) - (*Ey2D)(nx_d-2,j)   )
        +                    dt_ov_dy * ( (*Ex2D)(nx_d-2,j) - (*Ex2D)(nx_d-2,j-1) );
    }
    
//}// end parallel
}//END solveMaxwellFaraday



