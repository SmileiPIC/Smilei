#include "ElectroMagnBC2D_Trans_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn2D.h"
#include "Field2D.h"

using namespace std;

ElectroMagnBC2D_Trans_SM::ElectroMagnBC2D_Trans_SM( PicParams &params, LaserParams &laser_params ) 
    : ElectroMagnBC( params, laser_params )
{
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params.n_space[0]+1+2*params.oversize[0];
    nx_d = params.n_space[0]+2+2*params.oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params.n_space[1]+1+2*params.oversize[1];
    ny_d = params.n_space[1]+2+2*params.oversize[1];

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = params.cell_length[0];
    dt_ov_dx = dt/dx;
    dx_ov_dt = 1.0/dt_ov_dx;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = params.cell_length[1];
    dt_ov_dy = dt/dy;
    dy_ov_dt = 1.0/dt_ov_dy;
    
    Bz_yvalmin_Trans.resize(nx_p);
    Bz_yvalmax_Trans.resize(nx_p);
    By_yvalmin_Trans.resize(nx_p+1);
    By_yvalmax_Trans.resize(nx_p+1);
    Bx_yvalmin_Trans.resize(nx_p);
    Bx_yvalmax_Trans.resize(nx_p);

    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------

    // South boundary
    double theta  = 0.0; 
    double factor = 1.0 / (cos(theta) + dt_ov_dy );
    Alpha_SM_S    = 2.0                     * factor;
    Beta_SM_S     = - (cos(theta)-dt_ov_dy) * factor;
    Delta_SM_S    = - (sin(theta)+dt_ov_dx) * factor;
    Epsilon_SM_S  = - (sin(theta)-dt_ov_dx) * factor;
    // North boundary
    theta  = M_PI; 
    factor = 1.0 / (cos(theta) - dt_ov_dy);
    Alpha_SM_N    = 2.0                     * factor;
    Beta_SM_N     = - (cos(theta)+dt_ov_dy) * factor;
    Delta_SM_N    = - (sin(theta)+dt_ov_dx) * factor;
    Epsilon_SM_N  = - (sin(theta)-dt_ov_dx) * factor;

}

ElectroMagnBC2D_Trans_SM::~ElectroMagnBC2D_Trans_SM()
{
}


void ElectroMagnBC2D_Trans_SM::save_fields_BC2D_Trans(Field* my_field) {
    Field2D* field2D=static_cast<Field2D*>(my_field);
    
    for (unsigned int j=0 ;  j<nx_p ; j++) {
    	if (field2D->name=="Bz"){   
       		Bz_yvalmin_Trans[j]=(*field2D)(j,0);
		Bz_yvalmax_Trans[j]=(*field2D)(j,ny_d-1);
		}
		
    	if (field2D->name=="By"){   
       		By_yvalmin_Trans[j]=(*field2D)(j,0);
		By_yvalmax_Trans[j]=(*field2D)(j,ny_p-1);
		}
		
    	if (field2D->name=="Bx"){  
        	Bx_yvalmin_Trans[j]=(*field2D)(j,0);
    		Bx_yvalmax_Trans[j]=(*field2D)(j,ny_d-1);
    		}
	}
    	if (field2D->name=="By"){  
		By_yvalmax_Trans[nx_p]=(*field2D)(nx_p,ny_p-1);
	}
	
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_Trans_SM::apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi)
{
    // Static cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
//    Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
    Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
    Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);

    // -----------------------------------------
    // Silver-Mueller boundary conditions (South)
    // -----------------------------------------
    if ( smpi->isSouthern() ) {
   // MESSAGE("BC south : " << Bx_yvalmin_Trans[0])
        // for Bx^(p,d)
        for (unsigned int j=0 ; j<nx_p ; j++) {

            (*Bx2D)(j,0) = -Alpha_SM_S   * (*Ez2D)(j,0)
                           +              Beta_SM_S    *( (*Bx2D)(j,1)-Bx_yvalmin_Trans[j])
                           +              Delta_SM_S   *( (*By2D)(j+1,0)-By_yvalmin_Trans[j+1])
                           +              Epsilon_SM_S *( (*By2D)(j,0)-By_yvalmin_Trans[j])
			   +		  Bx_yvalmin_Trans[j];
        }
        // for Bz^(d,d)
        for (unsigned int j=0 ; j<nx_d ; j++) {

            (*Bz2D)(j,0) = Alpha_SM_S * (*Ex2D)(j,0)
                           +               Beta_SM_S  *( (*Bz2D)(j,1)-Bz_yvalmin_Trans[j])
			   +		   Bz_yvalmin_Trans[j];
        }
    }//if South

    // -----------------------------------------
    // Silver-Mueller boundary conditions (North)
    // -----------------------------------------
   // MESSAGE("BC north : " << Bx_yvalmax_Trans[0])
    if ( smpi->isNorthern() ) {
        // for Bx^(p,d)
        for (unsigned int j=0 ; j<nx_p ; j++) {
            (*Bx2D)(j,ny_d-1) = -Alpha_SM_N   * (*Ez2D)(j,ny_p-1)
                                +                   Beta_SM_N    *( (*Bx2D)(j,ny_d-2) -Bx_yvalmax_Trans[j])
                                +                   Delta_SM_N   *( (*By2D)(j+1,ny_p-1) -By_yvalmax_Trans[j+1])
                                +                   Epsilon_SM_N *( (*By2D)(j,ny_p-1) -By_yvalmax_Trans[j])
				+		    Bx_yvalmax_Trans[j];
        }
        // for Bz^(d,d)
        for (unsigned int j=0 ; j<nx_d ; j++) {

            (*Bz2D)(j,ny_d-1) = Alpha_SM_N * (*Ex2D)(j,ny_p-1)
                                +                    Beta_SM_N  *( (*Bz2D)(j,ny_d-2)- Bz_yvalmax_Trans[j])
			        +		     Bz_yvalmax_Trans[j];
        }
    }//if North

}// END apply

