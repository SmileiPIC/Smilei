#include "ElectroMagnBC2D_Long_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn2D.h"
#include "Field2D.h"

using namespace std;

ElectroMagnBC2D_Long_SM::ElectroMagnBC2D_Long_SM( PicParams &params, LaserParams &laser_params )
    : ElectroMagnBC( params, laser_params )
{
    // conversion factor from degree to radian
    conv_deg2rad = M_PI/180.0;
    
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
    
    Bz_xvalmin_Long.resize(ny_p);
    Bz_xvalmax_Long.resize(ny_p);
    By_xvalmin_Long.resize(ny_p);
    By_xvalmax_Long.resize(ny_p);
    Bx_xvalmin_Long.resize(ny_p+1);
    Bx_xvalmax_Long.resize(ny_p+1);

    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------

    // West boundary
    double theta  = 0.0*conv_deg2rad; //0.0; //! \todo Introduce in parameters for Boundary cond., e.g., params.EMBoundary->theta_W
    double factor = 1.0 / (cos(theta) + dt_ov_dx);
    Alpha_SM_W    = 2.0                     * factor;
    Beta_SM_W     = - (cos(theta)-dt_ov_dx) * factor;
    Gamma_SM_W    = 4.0 * cos(theta)        * factor;
    Delta_SM_W    = - (sin(theta)+dt_ov_dy) * factor;
    Epsilon_SM_W  = - (sin(theta)-dt_ov_dy) * factor;

    // East boundary
    theta         = M_PI;
    factor        = 1.0 / (cos(theta) - dt_ov_dx);
    Alpha_SM_E    = 2.0                      * factor;
    Beta_SM_E     = - (cos(theta)+dt_ov_dx)  * factor;
    Gamma_SM_E    = 4.0 * cos(theta)         * factor;
    Delta_SM_E    = - (sin(theta)+dt_ov_dy)  * factor;
    Epsilon_SM_E  = - (sin(theta)-dt_ov_dy)  * factor;

}

ElectroMagnBC2D_Long_SM::~ElectroMagnBC2D_Long_SM()
{
}

void ElectroMagnBC2D_Long_SM::save_fields_BC2D_Long(Field* my_field) {
    Field2D* field2D=static_cast<Field2D*>(my_field);
    
    for (unsigned int j=0 ; j<ny_p ; j++) {
    	if (field2D->name=="Bz"){   
       		Bz_xvalmin_Long[j]=(*field2D)(0,j);
		Bz_xvalmax_Long[j]=(*field2D)(nx_d-1,j);
		}
   	if (field2D->name=="By"){   
      		By_xvalmin_Long[j]=(*field2D)(0,j);
    		By_xvalmax_Long[j]=(*field2D)(nx_d-1,j);
    		}
   	if (field2D->name=="Bx"){   
      		Bx_xvalmin_Long[j]=(*field2D)(0,j);
    		Bx_xvalmax_Long[j]=(*field2D)(nx_p-1,j);
    		}
	}
   	if (field2D->name=="Bx"){ 
    		Bx_xvalmax_Long[ny_p]=(*field2D)(nx_p-1,ny_p);
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_Long_SM::apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi)
{
    // Static cast of the fields
//    Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
    Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
    Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);


    // ! \todo Transverse profile & incidence angle is not yet introduced (MG)
    // -----------------------------------------
    // Laser temporal profile
    // -----------------------------------------
    double byW, bzW, byE, bzE;
//    double dfa; //Distance from axis
//    dfa = 0.0;



    // If !periodic
    //   BC : Bx(i=0...nx_p, 0) & Bx(i=0...nx_p, ny_d-1)
    //   BC : Bz(i=0...nx_d-1, 0) & Bz(i=0...nx_d-1, ny_d-1)
    // else
    //   Ez(i,-1)/Ez(i,ny_d-1) defined -> OK
    //   Ex(i,-1)/Ex(i,ny_d-1) defined -> OK

    // BC : By(0, j=0...ny_p  ) & By(nx_d-1, j=0...ny_p  )	    -> TO DO
    // BC : Bz(O, j=0...ny_d-1) & Bz(nx_d-1, j=0...ny_d-1)		-> TO DO


    // -----------------------------------------
    // Silver-Mueller boundary conditions (West)
    // -----------------------------------------
    if ( smpi->isWestern() ) {
        
        // for By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {

            byW = 0.;
            double yp     = smpi->getDomainLocalMin(1) + ((double)j)     * dy;
            for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
                if (laser_[ilaser]->laser_struct.boxSide == "west") {
                    
                    if ( (laser_[ilaser]->laser_struct.isFocused)||(laser_[ilaser]->laser_struct.angle!=0) ) {
                        byW += 0.0;
                    }
                    else {
                        byW += laser_[ilaser]->a0_delta_y_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual)
                            *  laser_[ilaser]->transverse_profile2D(time_dual,yp);
                    }//isFocused or angle!=0
                    
                }
            }//ilaser

            (*By2D)(0,j) = Alpha_SM_W   * (*Ez2D)(0,j)
                           +              Beta_SM_W    *( (*By2D)(1,j)-By_xvalmin_Long[j])
                           +              Gamma_SM_W   * byW
                           +              Delta_SM_W   *( (*Bx2D)(0,j+1)-Bx_xvalmin_Long[j+1] )
                           +              Epsilon_SM_W *( (*Bx2D)(0,j)-Bx_xvalmin_Long[j] )
			   +		  By_xvalmin_Long[j];
        }
        
        // for Bz^(d,d)
        for (unsigned int j=0 ; j<ny_d ; j++) {

            bzW = 0.;
            double yd     = smpi->getDomainLocalMin(1) + ((double)j-0.5) * dy;
            
            for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
                if (laser_[ilaser]->laser_struct.boxSide == "west") {
                    
                    if ( (laser_[ilaser]->laser_struct.isFocused) || (laser_[ilaser]->laser_struct.angle!=0) ) {
                        double delay   = laser_[ilaser]->laser_struct.delay;
                        double xfoc    = laser_[ilaser]->laser_struct.focus[0];
                        double yfoc    = laser_[ilaser]->laser_struct.focus[1];
                        double theta   = laser_[ilaser]->laser_struct.angle * conv_deg2rad;
                        double zeta    = -xfoc*cos(theta) + (yd-yfoc)*sin(theta);
                        double rho     =  xfoc*sin(theta) + (yd-yfoc)*cos(theta);
                        double tau     = time_dual - yd*sin(theta) - delay;
                        double bwaist  = 0.5/sqrt(log(2.0)) * laser_[ilaser]->laser_struct.profile_transv.double_params[0];
                        double z2ovLr2 = pow(zeta,2)/pow(bwaist,4);
                        double waist   = bwaist * sqrt( 1.0 + z2ovLr2 );
                        double curvRad = 1000.0 * laser_[ilaser]->laser_struct.profile_transv.double_params[0];
                        if (zeta!=0)
                               curvRad = zeta* ( 1.0 + 1.0/z2ovLr2 );
                        double gouyPhs = -0.5 * atan( sqrt(z2ovLr2) );
                        double phi     = -0.5 * pow(rho,2)/curvRad + gouyPhs;
                        bzW += laser_[ilaser]->laser_struct.a0 * sin(tau+phi) * laser_[ilaser]->time_profile(tau)
                            *  laser_[ilaser]->transverse_profile2D(time_dual,rho/waist);
                    }
                    else {
                        bzW += laser_[ilaser]->a0_delta_z_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual)
                            *  laser_[ilaser]->transverse_profile2D(time_dual,yd);
                    }//isFocused or angle!=0
                }
            }//ilaser

            (*Bz2D)(0,j) = -Alpha_SM_W * (*Ey2D)(0,j)
                           +               Beta_SM_W  *( (*Bz2D)(1,j)- Bz_xvalmin_Long[j])
                           +               Gamma_SM_W * bzW
		 	   + 	   	   Bz_xvalmin_Long[j];
        }
    }//if West

    // -----------------------------------------
    // Silver-Mueller boundary conditions (East)
    // -----------------------------------------
    if ( smpi->isEastern() ) {
        // for By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            
            byE = 0.;
            double yp     = smpi->getDomainLocalMin(1) + ((double)j)     * dy;
            for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
                    // Incident field (west boundary)
                if (laser_[ilaser]->laser_struct.boxSide == "east") {
                    byE += laser_[ilaser]->a0_delta_y_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual) * laser_[ilaser]->transverse_profile2D(time_dual,yp);
                }
            }//ilaser

            (*By2D)(nx_d-1,j) = Alpha_SM_E   * (*Ez2D)(nx_p-1,j)
                                +                   Beta_SM_E    *( (*By2D)(nx_d-2,j) -By_xvalmax_Long[j])
                                +                   Gamma_SM_E   * byE
                                +                   Delta_SM_E   *( (*Bx2D)(nx_p-1,j+1) -Bx_xvalmax_Long[j+1])// Check x-index
                                +                   Epsilon_SM_E *( (*Bx2D)(nx_p-1,j) -Bx_xvalmax_Long[j])
				+ 	 	    By_xvalmax_Long[j];
        }
        
        // for Bz^(d,d)
        for (unsigned int j=0 ; j<ny_d ; j++) {

            bzE = 0.;
            double yd     = smpi->getDomainLocalMin(1) + ((double)j-0.5) * dy;
            for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
                if (laser_[ilaser]->laser_struct.boxSide == "east") {
                    // Incident field (east boundary)
                    bzE += laser_[ilaser]->a0_delta_z_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual) * laser_[ilaser]->transverse_profile2D(time_dual,yd);
                }
            }//ilaser

            (*Bz2D)(nx_d-1,j) = -Alpha_SM_E * (*Ey2D)(nx_p-1,j)
                                +                    Beta_SM_E  *( (*Bz2D)(nx_d-2,j) -Bz_xvalmax_Long[j])
                                +                    Gamma_SM_E * bzE
				+ 	 	     Bz_xvalmax_Long[j];
        }
    }//if East


}// END apply

