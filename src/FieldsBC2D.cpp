#include "FieldsBC2D.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn2D.h"
#include "Field2D.h"
#include "Laser.h"
#include "Tools.h"

using namespace std;

FieldsBC2D::FieldsBC2D( PicParams *params )
    : FieldsBC( params )
{
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params->n_space[0]+1+2*params->oversize[0];
    nx_d = params->n_space[0]+2+2*params->oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params->n_space[1]+1+2*params->oversize[1];
    ny_d = params->n_space[1]+2+2*params->oversize[1];

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = params->cell_length[0];
    dt_ov_dx = dt/dx;
    dx_ov_dt = 1.0/dt_ov_dx;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = params->cell_length[1];
    dt_ov_dy = dt/dy;
    dy_ov_dt = 1.0/dt_ov_dy;

    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------

    // West boundary
    double theta  = 0.0; //! \todo Introduce in parameters for Boundary cond., e.g., params->EMBoundary->theta_W
    double factor = 1.0 / (cos(theta) + dt_ov_dx);
    Alpha_SM_W    = 2.0                     * factor;
    Beta_SM_W     = - (cos(theta)-dt_ov_dx) * factor;
    Gamma_SM_W    = 4.0 * cos(theta)        * factor;
    Delta_SM_W    = - (sin(theta)+dt_ov_dy) * factor;
    Epsilon_SM_W  = - (sin(theta)-dt_ov_dy) * factor;
    MESSAGE("WEST : " << Alpha_SM_W << Beta_SM_W << Gamma_SM_W);

    // East boundary
    theta         = M_PI;
    factor        = 1.0 / (cos(theta) - dt_ov_dx);
    Alpha_SM_E    = 2.0                      * factor;
    Beta_SM_E     = - (cos(theta)+dt_ov_dx)  * factor;
    Gamma_SM_E    = 4.0 * cos(theta)         * factor;
    Delta_SM_E    = - (sin(theta)+dt_ov_dy)  * factor;
    Epsilon_SM_E  = - (sin(theta)-dt_ov_dy)  * factor;
    MESSAGE("EAST : " << Alpha_SM_E << Beta_SM_E << Gamma_SM_E);
}

FieldsBC2D::~FieldsBC2D()
{

}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void FieldsBC2D::apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi)
{
    // Static cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
    Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
    Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);


    // ! \todo Transverse profile & incidence angle is not yet introduced (MG)
    // -----------------------------------------
    // Laser temporal profile
    // -----------------------------------------
    double byW=0.0, bzW=0.0, byE=0.0, bzE=0.0;
    double dfa; //Distance from axis
    dfa = 0.0;



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
    if ( smpi->isWester() ) {
        // for By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {

            //dfa = smpi->getDomainLocalMin(1)+j*params.cell_length[1]-params.sim_length[1]/2. ; //dfa is algebric.
            for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
                if (laser_[ilaser]->laser_struct.angle == 0) {
                    // Incident field (west boundary)
                    byW += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual) * laser_[ilaser]->y_profile(dfa);
                }
            }//ilaser

            (*By2D)(0,j) = Alpha_SM_W   * (*Ez2D)(0,j)
                           +              Beta_SM_W    * (*By2D)(1,j)
                           +              Gamma_SM_W   * byW
                           +              Delta_SM_W   * (*Bx2D)(0,j+1)
                           +              Epsilon_SM_W * (*Bx2D)(0,j);
        }
        // for Bz^(d,d)
        for (unsigned int j=0 ; j<ny_d ; j++) {

            //dfa = smpi->getDomainLocalMin(1)+j*params.cell_length[1]-params->sim_length[1]/2. ; //dfa is algebric.
            for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
                if (laser_[ilaser]->laser_struct.angle == 0) {
                    // Incident field (west boundary)
                    bzW += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual) * laser_[ilaser]->y_profile(dfa);
                }
            }//ilaser

            (*Bz2D)(0,j) = -Alpha_SM_W * (*Ey2D)(0,j)
                           +               Beta_SM_W  * (*Bz2D)(1,j)
                           +               Gamma_SM_W * bzW;
        }
    }//if West

    // -----------------------------------------
    // Silver-Mueller boundary conditions (East)
    // -----------------------------------------
    if ( smpi->isEaster() ) {
        // for By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            
            //dfa = smpi->getDomainLocalMin(1)+j*params.cell_length[1]-params->sim_length[1]/2. ; //dfa is algebric.
            for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
                    // Incident field (west boundary)
                if (laser_[ilaser]->laser_struct.angle == 180) {
                    byE += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual) * laser_[ilaser]->y_profile(dfa);
                }
            }//ilaser

            (*By2D)(nx_d-1,j) = Alpha_SM_E   * (*Ez2D)(nx_p-1,j)
                                +                   Beta_SM_E    * (*By2D)(nx_d-2,j)
                                +                   Gamma_SM_E   * byE
                                +                   Delta_SM_E   * (*Bx2D)(nx_p-1,j+1) // Check x-index
                                +                   Epsilon_SM_E * (*Bx2D)(nx_p-1,j);
        }
        // for Bz^(d,d)
        for (unsigned int j=0 ; j<ny_d ; j++) {

            //dfa = smpi->getDomainLocalMin(1)+j*params.cell_length[1]-params->sim_length[1]/2. ; //dfa is algebric.
            for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
                if (laser_[ilaser]->laser_struct.angle == 180) {
                    // Incident field (east boundary)
                    bzE += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual) * laser_[ilaser]->y_profile(dfa);
                }
            }//ilaser

            (*Bz2D)(nx_d-1,j) = -Alpha_SM_E * (*Ey2D)(nx_p-1,j)
                                +                    Beta_SM_E  * (*Bz2D)(nx_d-2,j)
                                +                    Gamma_SM_E * bzE;
        }
    }//if East

}// END apply

