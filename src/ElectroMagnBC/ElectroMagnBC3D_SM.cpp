#include "ElectroMagnBC3D_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field3D.h"
#include "Tools.h"
#include "Laser.h"

using namespace std;

ElectroMagnBC3D_SM::ElectroMagnBC3D_SM( Params &params, Patch* patch )
  : ElectroMagnBC( params, patch )
{
    // conversion factor from degree to radian
    conv_deg2rad = M_PI/180.0;
    
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params.n_space[0]+1+2*params.oversize[0];
    nx_d = nx_p+1;
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params.n_space[1]+1+2*params.oversize[1];
    ny_d = ny_p+1;
    // number of nodes of the primal and dual grid in the z-direction
    nz_p = params.n_space[2]+1+2*params.oversize[2];
    nz_d = nz_p+1;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = params.cell_length[0];
    dt_ov_dx = dt/dx;
    dx_ov_dt = 1.0/dt_ov_dx;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = params.cell_length[1];
    dt_ov_dy = dt/dy;
    dy_ov_dt = 1.0/dt_ov_dy;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the z-direction)
    dz       = params.cell_length[2];
    dt_ov_dz = dt/dz;
    dz_ov_dt = 1.0/dt_ov_dz;


    std::vector<unsigned int> dims(2,0);

    // BCs at the x-border
    dims = { ny_d, nz_d }; // Bx^(p,d,d)
    Bx_xvalmin_Long = new Field2D(dims); 
    Bx_xvalmax_Long = new Field2D(dims);
    dims = { ny_p, nz_d }; // By^(d,p,d)
    By_xvalmin_Long = new Field2D(dims);
    By_xvalmax_Long = new Field2D(dims);
    dims = { ny_d, nz_p }; // Bz^(d,d,p)
    Bz_xvalmin_Long = new Field2D(dims);
    Bz_xvalmax_Long = new Field2D(dims);
    
    // BCs in the y-border
    dims = { nx_p, nz_d }; // Bx^(p,d,d)
    Bx_yvalmin_TransY = new Field2D(dims);
    Bx_yvalmax_TransY = new Field2D(dims);
    dims = { nx_d, nz_d }; // By^(d,p,d)
    By_yvalmin_TransY = new Field2D(dims);
    By_yvalmax_TransY = new Field2D(dims);
    dims = { nx_d, nz_p }; // Bz^(d,d,p)
    Bz_yvalmin_TransY = new Field2D(dims);
    Bz_yvalmax_TransY = new Field2D(dims);
    
    // BCs in the z-border
    dims = { nx_p, ny_d }; // Bx^(p,d,d)
    Bx_zvalmin_TransZ = new Field2D(dims);
    Bx_zvalmax_TransZ = new Field2D(dims);
    dims = { nx_d, ny_p }; // By^(d,p,d)
    By_zvalmin_TransZ = new Field2D(dims);
    By_zvalmax_TransZ = new Field2D(dims);
    dims = { nx_d, ny_d }; // Bz^(d,d,p)
    Bz_zvalmin_TransZ = new Field2D(dims);
    Bz_zvalmax_TransZ = new Field2D(dims);
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    //! \todo (MG) Check optimal angle for Silver-Muller BCs
    
    // West boundary
    double theta  = 0.0*conv_deg2rad; //0.0;
    double factor = 1.0 / (cos(theta) + dt_ov_dx);
    Alpha_SM_W    = 2.0                     * factor;
    Beta_SM_W     = - (cos(theta)-dt_ov_dx) * factor;
    Gamma_SM_W    = 4.0 * cos(theta)        * factor;
    Delta_SM_W    = - (sin(theta)+dt_ov_dy) * factor;
    Epsilon_SM_W  = - (sin(theta)-dt_ov_dy) * factor;
    Zeta_SM_W     = - dt_ov_dz              * factor;
    Eta_SM_W      =   dt_ov_dz              * factor;
    
    // East boundary
    theta         = M_PI;
    factor        = 1.0 / (cos(theta) - dt_ov_dx);
    Alpha_SM_E    = 2.0                      * factor;
    Beta_SM_E     = - (cos(theta)+dt_ov_dx)  * factor;
    Gamma_SM_E    = 4.0 * cos(theta)         * factor;
    Delta_SM_E    = - (sin(theta)+dt_ov_dy)  * factor;
    Epsilon_SM_E  = - (sin(theta)-dt_ov_dy)  * factor;
    Zeta_SM_E     = - dt_ov_dz              * factor;
    Eta_SM_E      =   dt_ov_dz              * factor;
    
    // South boundary
    theta  = 0.0;
    factor = 1.0 / (cos(theta) + dt_ov_dy );
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

#ifdef _PATCH3D_TODO
    // Bottom boundary
    theta  = 0.0;
    factor = 1.0 / (cos(theta) + dt_ov_dz );
    Alpha_SM_B    = 0.;
    Beta_SM_B     = 0.;
    Delta_SM_B    = 0.;
    Epsilon_SM_B  = 0.;
    
    // Top boundary
    theta  = M_PI;
    factor = 1.0 / (cos(theta) - dt_ov_dz);
    Alpha_SM_T    = 0.;
    Beta_SM_T     = 0.;
    Delta_SM_T    = 0.;
    Epsilon_SM_T  = 0.;
#endif

}

ElectroMagnBC3D_SM::~ElectroMagnBC3D_SM()
{
}


// Magnetic field Bx^(p,d,d)
// Magnetic field By^(d,p,d)
// Magnetic field Bz^(d,d,p)

void ElectroMagnBC3D_SM::save_fields_BC3D_Long(Field* my_field) {
    Field3D* field3D=static_cast<Field3D*>(my_field);
    
    if (field3D->name=="Bx"){
        field3D->extract_slice_yz(0,      Bx_xvalmin_Long);
        field3D->extract_slice_yz(nx_p-1, Bx_xvalmax_Long);
    }
    
    if (field3D->name=="By"){
        field3D->extract_slice_yz(0,      By_xvalmin_Long);
        field3D->extract_slice_yz(nx_d-1, By_xvalmax_Long);
    }
    
    if (field3D->name=="Bz"){
        field3D->extract_slice_yz(0,      Bz_xvalmin_Long);
        field3D->extract_slice_yz(nx_d-1, Bz_xvalmax_Long);
    }
    
}

void ElectroMagnBC3D_SM::save_fields_BC3D_TransY(Field* my_field) {
    Field3D* field3D=static_cast<Field3D*>(my_field);
  
    if (field3D->name=="Bx"){
        field3D->extract_slice_xz(0,      Bx_yvalmin_TransY);
        field3D->extract_slice_xz(ny_d-1, Bx_yvalmax_TransY);
    }
    
    if (field3D->name=="By"){
        field3D->extract_slice_xz(0,      By_yvalmin_TransY);
        field3D->extract_slice_xz(ny_p-1, By_yvalmax_TransY);
    }
    
    if (field3D->name=="Bz"){
        field3D->extract_slice_xz(0,      Bz_yvalmin_TransY);
        field3D->extract_slice_xz(ny_d-1, Bz_yvalmax_TransY);
    }
    
}

void ElectroMagnBC3D_SM::save_fields_BC3D_TransZ(Field* my_field) {
    Field3D* field3D=static_cast<Field3D*>(my_field);
  
    if (field3D->name=="Bx"){
        field3D->extract_slice_xy(0,      Bx_zvalmin_TransZ);
        field3D->extract_slice_xy(nz_d-1, Bx_zvalmax_TransZ);
    }
    
    if (field3D->name=="By"){
        field3D->extract_slice_xy(0,      By_zvalmin_TransZ);
        field3D->extract_slice_xy(nz_d-1, By_zvalmax_TransZ);
    }
    
    if (field3D->name=="Bz"){
        field3D->extract_slice_xy(0,      Bz_zvalmin_TransZ);
        field3D->extract_slice_xy(nz_p-1, Bz_zvalmax_TransZ);
    }
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_SM::apply_xmin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isWestern() ) {
        
        // Static cast of the fields
        //Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        vector<double> pos(2);

        // for By^(d,p,d) 
        pos[0] = patch->getDomainLocalMin(1) - EMfields->oversize[1]*dy;
        for (unsigned int j=0 ; j<ny_p ; j++) {
             pos[0] += dy;
             pos[1] = patch->getDomainLocalMin(2) - (0.5 + EMfields->oversize[2])*dz;
             for (unsigned int k=0 ; k<nz_d ; k++) {
                 pos[1] += dz;
                 // Lasers
                 double byW = 0.;
                 for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                     byW += vecLaser[ilaser]->getAmplitude0(pos, time_dual, j, k);
                 }

                 (*By3D)(0,j,k) = Alpha_SM_W   * (*Ez3D)(0,j,k)
                 +              Beta_SM_W    *( (*By3D)(1,j,k)-(*By_xvalmin_Long)(j,k))
                 +              Gamma_SM_W   * byW
                 +              Delta_SM_W   *( (*Bx3D)(0,j+1,k)-(*Bx_xvalmin_Long)(j+1,k) )
                 +              Epsilon_SM_W *( (*Bx3D)(0,j,k)-(*Bx_xvalmin_Long)(j,k) )
                 +              (*By_xvalmin_Long)(j,k);
             }// k  ---end compute By
         }//j  ---end compute By

        // for Bz^(d,d,p)
        pos[0] = patch->getDomainLocalMin(1) - (0.5 + EMfields->oversize[1])*dy;
        for (unsigned int j=0 ; j<ny_d ; j++) {
             pos[0] += dy;
             pos[1] = patch->getDomainLocalMin(2) - EMfields->oversize[2]*dz;
             for (unsigned int k=0 ; k<nz_p ; k++) {
                 pos[1] += dz;
                 // Lasers
                 double bzW = 0.;
                 for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                     bzW += vecLaser[ilaser]->getAmplitude0(pos, time_dual, j, k);
                 }
                 (*Bz3D)(0,j,k) = - Alpha_SM_W   * (*Ey3D)(0,j,k)
                 +              Beta_SM_W    *( (*Bz3D)(1,j,k)-(*Bz_xvalmin_Long)(j,k))
                 +              Gamma_SM_W   * bzW
                 +              Zeta_SM_W   *( (*Bx3D)(0,j,k+1)-(*Bx_xvalmin_Long)(j,k+1) )
                 +              Eta_SM_W *( (*Bx3D)(0,j,k)-(*Bx_xvalmin_Long)(j,k) )
                 +              (*Bz_xvalmin_Long)(j,k);

             }// k  ---end compute Bz
         }//j  ---end compute Bz

        
    }//if Western

}
// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_SM::apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isEastern() ) {
        
        // Static cast of the fields
        //Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        vector<double> pos(2);

        // for By^(d,p,d)
        pos[0] = patch->getDomainLocalMin(1) - EMfields->oversize[1]*dy;
        for (unsigned int j=0 ; j<ny_p ; j++) {
             pos[0] += dy;
             pos[1] = patch->getDomainLocalMin(2) - (0.5 + EMfields->oversize[2])*dz;
             for (unsigned int k=0 ; k<nz_d ; k++) {
                 pos[1] += dz;
                // Lasers
                double byE = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    byE += vecLaser[ilaser]->getAmplitude0(pos, time_dual, j, 0);
                }
            
                (*By3D)(nx_d-1,j,k) = Alpha_SM_E   * (*Ez3D)(nx_p-1,j,k)
                +                   Beta_SM_E    *( (*By3D)(nx_d-2,j,k) -(*By_xvalmax_Long)(j,k))
                +                   Gamma_SM_E   * byE
                +                   Delta_SM_E   *( (*Bx3D)(nx_p-1,j+1,k) -(*Bx_xvalmax_Long)(j+1,k))// Check x-index
                +                   Epsilon_SM_E *( (*Bx3D)(nx_p-1,j,k) -(*Bx_xvalmax_Long)(j,k))
                +                   (*By_xvalmax_Long)(j,k);
            
            }//k  ---end compute By
        }//j  ---end compute By
        
        
        // for Bz^(d,d,p)
        pos[0] = patch->getDomainLocalMin(1) - (0.5+EMfields->oversize[1])*dy;
        for (unsigned int j=0 ; j<ny_d ; j++) {
            pos[1] = patch->getDomainLocalMin(2) - EMfields->oversize[2]*dz;
            pos[0] += dy;
            for (unsigned int k=0 ; k<nz_d ; k++) {
                pos[1] += dz;
                // Lasers
                double bzE = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    bzE += vecLaser[ilaser]->getAmplitude1(pos, time_dual, j, 0);
                }
                
                (*Bz3D)(nx_d-1,j,k) = -Alpha_SM_E * (*Ey3D)(nx_p-1,j,k)
                +                    Beta_SM_E  *( (*Bz3D)(nx_d-2,j,k) -(*Bz_xvalmax_Long)(j,k))
                +                    Gamma_SM_E * bzE
                +                    Zeta_SM_E   *( (*Bx3D)(0,j,k+1)-(*Bx_xvalmax_Long)(j,k+1) )
                +                    Eta_SM_E *( (*Bx3D)(0,j,k)-(*Bx_xvalmax_Long)(j,k) )
                +                    (*Bz_xvalmax_Long)(j,k);
            
            }//k  ---end compute Bz
        }//j  ---end compute Bz
        
        
    }//if Eastern

}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_SM::apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isSouthern() ) {
        
        // Static cast of the fields
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
#ifdef _PATCH3D_TODO
#endif

        
    }//if Southern
    
}
// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_SM::apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isNorthern() ) {
        
        // Static cast of the fields
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
#ifdef _PATCH3D_TODO
#endif
 
        
    }//if Northern
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_SM::apply_zmin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isBottom() ) {
        
        // Static cast of the fields
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
#ifdef _PATCH3D_TODO
#endif

    } //if Bottom
    
}
// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_SM::apply_zmax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isTop() ) {
        
        // Static cast of the fields
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        
#ifdef _PATCH3D_TODO
#endif

    } //if Top
    
}


