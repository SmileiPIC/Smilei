#include "ElectroMagnBC3D_BM.h"

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

ElectroMagnBC3D_BM::ElectroMagnBC3D_BM( Params &params, Patch* patch, unsigned int _min_max )
: ElectroMagnBC( params, patch, _min_max )
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
    
    Bx_val = By_val = Bz_val = nullptr;
    
    if (min_max==0 && patch->isXmin() ) {
        // BCs at the x-border min
        dims = { ny_d, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D(dims); Bx_val->put_to(0.);
        dims = { ny_p, nz_d }; // By^(d,p,d)
        By_val = new Field2D(dims); By_val->put_to(0.);
        dims = { ny_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D(dims); Bz_val->put_to(0.);
    }
    else if (min_max==1 && patch->isXmax() ) {
        // BCs at the x-border max
        dims = { ny_d, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D(dims); Bx_val->put_to(0.);
        dims = { ny_p, nz_d }; // By^(d,p,d)
        By_val = new Field2D(dims); By_val->put_to(0.);
        dims = { ny_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D(dims); Bz_val->put_to(0.);
    }
    else if (min_max==2 && patch->isYmin() ) {
        // BCs in the y-border min
        dims = { nx_p, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D(dims); Bx_val->put_to(0.);
        dims = { nx_d, nz_d }; // By^(d,p,d)
        By_val = new Field2D(dims); By_val->put_to(0.);
        dims = { nx_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D(dims); Bz_val->put_to(0.);
    }
    else if (min_max==3 && patch->isYmax() ) {
        // BCs in the y-border mix
        dims = { nx_p, nz_d }; // Bx^(p,d,d)
        Bx_val = new Field2D(dims); Bx_val->put_to(0.);
        dims = { nx_d, nz_d }; // By^(d,p,d)
        By_val = new Field2D(dims); By_val->put_to(0.);
        dims = { nx_d, nz_p }; // Bz^(d,d,p)
        Bz_val = new Field2D(dims); Bz_val->put_to(0.);
    }
    else if (min_max==4 && patch->isZmin() ) {
        // BCs in the z-border min
        dims = { nx_p, ny_d }; // Bx^(p,d,d)
        Bx_val = new Field2D(dims); Bx_val->put_to(0.);
        dims = { nx_d, ny_p }; // By^(d,p,d)
        By_val = new Field2D(dims); By_val->put_to(0.);
        dims = { nx_d, ny_d }; // Bz^(d,d,p)
        Bz_val = new Field2D(dims); Bz_val->put_to(0.);
    }
    else if (min_max==5 && patch->isZmax() ) {
        // BCs in the z-border max
        dims = { nx_p, ny_d }; // Bx^(p,d,d)
        Bx_val = new Field2D(dims); Bx_val->put_to(0.);
        dims = { nx_d, ny_p }; // By^(d,p,d)
        By_val = new Field2D(dims); By_val->put_to(0.);
        dims = { nx_d, ny_d }; // Bz^(d,d,p)
        Bz_val = new Field2D(dims); Bz_val->put_to(0.);
    }
    
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    //! \todo (MG) Check optimal angle for Silver-Muller BCs
    
    // Xmin boundary
    double theta  = 0.0*conv_deg2rad; //0.0;
    double factor = 1.0 / (cos(theta) + dt_ov_dx);
    Alpha_BM_W    = 2.0                     * factor;
    Beta_BM_W     = - (cos(theta)-dt_ov_dx) * factor;
    Gamma_BM_W    = 4.0 * cos(theta)        * factor;
    Delta_BM_W    = - (sin(theta)+dt_ov_dy) * factor;
    Epsilon_BM_W  = - (sin(theta)-dt_ov_dy) * factor;
    Zeta_BM_W     = - dt_ov_dz              * factor;
    Eta_BM_W      =   dt_ov_dz              * factor;
    
    // Xmax boundary
    theta         = M_PI;
    factor        = 1.0 / (cos(theta) - dt_ov_dx);
    Alpha_BM_E    = 2.0                      * factor;
    Beta_BM_E     = - (cos(theta)+dt_ov_dx)  * factor;
    Gamma_BM_E    = 4.0 * cos(theta)         * factor;
    Delta_BM_E    = - (sin(theta)+dt_ov_dy)  * factor;
    Epsilon_BM_E  = - (sin(theta)-dt_ov_dy)  * factor;
    Zeta_BM_E     = - dt_ov_dz              * factor;
    Eta_BM_E      =   dt_ov_dz              * factor;
    
    // Ymin boundary
    theta  = 0.0;
    factor = 1.0 / (cos(theta) + dt_ov_dy );
    Alpha_BM_S    = 2.0                     * factor;
    Beta_BM_S     = - (cos(theta)-dt_ov_dy) * factor;
    Delta_BM_S    = - (sin(theta)+dt_ov_dz) * factor;
    Epsilon_BM_S  = - (sin(theta)-dt_ov_dz) * factor;
    Zeta_BM_S     = - dt_ov_dx              * factor;
    Eta_BM_S      =   dt_ov_dx              * factor;
    
    // Ymax boundary
    theta  = M_PI;
    factor = 1.0 / (cos(theta) - dt_ov_dy);
    Alpha_BM_N    = 2.0                     * factor;
    Beta_BM_N     = - (cos(theta)+dt_ov_dy) * factor;
    Delta_BM_N    = - (sin(theta)+dt_ov_dz) * factor;
    Epsilon_BM_N  = - (sin(theta)-dt_ov_dz) * factor;
    Zeta_BM_N     = - dt_ov_dx              * factor;
    Eta_BM_N      =   dt_ov_dx              * factor;
    
    // Zmin boundary
    theta  = 0.0;
    factor = 1.0 / (cos(theta) + dt_ov_dz);
    Alpha_BM_B    = 2.0                     * factor;
    Beta_BM_B     = - (cos(theta)-dt_ov_dz) * factor;
    Delta_BM_B    = - (sin(theta)+dt_ov_dx) * factor;
    Epsilon_BM_B  = - (sin(theta)-dt_ov_dx) * factor;
    Zeta_BM_B     = - dt_ov_dy              * factor;
    Eta_BM_B      =   dt_ov_dy              * factor;
    
    // Zmax boundary
    theta         = M_PI;
    factor        = 1.0 / (cos(theta) - dt_ov_dz);
    Alpha_BM_T    = 2.0                      * factor;
    Beta_BM_T     = - (cos(theta)+dt_ov_dz)  * factor;
    Delta_BM_T    = - (sin(theta)+dt_ov_dx)  * factor;
    Epsilon_BM_T  = - (sin(theta)-dt_ov_dx)  * factor;
    Zeta_BM_T     = - dt_ov_dy              * factor;
    Eta_BM_T      =   dt_ov_dy              * factor;
    
}

ElectroMagnBC3D_BM::~ElectroMagnBC3D_BM()
{
    if (Bx_val) delete Bx_val ;
    if (By_val) delete By_val ;
    if (Bz_val) delete Bz_val ;
}


// Magnetic field Bx^(p,d,d)
// Magnetic field By^(d,p,d)
// Magnetic field Bz^(d,d,p)

void ElectroMagnBC3D_BM::save_fields(Field* my_field, Patch* patch) {
    Field3D* field3D=static_cast<Field3D*>(my_field);
    
    if (min_max==0 && patch->isXmin() ) {
        
        if (field3D->name=="Bx"){
            field3D->extract_slice_yz(0,      Bx_val);
        }
        else if (field3D->name=="By"){
            field3D->extract_slice_yz(0,      By_val);
        }
        else if (field3D->name=="Bz"){
            field3D->extract_slice_yz(0,      Bz_val);
        }
    }
    else if (min_max==1 && patch->isXmax() ) {
        if (field3D->name=="Bx"){
            field3D->extract_slice_yz(0,      Bx_val);
        }
        else if (field3D->name=="By"){
            field3D->extract_slice_yz(0,      By_val);
        }
        else if (field3D->name=="Bz"){
            field3D->extract_slice_yz(0,      Bz_val);
        }
    }
    else if (min_max==2 && patch->isYmin() ) {
        if (field3D->name=="Bx"){
            field3D->extract_slice_xz(0,      Bx_val);
        }
        else if (field3D->name=="By"){
            field3D->extract_slice_xz(0,      By_val);
        }
        else if (field3D->name=="Bz"){
            field3D->extract_slice_xz(0,      Bz_val);
        }
    }
    else if (min_max==3 && patch->isYmax() ) {
        if (field3D->name=="Bx"){
            field3D->extract_slice_xz(ny_d-1, Bx_val);
        }
        else if (field3D->name=="By"){
            field3D->extract_slice_xz(ny_p-1, By_val);
        }
        else if (field3D->name=="Bz"){
            field3D->extract_slice_xz(ny_d-1, Bz_val);
        }
    }
    else if (min_max==4 && patch->isZmin() ) {
        
        if (field3D->name=="Bx"){
            field3D->extract_slice_xy(0,      Bx_val);
        }
        else if (field3D->name=="By"){
            field3D->extract_slice_xy(0,      By_val);
        }
        else if (field3D->name=="Bz"){
            field3D->extract_slice_xy(0,      Bz_val);
        }
    }
    else if (min_max==5 && patch->isZmax() ) {
        
        if (field3D->name=="Bx"){
            field3D->extract_slice_xy(nz_d-1, Bx_val);
        }
        else if (field3D->name=="By"){
            field3D->extract_slice_xy(nz_d-1, By_val);
        }
        else if (field3D->name=="Bz"){
            field3D->extract_slice_xy(nz_p-1, Bz_val);
        }
    }
}


void ElectroMagnBC3D_BM::disableExternalFields()
{
    delete Bx_val;
    Bx_val = NULL;
    delete By_val;
    By_val = NULL;
    delete Bz_val;
    Bz_val = NULL;
}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_BM::apply(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( min_max==0 && patch->isXmin() ) {
        
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
                
                (*By3D)(0,j,k) = Alpha_BM_W   * (*Ez3D)(0,j,k)
                +              Beta_BM_W    *( (*By3D)(1,j,k)-(*By_val)(j,k))
                +              Gamma_BM_W   * byW
                +              Delta_BM_W   *( (*Bx3D)(0,j+1,k)-(*Bx_val)(j+1,k) )
                +              Epsilon_BM_W *( (*Bx3D)(0,j,k)-(*Bx_val)(j,k) )
                +              (*By_val)(j,k);
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
                    bzW += vecLaser[ilaser]->getAmplitude1(pos, time_dual, j, k);
                }
                
                (*Bz3D)(0,j,k) = - Alpha_BM_W   * (*Ey3D)(0,j,k)
                +              Beta_BM_W    *( (*Bz3D)(1,j,k)-(*Bz_val)(j,k))
                +              Gamma_BM_W   * bzW
                +              Zeta_BM_W   *( (*Bx3D)(0,j,k+1)-(*Bx_val)(j,k+1) )
                +              Eta_BM_W *( (*Bx3D)(0,j,k)-(*Bx_val)(j,k) )
                +              (*Bz_val)(j,k);
                
            }// k  ---end compute Bz
        }//j  ---end compute Bz
    }
    else if (min_max==1 && patch->isXmax() ) {
        
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
                    byE += vecLaser[ilaser]->getAmplitude0(pos, time_dual, j, k);
                }
                
                (*By3D)(nx_d-1,j,k) = Alpha_BM_E   * (*Ez3D)(nx_p-1,j,k)
                +                   Beta_BM_E    *( (*By3D)(nx_d-2,j,k) -(*By_val)(j,k))
                +                   Gamma_BM_E   * byE
                +                   Delta_BM_E   *( (*Bx3D)(nx_p-1,j+1,k) -(*Bx_val)(j+1,k))// Check x-index
                +                   Epsilon_BM_E *( (*Bx3D)(nx_p-1,j,k) -(*Bx_val)(j,k))
                +                   (*By_val)(j,k);
                
            }//k  ---end compute By
        }//j  ---end compute By
        
        // for Bz^(d,d,p)
        pos[0] = patch->getDomainLocalMin(1) - (0.5+EMfields->oversize[1])*dy;
        for (unsigned int j=0 ; j<ny_d ; j++) {
            pos[1] = patch->getDomainLocalMin(2) - EMfields->oversize[2]*dz;
            pos[0] += dy;
            for (unsigned int k=0 ; k<nz_p ; k++) {
                pos[1] += dz;
                // Lasers
                double bzE = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    bzE += vecLaser[ilaser]->getAmplitude1(pos, time_dual, j, k);
                }
                
                (*Bz3D)(nx_d-1,j,k) = -Alpha_BM_E * (*Ey3D)(nx_p-1,j,k)
                +                    Beta_BM_E  *( (*Bz3D)(nx_d-2,j,k) -(*Bz_val)(j,k))
                +                    Gamma_BM_E * bzE
                +                    Zeta_BM_E   *( (*Bx3D)(nx_p-1,j,k+1)-(*Bx_val)(j,k+1) )
                +                    Eta_BM_E *( (*Bx3D)(nx_p-1,j,k)-(*Bx_val)(j,k) )
                +                    (*Bz_val)(j,k);
            }//k  ---end compute Bz
        }//j  ---end compute Bz
    }
    else if (min_max==2 && patch->isYmin() ) {
        
        // Static cast of the fields
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        // for Bx^(p,d,d)
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int k=0 ; k<nz_d ; k++) {
                (*Bx3D)(i,0,k) = - Alpha_BM_S   * (*Ez3D)(i,0,k)
                +              Beta_BM_S    *( (*Bx3D)(i,1,k)-(*Bx_val)(i,k))
                +              Zeta_BM_S   *( (*By3D)(i+1,0,k)-(*By_val)(i+1,k) )
                +              Eta_BM_S *( (*By3D)(i,0,k)-(*By_val)(i,k) )
                +              (*Bx_val)(i,k);
            }// k  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bz^(d,d,p)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int k=0 ; k<nz_p ; k++) {
                (*Bz3D)(i,0,k) = Alpha_BM_S   * (*Ex3D)(i,0,k)
                +              Beta_BM_S    *( (*Bz3D)(i,1,k)-(*Bz_val)(i,k))
                +              Delta_BM_S   *( (*By3D)(i,0,k+1)-(*By_val)(i,k+1) )
                +              Epsilon_BM_S *( (*By3D)(i,0,k)-(*By_val)(i,k) )
                +              (*Bz_val)(i,k);
            }// k  ---end compute Bz
        }//i  ---end compute Bz
        
    }
    else if (min_max==3 && patch->isYmax() ) {
        
        // Static cast of the fields
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        // for Bx^(p,d,d)
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int k=0 ; k<nz_d ; k++) {
                
                (*Bx3D)(i,ny_d-1,k) = -Alpha_BM_N * (*Ez3D)(i,ny_p-1,k)
                +                    Beta_BM_N  *( (*Bx3D)(i,ny_d-2,k) -(*Bx_val)(i,k))
                +                    Zeta_BM_N   *( (*By3D)(i+1,ny_p-1,k)-(*By_val)(i+1,k) )
                +                    Eta_BM_N *( (*By3D)(i,ny_p-1,k)-(*By_val)(i,k) )
                +                    (*Bx_val)(i,k);
                
            }//k  ---end compute Bz
        }//j  ---end compute Bz
        
        // for Bz^(d,d,p)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int k=0 ; k<nz_p ; k++) {
                
                (*Bz3D)(i,ny_d-1,k) = Alpha_BM_N   * (*Ex3D)(i,ny_p-1,k)
                +                   Beta_BM_N    *( (*Bz3D)(i,ny_d-2,k) -(*Bz_val)(i,k))
                +                   Delta_BM_N   *( (*By3D)(i,ny_p-1,k+1) -(*By_val)(i,k+1))
                +                   Epsilon_BM_N *( (*By3D)(i,ny_p-1,k) -(*By_val)(i,k))
                +                   (*Bz_val)(i,k);
                
            }//k  ---end compute Bz
        }//j  ---end compute Bz
    }
    else if (min_max==4 && patch->isZmin() ) {
        
        // Static cast of the fields
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        //Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        // for Bx^(p,d,d)
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                
                (*Bx3D)(i,j,0) = Alpha_BM_B   * (*Ey3D)(i,j,0)
                +              Beta_BM_B    *( (*Bx3D)(i,j,1)-(*Bx_val)(i,j))
                +              Delta_BM_B   *( (*Bz3D)(i+1,j,0)-(*Bz_val)(i+1,j) )
                +              Epsilon_BM_B *( (*Bz3D)(i,j,0)-(*Bz_val)(i,j) )
                +              (*Bx_val)(i,j);
            }// j  ---end compute Bx
        }//i  ---end compute Bx
        
        // for By^(d,p,d)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                
                (*By3D)(i,j,0) = - Alpha_BM_B   * (*Ex3D)(i,j,0)
                +              Beta_BM_B    *( (*By3D)(i,j,1)-(*By_val)(i,j))
                +              Zeta_BM_B   *( (*Bz3D)(i,j+1,0)-(*Bz_val)(i,j+1) )
                +              Eta_BM_B *( (*Bz3D)(i,j,0)-(*Bz_val)(i,j) )
                +              (*By_val)(i,j);
                
            }// j  ---end compute By
        }//i  ---end compute By
        
    }
    else if (min_max==5 && patch->isZmax() ) {
        
        // Static cast of the fields
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        //Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        
        // for Bx^(p,d,d)
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                
                (*Bx3D)(i,j,nz_d-1) = Alpha_BM_T   * (*Ey3D)(i,j,nz_p-1)
                +                   Beta_BM_T    *( (*Bx3D)(i,j,nz_d-2) -(*Bx_val)(i,j))
                +                   Delta_BM_T   *( (*Bz3D)(i+1,j,nz_p-1) -(*Bz_val)(i+1,j))
                +                   Epsilon_BM_T *( (*Bz3D)(i,j,nz_p-1) -(*Bz_val)(i,j))
                +                   (*Bx_val)(i,j);
                
            }//j  ---end compute Bx
        }//i  ---end compute Bx
        
        
        // for By^(d,p,d)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                
                (*By3D)(i,j,nz_d-1) = -Alpha_BM_T * (*Ex3D)(i,j,nz_p-1)
                +                    Beta_BM_T  *( (*By3D)(i,j,nz_d-2) -(*By_val)(i,j))
                +                    Zeta_BM_T   *( (*Bz3D)(i,j+1,nz_p-1)-(*Bz_val)(i,j+1) )
                +                    Eta_BM_T *( (*Bz3D)(i,j,nz_p-1)-(*Bz_val)(i,j) )
                +                    (*By_val)(i,j);
                
            }//j  ---end compute By
        }//i  ---end compute By
        
    }
}
