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

ElectroMagnBC3D_SM::ElectroMagnBC3D_SM( Params &params, Patch* patch, unsigned int _min_max )
: ElectroMagnBC( params, patch, _min_max )
{
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params.n_space[0]*params.global_factor[0]+1+2*params.oversize[0];
    nx_d = nx_p+1;
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params.n_space[1]*params.global_factor[1]+1+2*params.oversize[1];
    ny_d = ny_p+1;
    // number of nodes of the primal and dual grid in the z-direction
    nz_p = params.n_space[2]*params.global_factor[2]+1+2*params.oversize[2];
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
    double theta  = params.EM_BCs_theta[0][0];
    double factor = 1.0 / (cos(theta) + dt_ov_dx);
    Alpha_SM_W    = 2.0                     * factor;
    Beta_SM_W     = - (cos(theta)-dt_ov_dx) * factor;
    Gamma_SM_W    = 4.0 * cos(theta)        * factor;
    Delta_SM_W    = - (sin(theta)+dt_ov_dy) * factor;
    Epsilon_SM_W  = - (sin(theta)-dt_ov_dy) * factor;
    Zeta_SM_W     = - dt_ov_dz              * factor;
    Eta_SM_W      =   dt_ov_dz              * factor;
    
    // Xmax boundary
    theta  = params.EM_BCs_theta[0][1];
    factor        = 1.0 / (cos(theta) - dt_ov_dx);
    Alpha_SM_E    = 2.0                      * factor;
    Beta_SM_E     = - (cos(theta)+dt_ov_dx)  * factor;
    Gamma_SM_E    = 4.0 * cos(theta)         * factor;
    Delta_SM_E    = - (sin(theta)+dt_ov_dy)  * factor;
    Epsilon_SM_E  = - (sin(theta)-dt_ov_dy)  * factor;
    Zeta_SM_E     = - dt_ov_dz              * factor;
    Eta_SM_E      =   dt_ov_dz              * factor;
    
    // Ymin boundary
    theta  = params.EM_BCs_theta[1][0];
    factor = 1.0 / (cos(theta) + dt_ov_dy );
    Alpha_SM_S    = 2.0                     * factor;
    Beta_SM_S     = - (cos(theta)-dt_ov_dy) * factor;
    Delta_SM_S    = - (sin(theta)+dt_ov_dz) * factor;
    Epsilon_SM_S  = - (sin(theta)-dt_ov_dz) * factor;
    Zeta_SM_S     = - dt_ov_dx              * factor;
    Eta_SM_S      =   dt_ov_dx              * factor;
    
    // Ymax boundary
    theta  = params.EM_BCs_theta[1][1];
    factor = 1.0 / (cos(theta) - dt_ov_dy);
    Alpha_SM_N    = 2.0                     * factor;
    Beta_SM_N     = - (cos(theta)+dt_ov_dy) * factor;
    Delta_SM_N    = - (sin(theta)+dt_ov_dz) * factor;
    Epsilon_SM_N  = - (sin(theta)-dt_ov_dz) * factor;
    Zeta_SM_N     = - dt_ov_dx              * factor;
    Eta_SM_N      =   dt_ov_dx              * factor;
    
    // Zmin boundary
    theta  = params.EM_BCs_theta[2][0];
    factor = 1.0 / (cos(theta) + dt_ov_dz);
    Alpha_SM_B    = 2.0                     * factor;
    Beta_SM_B     = - (cos(theta)-dt_ov_dz) * factor;
    Delta_SM_B    = - (sin(theta)+dt_ov_dx) * factor;
    Epsilon_SM_B  = - (sin(theta)-dt_ov_dx) * factor;
    Zeta_SM_B     = - dt_ov_dy              * factor;
    Eta_SM_B      =   dt_ov_dy              * factor;
    
    // Zmax boundary
    theta  = params.EM_BCs_theta[2][1];
    factor        = 1.0 / (cos(theta) - dt_ov_dz);
    Alpha_SM_T    = 2.0                      * factor;
    Beta_SM_T     = - (cos(theta)+dt_ov_dz)  * factor;
    Delta_SM_T    = - (sin(theta)+dt_ov_dx)  * factor;
    Epsilon_SM_T  = - (sin(theta)-dt_ov_dx)  * factor;
    Zeta_SM_T     = - dt_ov_dy              * factor;
    Eta_SM_T      =   dt_ov_dy              * factor;
    
}

ElectroMagnBC3D_SM::~ElectroMagnBC3D_SM()
{
    if (Bx_val) delete Bx_val ;
    if (By_val) delete By_val ;
    if (Bz_val) delete Bz_val ;
}


// Magnetic field Bx^(p,d,d)
// Magnetic field By^(d,p,d)
// Magnetic field Bz^(d,d,p)

void ElectroMagnBC3D_SM::save_fields(Field* my_field, Patch* patch) {
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


void ElectroMagnBC3D_SM::disableExternalFields()
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
void ElectroMagnBC3D_SM::apply(ElectroMagn* EMfields, double time_dual, Patch* patch)
{

    // Static cast of the fields
   Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
   Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
   Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
   Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
   Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
   Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
   vector<double> pos(2);

    if ( min_max==0 && patch->isXmin() ) {
        
        // for By^(d,p,d)
        for (unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++) {
            pos[0] = patch->getDomainLocalMin(1) + (j - EMfields->oversize[1])*dy;
            for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {
                pos[1] = patch->getDomainLocalMin(2) + (k -0.5 - EMfields->oversize[2])*dz;
                // Lasers
                double byW = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    byW += vecLaser[ilaser]->getAmplitude0(pos, time_dual, j, k);
                }
                
                (*By3D)(0,j,k) = Alpha_SM_W   * (*Ez3D)(0,j,k)
                +              Beta_SM_W    *( (*By3D)(1,j,k)-(*By_val)(j,k))
                +              Gamma_SM_W   * byW
                +              Delta_SM_W   *( (*Bx3D)(0,j+1,k)-(*Bx_val)(j+1,k) )
                +              Epsilon_SM_W *( (*Bx3D)(0,j,k)-(*Bx_val)(j,k) )
                +              (*By_val)(j,k);
            }// k  ---end compute By
        }//j  ---end compute By
        
        // for Bz^(d,d,p)
         for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {
            pos[0] = patch->getDomainLocalMin(1) + (j - 0.5 - EMfields->oversize[1])*dy;
            for (unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++) {
                pos[1] = patch->getDomainLocalMin(2) + (k - EMfields->oversize[2])*dz;
                // Lasers
                double bzW = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    bzW += vecLaser[ilaser]->getAmplitude1(pos, time_dual, j, k);
                }
                
                (*Bz3D)(0,j,k) = - Alpha_SM_W   * (*Ey3D)(0,j,k)
                +              Beta_SM_W    *( (*Bz3D)(1,j,k)-(*Bz_val)(j,k))
                +              Gamma_SM_W   * bzW
                +              Zeta_SM_W   *( (*Bx3D)(0,j,k+1)-(*Bx_val)(j,k+1) )
                +              Eta_SM_W *( (*Bx3D)(0,j,k)-(*Bx_val)(j,k) )
                +              (*Bz_val)(j,k);
                
            }// k  ---end compute Bz
        }//j  ---end compute Bz

        if (patch->isZmax()){ // Xmin/Zmax
            // Bz[0,j,nz_p-1] + beta(-x)Bx[0,j,nz_p-1] = S(-x)
            EMfields->beta_edge[3] = - Zeta_SM_W;
            EMfields->S_edge[3].resize(ny_d);
            unsigned int k = nz_p - 1 ;
            pos[1] = patch->getDomainLocalMin(2) + (k - EMfields->oversize[2])*dz;
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {
                pos[0] = patch->getDomainLocalMin(1) + (j - 0.5 - EMfields->oversize[1])*dy;
                // Lasers
                double bzW = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    bzW += vecLaser[ilaser]->getAmplitude1(pos, time_dual, j, k);
                }

                EMfields->S_edge[3][j] =  - Alpha_SM_W   * (*Ey3D)(0,j,k)
                +              Beta_SM_W   *( (*Bz3D)(1,j,k) -(*Bz_val)(j,k))
                +              Gamma_SM_W  * bzW
                +              Zeta_SM_W   *(                -(*Bx_val)(j,k+1) )
                +              Eta_SM_W    *( (*Bx3D)(0,j,k) -(*Bx_val)(j,k) )
                +              (*Bz_val)(j,k);
            }
        }

        if (patch->isZmin()){ // Xmin/Zmin
            // Bz[0,j,0] + beta(-x)Bx[0,j,0] = S(-x)
            EMfields->beta_edge[2] = - Eta_SM_W;
            EMfields->S_edge[2].resize(ny_d);
            unsigned int k = 0;
            pos[1] = patch->getDomainLocalMin(2) + (k - EMfields->oversize[2])*dz;
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {
                pos[0] = patch->getDomainLocalMin(1) + (j - 0.5 - EMfields->oversize[1])*dy;
                // Lasers
                double bzW = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    bzW += vecLaser[ilaser]->getAmplitude1(pos, time_dual, j, k);
                }

                EMfields->S_edge[2][j] =  - Alpha_SM_W   * (*Ey3D)(0,j,k)
                +              Beta_SM_W    *( (*Bz3D)(1,j,k)-(*Bz_val)(j,k))
                +              Gamma_SM_W   * bzW
                +              Zeta_SM_W   *( (*Bx3D)(0,j,k+1)-(*Bx_val)(j,k+1) )
                +              Eta_SM_W *(               -(*Bx_val)(j,k) )
                +              (*Bz_val)(j,k);
            }
        }

       if (patch->isYmax()){ // Xmin/Ymax
            // By[0,ny_p-1,k] + beta(-x)Bx[0,ny_p,k] = S(-x)
            EMfields->beta_edge[1] = - Delta_SM_W;
            EMfields->S_edge[1].resize(nz_d);
            unsigned int j = ny_p - 1 ;
            pos[0] = patch->getDomainLocalMin(1) + (j - EMfields->oversize[1])*dy;
            for (unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++) {
                pos[1] = patch->getDomainLocalMin(2) + (k - 0.5 - EMfields->oversize[2])*dz;
                // Lasers
                double byW = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    byW += vecLaser[ilaser]->getAmplitude0(pos, time_dual, j, k);
                }
                EMfields->S_edge[1][k] =  Alpha_SM_W   * (*Ez3D)(0,j,k)
                +              Beta_SM_W    *( (*By3D)(1,j,k)-(*By_val)(j,k))
                +              Gamma_SM_W   * byW
                +              Delta_SM_W   *(               -(*Bx_val)(j+1,k) )
                +              Epsilon_SM_W *( (*Bx3D)(0,j,k)-(*Bx_val)(j,k) )
                +              (*By_val)(j,k);
            }
        }
 
        if (patch->isYmin()){ // Xmin/Ymin
            // By[0,0,k] + beta(-x)Bx[0,0,k] = S(-x)
            EMfields->beta_edge[0] = - Epsilon_SM_W;
            EMfields->S_edge[0].resize(nz_d);
            unsigned int j = 0;
            pos[0] = patch->getDomainLocalMin(1) + (j - EMfields->oversize[1])*dy;
            for (unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++) {
                pos[1] = patch->getDomainLocalMin(2) + (k - 0.5 - EMfields->oversize[2])*dz;
                // Lasers
                double byW = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    byW += vecLaser[ilaser]->getAmplitude0(pos, time_dual, j, k);
                }
                EMfields->S_edge[0][k] =  Alpha_SM_W   * (*Ez3D)(0,j,k)
                +              Beta_SM_W    *( (*By3D)(1,j,k)  -(*By_val)(j,k))
                +              Gamma_SM_W   * byW
                +              Delta_SM_W   *( (*Bx3D)(0,j+1,k)-(*Bx_val)(j+1,k) )
                +              Epsilon_SM_W *(                 -(*Bx_val)(j,k) )
                +              (*By_val)(j,k);
            }
        }
    }
    else if (min_max==1 && patch->isXmax() ) {
        
        // for By^(d,p,d)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            pos[0] = patch->getDomainLocalMin(1) + (j - EMfields->oversize[1])*dy;
            for (unsigned int k=0 ; k<nz_d-patch->isZmax() ; k++) {
                pos[1] = patch->getDomainLocalMin(2) + (k - 0.5 - EMfields->oversize[2])*dz;
                // Lasers
                double byE = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    byE += vecLaser[ilaser]->getAmplitude0(pos, time_dual, j, k);
                }
                
                (*By3D)(nx_d-1,j,k) = Alpha_SM_E   * (*Ez3D)(nx_p-1,j,k)
                +                   Beta_SM_E    *( (*By3D)(nx_d-2,j,k) -(*By_val)(j,k))
                +                   Gamma_SM_E   * byE
                +                   Delta_SM_E   *( (*Bx3D)(nx_p-1,j+1,k) -(*Bx_val)(j+1,k))// Check x-index
                +                   Epsilon_SM_E *( (*Bx3D)(nx_p-1,j,k) -(*Bx_val)(j,k))
                +                   (*By_val)(j,k);
                
            }//k  ---end compute By
        }//j  ---end compute By
        
        // for Bz^(d,d,p)
        for (unsigned int j=0 ; j<ny_d ; j++) {
            pos[0] = patch->getDomainLocalMin(1) + (j - 0.5 - EMfields->oversize[1])*dy;
            for (unsigned int k=0 ; k<nz_p-patch->isZmax() ; k++) {
                pos[1] = patch->getDomainLocalMin(2) + (k - EMfields->oversize[2])*dz;
                // Lasers
                double bzE = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    bzE += vecLaser[ilaser]->getAmplitude1(pos, time_dual, j, k);
                }
                
                (*Bz3D)(nx_d-1,j,k) = -Alpha_SM_E * (*Ey3D)(nx_p-1,j,k)
                +                    Beta_SM_E  *( (*Bz3D)(nx_d-2,j,k) -(*Bz_val)(j,k))
                +                    Gamma_SM_E * bzE
                +                    Zeta_SM_E   *( (*Bx3D)(nx_p-1,j,k+1)-(*Bx_val)(j,k+1) )
                +                    Eta_SM_E *( (*Bx3D)(nx_p-1,j,k)-(*Bx_val)(j,k) )
                +                    (*Bz_val)(j,k);
            }//k  ---end compute Bz
        }//j  ---end compute Bz

         if (patch->isZmax()){ // Xmax/Zmax
            // Bz[nx_p,j,nz_p-1] + beta(+x)Bx[nx_p-1,j,nz_p-1] = S(+x)
            EMfields->beta_edge[7] = - Zeta_SM_E;
            EMfields->S_edge[7].resize(ny_d);
            unsigned int k = nz_p - 1 ;
            pos[1] = patch->getDomainLocalMin(2) + (k - EMfields->oversize[2])*dz;
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {
                pos[0] = patch->getDomainLocalMin(1) + (j - 0.5 - EMfields->oversize[1])*dy;
                // Lasers
                double bzE = 0.;
                for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
                    bzE += vecLaser[ilaser]->getAmplitude1(pos, time_dual, j, k);
                }

                EMfields->S_edge[7][j] = -Alpha_SM_E * (*Ey3D)(nx_p-1,j,k)
                +                    Beta_SM_E  *( (*Bz3D)(nx_d-2,j,k)   -(*Bz_val)(j,k))
                +                    Gamma_SM_E * bzE
                +                    Zeta_SM_E   *(                      -(*Bx_val)(j,k+1) )
                +                    Eta_SM_E *( (*Bx3D)(nx_p-1,j,k)     -(*Bx_val)(j,k) )
                +                    (*Bz_val)(j,k);
             }
        }

    }
    else if (min_max==2 && patch->isYmin() ) {
        
        // for Bx^(p,d,d)
        for (unsigned int i=patch->isXmin() ; i<nx_p ; i++) {
            for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  
                (*Bx3D)(i,0,k) = - Alpha_SM_S   * (*Ez3D)(i,0,k)
                +              Beta_SM_S    *( (*Bx3D)(i,1,k)-(*Bx_val)(i,k))
                +              Zeta_SM_S   *( (*By3D)(i+1,0,k)-(*By_val)(i+1,k) )
                +              Eta_SM_S *( (*By3D)(i,0,k)-(*By_val)(i,k) )
                +              (*Bx_val)(i,k);
            }// k  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bz^(d,d,p)
        for (unsigned int i=patch->isXmin() ; i<nx_d ; i++) {
            for (unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++) {  
                (*Bz3D)(i,0,k) = Alpha_SM_S   * (*Ex3D)(i,0,k)
                +              Beta_SM_S    *( (*Bz3D)(i,1,k)-(*Bz_val)(i,k))
                +              Delta_SM_S   *( (*By3D)(i,0,k+1)-(*By_val)(i,k+1) )
                +              Epsilon_SM_S *( (*By3D)(i,0,k)-(*By_val)(i,k) )
                +              (*Bz_val)(i,k);
            }// k  ---end compute Bz
        }//i  ---end compute Bz
         if (patch->isXmin()){ // Ymin/Xmin
            // Bx[0,0,k] + beta(-y)By[0,0,k] = S(-y)
            EMfields->beta_edge[8] = - Eta_SM_S;
            EMfields->S_edge[8].resize(nz_d);
            unsigned int i = 0 ;
            for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  

                EMfields->S_edge[8][k] = - Alpha_SM_S   * (*Ez3D)(i,0,k)
                +              Beta_SM_S    *( (*Bx3D)(i,1,k) -(*Bx_val)(i,k))
                +              Zeta_SM_S   *( (*By3D)(i+1,0,k)-(*By_val)(i+1,k) )
                +              Eta_SM_S *(                    -(*By_val)(i,k) )
                +              (*Bx_val)(i,k);
            }
        }

        if (patch->isXmax()){ // Ymin/Xmax
            // Bx[nx_p-1,0,k] + beta(-y)By[nx_p,0,k] = S(-y)
            EMfields->beta_edge[9] = - Zeta_SM_S;
            EMfields->S_edge[9].resize(nz_d);
            unsigned int i = nx_p - 1 ;
            for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  

                EMfields->S_edge[9][k] = - Alpha_SM_S   * (*Ez3D)(i,0,k)
                +              Beta_SM_S    *( (*Bx3D)(i,1,k) -(*Bx_val)(i,k))
                +              Zeta_SM_S   *(                 -(*By_val)(i+1,k) )
                +              Eta_SM_S *( (*By3D)(i,0,k)     -(*By_val)(i,k) )
                +              (*Bx_val)(i,k);
            }
        }

        if (patch->isZmax()){ // Ymin/Zmax
            // Bz[i,0,nz_p-1] + beta(-y)By[i,0,nz_p] = S(-y)
            EMfields->beta_edge[11] = - Delta_SM_S;
            EMfields->S_edge[11].resize(nx_d);
            unsigned int k = nz_p - 1 ;
            for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  

                EMfields->S_edge[11][i] = Alpha_SM_S   * (*Ex3D)(i,0,k)
                +              Beta_SM_S    *( (*Bz3D)(i,1,k)  -(*Bz_val)(i,k))
                +              Delta_SM_S   *(                 -(*By_val)(i,k+1) )
                +              Epsilon_SM_S *( (*By3D)(i,0,k)  -(*By_val)(i,k) )
                +              (*Bz_val)(i,k);
            }
        }

        if (patch->isZmin()){ // Ymin/Zmin
            // Bz[i,0,0] + beta(-y)By[i,0,0] = S(-y)
            EMfields->beta_edge[10] = - Epsilon_SM_S;
            EMfields->S_edge[10].resize(nx_d);
            unsigned int k = 0;
            for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  

                EMfields->S_edge[10][i] = Alpha_SM_S   * (*Ex3D)(i,0,k)
                +              Beta_SM_S    *( (*Bz3D)(i,1,k)-(*Bz_val)(i,k))
                +              Delta_SM_S   *( (*By3D)(i,0,k+1)-(*By_val)(i,k+1) )
                +              Epsilon_SM_S *(               -(*By_val)(i,k) )
                +              (*Bz_val)(i,k);
            }
        }
    }
    else if (min_max==3 && patch->isYmax() ) {
        
        // for Bx^(p,d,d)
        for (unsigned int i=patch->isXmin() ; i<nx_p ; i++) {
            for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {
                
                (*Bx3D)(i,ny_d-1,k) = -Alpha_SM_N * (*Ez3D)(i,ny_p-1,k)
                +                    Beta_SM_N  *( (*Bx3D)(i,ny_d-2,k) -(*Bx_val)(i,k))
                +                    Zeta_SM_N   *( (*By3D)(i+1,ny_p-1,k)-(*By_val)(i+1,k) )
                +                    Eta_SM_N *( (*By3D)(i,ny_p-1,k)-(*By_val)(i,k) )
                +                    (*Bx_val)(i,k);
                
            }//k  ---end compute Bz
        }//j  ---end compute Bz
        
        // for Bz^(d,d,p)
        for (unsigned int i=patch->isXmin() ; i<nx_d ; i++) {
            for (unsigned int k=patch->isZmin() ; k<nz_p-patch->isZmax() ; k++) {
                
                (*Bz3D)(i,ny_d-1,k) = Alpha_SM_N   * (*Ex3D)(i,ny_p-1,k)
                +                   Beta_SM_N    *( (*Bz3D)(i,ny_d-2,k) -(*Bz_val)(i,k))
                +                   Delta_SM_N   *( (*By3D)(i,ny_p-1,k+1) -(*By_val)(i,k+1))
                +                   Epsilon_SM_N *( (*By3D)(i,ny_p-1,k) -(*By_val)(i,k))
                +                   (*Bz_val)(i,k);
                
            }//k  ---end compute Bz
        }//j  ---end compute Bz

        if (patch->isZmax()){ // Ymax/Zmax
            // Bz[i,ny_p,nz_p-1] + beta(+y)By[i,ny_p-1,nz_p] = S(+y)
            EMfields->beta_edge[15] = - Delta_SM_N;
            EMfields->S_edge[15].resize(nx_d);
            unsigned int k = nz_p - 1 ;
            for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {

                EMfields->S_edge[15][i] =  Alpha_SM_N   * (*Ex3D)(i,ny_p-1,k)
                +                   Beta_SM_N    *( (*Bz3D)(i,ny_d-2,k)   -(*Bz_val)(i,k))
                +                   Delta_SM_N   *(                       -(*By_val)(i,k+1))
                +                   Epsilon_SM_N *( (*By3D)(i,ny_p-1,k)   -(*By_val)(i,k))
                +                   (*Bz_val)(i,k);           
            }
        }

        if (patch->isZmin()){ // Ymax/Zmin
            // Bz[i,ny_p,0] + beta(+y)By[i,ny_p-1,0] = S(+y)
            EMfields->beta_edge[14] = - Epsilon_SM_N;
            EMfields->S_edge[14].resize(nx_d);
            unsigned int k = 0 ;
            for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {

                EMfields->S_edge[14][i] =  Alpha_SM_N   * (*Ex3D)(i,ny_p-1,k)
                +                   Beta_SM_N    *( (*Bz3D)(i,ny_d-2,k)   -(*Bz_val)(i,k))
                +                   Delta_SM_N   *( (*By3D)(i,ny_p-1,k+1) -(*By_val)(i,k+1))
                +                   Epsilon_SM_N *(                       -(*By_val)(i,k))
                +                   (*Bz_val)(i,k);
            }
        }

        if (patch->isXmin()){ // Ymax/Xmin
            // Bx[0,ny_p,k] + beta(+y)By[0,ny_p-1,k] = S(+y)
            EMfields->beta_edge[12] = - Eta_SM_N;
            EMfields->S_edge[12].resize(nz_d);
            unsigned int i = 0 ;
            for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  

                EMfields->S_edge[12][k] = -Alpha_SM_N * (*Ez3D)(i,ny_p-1,k)
                +                    Beta_SM_N  *( (*Bx3D)(i,ny_d-2,k)   -(*Bx_val)(i,k))
                +                    Zeta_SM_N   *( (*By3D)(i+1,ny_p-1,k)-(*By_val)(i+1,k) )
                +                    Eta_SM_N *(                         -(*By_val)(i,k) )
                +                    (*Bx_val)(i,k);
            }
        }


        if (patch->isXmax()){ // Ymax/Xmax
            // Bx[nx_p-1,ny_p,k] + beta(+y)By[nx_p,ny_p-1,k] = S(+y)
            EMfields->beta_edge[13] = - Zeta_SM_N;
            EMfields->S_edge[13].resize(nz_d);
            unsigned int i = nx_p - 1 ;
            for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  

                EMfields->S_edge[13][k] = -Alpha_SM_N * (*Ez3D)(i,ny_p-1,k)
                +                    Beta_SM_N  *( (*Bx3D)(i,ny_d-2,k)   -(*Bx_val)(i,k))
                +                    Zeta_SM_N   *(                      -(*By_val)(i+1,k) )
                +                    Eta_SM_N *( (*By3D)(i,ny_p-1,k)     -(*By_val)(i,k) )
                +                    (*Bx_val)(i,k);
            }
        }
    }
    else if (min_max==4 && patch->isZmin() ) {
        
        // for Bx^(p,d,d)
        for (unsigned int i=patch->isXmin() ; i<nx_p ; i++) {
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {
                
                (*Bx3D)(i,j,0) = Alpha_SM_B   * (*Ey3D)(i,j,0)
                +              Beta_SM_B    *( (*Bx3D)(i,j,1)-(*Bx_val)(i,j))
                +              Delta_SM_B   *( (*Bz3D)(i+1,j,0)-(*Bz_val)(i+1,j) )
                +              Epsilon_SM_B *( (*Bz3D)(i,j,0)-(*Bz_val)(i,j) )
                +              (*Bx_val)(i,j);
            }// j  ---end compute Bx
        }//i  ---end compute Bx
        
        // for By^(d,p,d)
        for (unsigned int i=patch->isXmin() ; i<nx_d ; i++) {
            for (unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++) {
                
                (*By3D)(i,j,0) = - Alpha_SM_B   * (*Ex3D)(i,j,0)
                +              Beta_SM_B    *( (*By3D)(i,j,1)-(*By_val)(i,j))
                +              Zeta_SM_B   *( (*Bz3D)(i,j+1,0)-(*Bz_val)(i,j+1) )
                +              Eta_SM_B *( (*Bz3D)(i,j,0)-(*Bz_val)(i,j) )
                +              (*By_val)(i,j);
                
            }// j  ---end compute By
        }//i  ---end compute By
        
        if (patch->isXmax()){ // Zmin/Xmax
            // Bx[nx_p-1,j,0] + beta(-z)Bz[nx_p,j,0] = S(-z)
            EMfields->beta_edge[17] = - Delta_SM_B;
            EMfields->S_edge[17].resize(ny_d);
            unsigned int i = nx_p - 1 ;
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {  

                EMfields->S_edge[17][j] = Alpha_SM_B   * (*Ey3D)(i,j,0)
                +              Beta_SM_B    *( (*Bx3D)(i,j,1)  -(*Bx_val)(i,j))
                +              Delta_SM_B   *(                 -(*Bz_val)(i+1,j) )
                +              Epsilon_SM_B *( (*Bz3D)(i,j,0)  -(*Bz_val)(i,j) )
                +              (*Bx_val)(i,j);
            }
        }

        if (patch->isYmax()){ // Zmin/Ymax
            // By[i,ny_p-1,0] + beta(-z)Bz[i,ny_p,0] = S(-z)
            EMfields->beta_edge[19] = - Zeta_SM_B;
            EMfields->S_edge[19].resize(nx_d);
            unsigned int j = ny_p - 1 ;
            for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  

                EMfields->S_edge[19][i] = - Alpha_SM_B   * (*Ex3D)(i,j,0)
                +              Beta_SM_B    *( (*By3D)(i,j,1) -(*By_val)(i,j))
                +              Zeta_SM_B   *(                 -(*Bz_val)(i,j+1) )
                +              Eta_SM_B    *( (*Bz3D)(i,j,0)  -(*Bz_val)(i,j) )
                +              (*By_val)(i,j);
            }
        }

        if (patch->isYmin()){ // Zmin/Ymin
            // By[i,0,0] + beta(-z)Bz[i,0,0] = S(-z)
            EMfields->beta_edge[18] = - Eta_SM_B;
            EMfields->S_edge[18].resize(nx_d);
            unsigned int j = 0;
            for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  

                EMfields->S_edge[18][i] = - Alpha_SM_B   * (*Ex3D)(i,j,0)
                +              Beta_SM_B    *( (*By3D)(i,j,1)-(*By_val)(i,j))
                +              Zeta_SM_B   *( (*Bz3D)(i,j+1,0)-(*Bz_val)(i,j+1) )
                +              Eta_SM_B *(               -(*Bz_val)(i,j) )
                +              (*By_val)(i,j);
            }
        }

        if (patch->isXmin()){ // Zmin/Xmin
            // Bx[0,j,0] + beta(-z)Bz[0,j,0] = S(-z)
            EMfields->beta_edge[16] = - Epsilon_SM_B;
            EMfields->S_edge[16].resize(ny_d);
            unsigned int i = 0;
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {  

                EMfields->S_edge[16][j] = Alpha_SM_B   * (*Ey3D)(i,j,0)
                +              Beta_SM_B    *( (*Bx3D)(i,j,1)-(*Bx_val)(i,j))
                +              Delta_SM_B   *( (*Bz3D)(i+1,j,0)-(*Bz_val)(i+1,j) )
                +              Epsilon_SM_B *(               -(*Bz_val)(i,j) )
                +              (*Bx_val)(i,j);
            }

        }
        
    }
    else if (min_max==5 && patch->isZmax() ) {
        
        // for Bx^(p,d,d)
        for (unsigned int i=patch->isXmin() ; i<nx_p-patch->isXmax() ; i++) {
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {
                
                (*Bx3D)(i,j,nz_d-1) = Alpha_SM_T   * (*Ey3D)(i,j,nz_p-1)
                +                   Beta_SM_T    *( (*Bx3D)(i,j,nz_d-2) -(*Bx_val)(i,j))
                +                   Delta_SM_T   *( (*Bz3D)(i+1,j,nz_p-1) -(*Bz_val)(i+1,j))
                +                   Epsilon_SM_T *( (*Bz3D)(i,j,nz_p-1) -(*Bz_val)(i,j))
                +                   (*Bx_val)(i,j);
                
            }//j  ---end compute Bx
        }//i  ---end compute Bx
        
        
        // for By^(d,p,d)
        for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {
            for (unsigned int j=patch->isYmin() ; j<ny_p-patch->isYmax() ; j++) {
                
                (*By3D)(i,j,nz_d-1) = -Alpha_SM_T * (*Ex3D)(i,j,nz_p-1)
                +                    Beta_SM_T  *( (*By3D)(i,j,nz_d-2) -(*By_val)(i,j))
                +                    Zeta_SM_T   *( (*Bz3D)(i,j+1,nz_p-1)-(*Bz_val)(i,j+1) )
                +                    Eta_SM_T *( (*Bz3D)(i,j,nz_p-1)-(*Bz_val)(i,j) )
                +                    (*By_val)(i,j);
                
            }//j  ---end compute By
        }//i  ---end compute By

        if (patch->isYmax()){ // Zmax/Ymax
            // By[i,ny_p-1,nz_p] + beta(+z)Bz[i,ny_p,nz_p-1] = S(+z)
            EMfields->beta_edge[23] = - Zeta_SM_T;
            EMfields->S_edge[23].resize(nx_d);
            unsigned int j = ny_p - 1;
            for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  

                EMfields->S_edge[23][i] = -Alpha_SM_T  *  (*Ex3D)(i,j,nz_p-1)
                +                      Beta_SM_T   *( (*By3D)(i,j,nz_d-2)   -(*By_val)(i,j))
                +                      Zeta_SM_T   *(                       -(*Bz_val)(i,j+1)) 
                +                      Eta_SM_T    *( (*Bz3D)(i,j,nz_p-1)   -(*Bz_val)(i,j)) 
                +                    (*By_val)(i,j);
            }
        }

        if (patch->isYmin()){ // Zmax/Ymin
            // By[i,0,nz_p] + beta(+z)Bz[i,0,nz_p-1] = S(+z)
            EMfields->beta_edge[22] = - Eta_SM_T;
            EMfields->S_edge[22].resize(nx_d);
            unsigned int j = 0;
            for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  

                EMfields->S_edge[22][i] = -Alpha_SM_T  *  (*Ex3D)(i,j,nz_p-1)
                +                      Beta_SM_T   *( (*By3D)(i,j,nz_d-2)   -(*By_val)(i,j))
                +                      Zeta_SM_T   *( (*Bz3D)(i,j+1,nz_p-1) -(*Bz_val)(i,j+1)) 
                +                      Eta_SM_T    *(                       -(*Bz_val)(i,j)) 
                +                    (*By_val)(i,j);
            }
        }

        if (patch->isXmin()){ // Zmax/Xmin
            // Bx[0,j,nz_p] + beta(+z)Bz[0,j,nz_p-1] = S(+z)
            EMfields->beta_edge[20] = - Epsilon_SM_T;
            EMfields->S_edge[20].resize(ny_d);
            unsigned int i = 0;
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {  

                EMfields->S_edge[20][j] = Alpha_SM_T   * (*Ey3D)(i,j,nz_p-1)
                +                   Beta_SM_T    *( (*Bx3D)(i,j,nz_d-2)   -(*Bx_val)(i,j))
                +                   Delta_SM_T   *( (*Bz3D)(i+1,j,nz_p-1) -(*Bz_val)(i+1,j))
                +                   Epsilon_SM_T *(                       -(*Bz_val)(i,j))
                +                   (*Bx_val)(i,j);
            }
        }

        if (patch->isXmax()){ // Zmax/Xmax
            // Bx[nx_p-1,j,nz_p] + beta(+z)Bz[nx_p,j,nz_p-1] = S(+z)
            EMfields->beta_edge[21] = - Delta_SM_T;
            EMfields->S_edge[21].resize(ny_d);
            unsigned int i = nx_p - 1;
            for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {  

                EMfields->S_edge[21][j] = Alpha_SM_T   * (*Ey3D)(i,j,nz_p-1)
                +                   Beta_SM_T    *( (*Bx3D)(i,j,nz_d-2)   -(*Bx_val)(i,j))
                +                   Delta_SM_T   *(                       -(*Bz_val)(i+1,j))
                +                   Epsilon_SM_T *( (*Bz3D)(i,j,nz_p-1)   -(*Bz_val)(i,j))
                +                   (*Bx_val)(i,j);
            }
        } 
    }

 //Now deal with the edges 
    //DO NOT SUPPORT EXTERNAL FIELDS. B*_val arrays do not have proper size.
    if (min_max == 5) {

        double one_ov_dbeta ;
        if (patch->isXmin()){

            unsigned int i = 0;
            if (patch->isYmin()){ 
                // Xmin/Ymin
                // edge 0 : By[0,0,k] + beta(-x)Bx[0,0,k] = S(-x)
                // edge 8 : Bx[0,0,k] + beta(-y)By[0,0,k] = S(-y)
                one_ov_dbeta = 1./(1. - EMfields->beta_edge[0]*EMfields->beta_edge[8]);
                for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  
                    (*By3D)(i,0,k  ) = ( EMfields->S_edge[0][k] - EMfields->beta_edge[0]* EMfields->S_edge[8][k]) * one_ov_dbeta ;
                    (*Bx3D)(i,0,k  ) =   EMfields->S_edge[8][k] - EMfields->beta_edge[8]*(*By3D)(i,0,k  ) ;
                }
            }//End Xmin Ymin edge

            if (patch->isYmax()){ 
                // Xmin/Ymax
                //edge 1 : By[0,ny_p-1,k] + beta(-x)Bx[0,ny_p,k] = S(-x)
                //edge 12 : Bx[0,ny_p,k] + beta(-y)By[0,ny_p-1,k] = S(-y)
                one_ov_dbeta = 1./(1. - EMfields->beta_edge[1]*EMfields->beta_edge[12]);
                for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  
                    (*By3D)(i,ny_p-1,k  ) = ( EMfields->S_edge[1][k] - EMfields->beta_edge[1]* EMfields->S_edge[12][k]) * one_ov_dbeta ;
                    (*Bx3D)(i,ny_p,k  ) =   EMfields->S_edge[12][k] - EMfields->beta_edge[12]*(*By3D)(i,ny_p-1,k  ) ;
                }
            }// End Xmin Ymax edge

            if (patch->isZmin()){ 
                // Xmin/Zmin
                // edge 2  : Bz[0,j,0] + beta(-x)Bx[0,j,0] = S(-x)
                // edge 16 : Bx[0,j,0] + beta(-z)Bz[0,j,0] = S(-z)
                double one_ov_dbeta = 1./(1. - EMfields->beta_edge[2]*EMfields->beta_edge[16]);
                for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {  
                    (*Bz3D)(i,j,0  ) = ( EMfields->S_edge[2][j] - EMfields->beta_edge[2]* EMfields->S_edge[16][j]) * one_ov_dbeta ;
                    (*Bx3D)(i,j,0  ) =   EMfields->S_edge[16][j] - EMfields->beta_edge[16]*(*Bz3D)(i,j,0  ) ;
                }
            } // End Xmin Zmin edge

            if (patch->isZmax()){ 
                // Xmin/Zmax
                // edge 3  : Bz[0,j,nz_p-1] + beta(-x)Bx[0,j,nz_p] = S(-x)
                // edge 20 : Bx[0,j,nz_p] + beta(+z)Bz[0,j,nz_p-1] = S(+z)
                one_ov_dbeta = 1./(1. - EMfields->beta_edge[3]*EMfields->beta_edge[20]);
                for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {  
                    (*Bz3D)(i,j,nz_p-1  ) = ( EMfields->S_edge[3][j] - EMfields->beta_edge[3]* EMfields->S_edge[20][j]) * one_ov_dbeta ;
                    (*Bx3D)(i,j,nz_p  ) =   EMfields->S_edge[20][j] - EMfields->beta_edge[20]*(*Bz3D)(i,j,nz_p-1  ) ;
                }
            }//End Xmin/Zmax edge
        } //End series of Xmin edges

        if (patch->isXmax()){

            unsigned int i = nx_p - 1;
            //if (patch->isYmin()){ 
            //    // Xmax/Ymin
            //    // edge 4 : By[nx_p,0,k] + beta(+x)Bx[nx_p-1,0,k] = S(+x)
            //    // edge 9 : Bx[nx_p-1,0,k] + beta(-y)By[nx_p-1,0,k] = S(-y)
            //    one_ov_dbeta = 1./(1. - EMfields->beta_edge[4]*EMfields->beta_edge[9]);
            //    for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  
            //        (*By3D)(i+1,0,k  ) = ( EMfields->S_edge[4][k] - EMfields->beta_edge[4]* EMfields->S_edge[9][k]) * one_ov_dbeta ;
            //        (*Bx3D)(i,0,k  ) =   EMfields->S_edge[9][k] - EMfields->beta_edge[9]*(*By3D)(i+1,0,k  ) ;
            //    }
            //}//End Xmax Ymin edge

            //if (patch->isYmax()){ 
            //    // Xmax/Ymax
            //    //edge 5 :  By[nx_p,ny_p-1,k] + beta(+x)Bx[nx_p-1,ny_p,k] = S(+x)
            //    //edge 13 : Bx[nx_p-1,ny_p,k] + beta(-y)By[nx_p,ny_p-1,k] = S(-y)
            //    one_ov_dbeta = 1./(1. - EMfields->beta_edge[5]*EMfields->beta_edge[13]);
            //    for (unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++) {  
            //        (*By3D)(i+1,ny_p-1,k  ) = ( EMfields->S_edge[5][k] - EMfields->beta_edge[5]* EMfields->S_edge[13][k]) * one_ov_dbeta ;
            //        (*Bx3D)(i,ny_p,k  ) =   EMfields->S_edge[13][k] - EMfields->beta_edge[13]*(*By3D)(i+1,ny_p-1,k  ) ;
            //    }
            //}// End Xmax Ymax edge

            //if (patch->isZmin()){ 
            //    // Xmax/Zmin
            //    // edge 6  : Bz[nx_p,j,0] + beta(+x)Bx[nx_p-1,j,0] = S(+x)
            //    // edge 17 : Bx[nx_p-1,j,0] + beta(-z)Bz[nx_p,j,0] = S(-z)
            //    double one_ov_dbeta = 1./(1. - EMfields->beta_edge[6]*EMfields->beta_edge[17]);
            //    for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {  
            //        (*Bz3D)(i+1,j,0  ) = ( EMfields->S_edge[6][j] - EMfields->beta_edge[6]* EMfields->S_edge[17][j]) * one_ov_dbeta ;
            //        (*Bx3D)(i,j,0  ) =   EMfields->S_edge[17][j] - EMfields->beta_edge[17]*(*Bz3D)(i+1,j,0  ) ;
            //    }
            //} // End Xmax Zmin edge

            if (patch->isZmax()){ 
                // Xmax/Zmax
                // edge 7  : Bz[nx_p,j,nz_p-1] + beta(+x)Bx[nx_p-1,j,nz_p] = S(+x)
                // edge 21 : Bx[nx_p-1,j,nz_p] + beta(+z)Bz[nx_p,j,nz_p-1] = S(+z)
                one_ov_dbeta = 1./(1. - EMfields->beta_edge[7]*EMfields->beta_edge[21]);
                for (unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++) {  
                    (*Bz3D)(i+1,j,nz_p-1  ) = ( EMfields->S_edge[7][j] - EMfields->beta_edge[7]* EMfields->S_edge[21][j]) * one_ov_dbeta ;
                    (*Bx3D)(i,j,nz_p  ) =   EMfields->S_edge[21][j] - EMfields->beta_edge[21]*(*Bz3D)(i+1,j,nz_p-1  ) ;
                }
            }//End Xmax/Zmax edge
        } //End series of Xmax edges

        if (patch->isYmin()){ 
            unsigned int j = 0;
            if (patch->isZmin()){  //Ymin/Zmin
                // edge 10 : Bz[i,0,0] + beta(-y)By[i,0,0] = S(-y)
                // edge 18 : By[i,0,0] + beta(-z)Bz[i,0,0] = S(-z)
                double one_ov_dbeta = 1./(1. - EMfields->beta_edge[10]*EMfields->beta_edge[18]);
                for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  
                    (*Bz3D)(i,j,0  ) = ( EMfields->S_edge[10][i] - EMfields->beta_edge[10]* EMfields->S_edge[18][i]) * one_ov_dbeta ;
                    (*By3D)(i,j,0  ) =   EMfields->S_edge[18][i] - EMfields->beta_edge[18]*(*Bz3D)(i,j,0  ) ;
                }
            }

            if (patch->isZmax()){ 
                // Ymin/Zmax
                //edge 11 : Bz[i,0,nz_p-1] + beta(-y)By[i,0,nz_p]   = S(-y)
                //edge 22 : By[i,0,nz_p]   + beta(+z)Bz[i,0,nz_p-1] = S(+z)
                one_ov_dbeta = 1./(1. - EMfields->beta_edge[11]*EMfields->beta_edge[22]);
                for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  
                    (*Bz3D)(i,j,nz_p-1  ) = ( EMfields->S_edge[11][i] - EMfields->beta_edge[11]* EMfields->S_edge[22][i]) * one_ov_dbeta ;
                    (*By3D)(i,j,nz_p  ) =   EMfields->S_edge[22][i] - EMfields->beta_edge[22]*(*Bz3D)(i,j,nz_p-1  ) ;

                }
            } //End Ymin /Zmax edge
        } //End series of Ymin edges

        if (patch->isYmax()){ 
            unsigned int j = ny_p - 1;
            if (patch->isZmin()){  //Ymax/Zmin
                // edge 14 : Bz[i,ny_p,0] + beta(+y)By[i,ny_p-1,0] = S(+y)
                // edge 19 : By[i,ny_p-1,0] + beta(-z)Bz[i,ny_p,0] = S(-z)
                double one_ov_dbeta = 1./(1. - EMfields->beta_edge[14]*EMfields->beta_edge[19]);
                for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  
                    (*Bz3D)(i,j+1,0  ) = ( EMfields->S_edge[14][i] - EMfields->beta_edge[14]* EMfields->S_edge[19][i]) * one_ov_dbeta ;
                    (*By3D)(i,j,0  )   =   EMfields->S_edge[19][i] - EMfields->beta_edge[19]*(*Bz3D)(i,j+1,0  ) ;
                }
            }//End Ymax /Zmin edge

            if (patch->isZmax()){ 
                // Ymax/Zmax
                //edge 15 : Bz[i,ny_p,nz_p-1] + beta(+y)By[i,ny_p-1,nz_p]   = S(+y)
                //edge 23 : By[i,ny_p-1,nz_p]   + beta(+z)Bz[i,ny_p,nz_p-1] = S(+z)
                one_ov_dbeta = 1./(1. - EMfields->beta_edge[15]*EMfields->beta_edge[23]);
                for (unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++) {  
                    (*Bz3D)(i,j+1,nz_p-1  ) = ( EMfields->S_edge[15][i] - EMfields->beta_edge[15]* EMfields->S_edge[23][i]) * one_ov_dbeta ;
                    (*By3D)(i,j,nz_p  ) =   EMfields->S_edge[23][i] - EMfields->beta_edge[23]*(*Bz3D)(i,j+1,nz_p-1  ) ;

                }
            } //End Ymax /Zmax edge
        } //End series of Ymax edges

    }
}
