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
    
    //! \todo (AB) Check optimal angle for buneman BC
    
    double theta  = 0.0*conv_deg2rad; //0.0;
    double phi    = 0.0*conv_deg2rad; //0.0;
    double cb =  cos(theta) * cos(phi) / (1.0 + cos(theta)*cos(phi));
    double ce =  1 - cb ;


    // Xmin boundary
    Alpha_BM_xmin    = (dt_ov_dx - 1.)  / (dt_ov_dx + 1.) ;
    Beta_BM_xmin     =  dt_ov_dy        / (dt_ov_dx + 1.) ;
    Gamma_BM_xmin    =  dt_ov_dz        / (dt_ov_dx + 1.) ;
    
    // Xmax boundary
    theta         = M_PI;
    Alpha_BM_xmax    = 2.0                      ;
    Beta_BM_xmax     = - (cos(theta)+dt_ov_dx)  ;
    Gamma_BM_xmax    = 4.0 * cos(theta)         ;
     
    // Ymin boundary
    theta  = 0.0;
    Alpha_BM_ymin    = 2.0                     ;
    Beta_BM_ymin     = - (cos(theta)-dt_ov_dy) ;
    
    // Ymax boundary
    theta  = M_PI;
    Alpha_BM_ymax    = 2.0                     ;
    Beta_BM_ymax     = - (cos(theta)+dt_ov_dy) ;
    
    // Zmin boundary
    theta  = 0.0;
    Alpha_BM_zmin    = 2.0                     ;
    Beta_BM_zmin     = - (cos(theta)-dt_ov_dz) ;
    
    // Zmax boundary
    theta         = M_PI;
    Alpha_BM_zmax    = 2.0                      ;
    Beta_BM_zmax     = - (cos(theta)+dt_ov_dz)  ;
    
}

ElectroMagnBC3D_BM::~ElectroMagnBC3D_BM()
{
    if (Bx_val) delete Bx_val ;
    if (By_val) delete By_val ;
    if (Bz_val) delete Bz_val ;
}

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
        Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
        //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
        //Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
        Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
        Field3D* Bx3D_old = static_cast<Field3D*>(EMfields->Bx_m);
        Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
        Field3D* By3D_old = static_cast<Field3D*>(EMfields->By_m);
        Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);
        Field3D* Bz3D_old = static_cast<Field3D*>(EMfields->Bz_m);
        
        
        // for By^(d,p,d)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            for (unsigned int k=0 ; k<nz_d ; k++) {
                (*By3D)(0,j,k) =  (*By3D_old)(1,j,k)
                +                 Alpha_BM_xmin    * ( (*By3D)(1,j,k  ) - (*By3D_old)(0,j,k) ) 
                -                 cb*Beta_BM_xmin  * ( (*Bx3D)(0,j+1,k) - (*Bx3D    )(0,j,k)  + (*Bx3D_old)(0,j+1,k) - (*Bx3D_old)(0,j,k)  )
                -                 ce*Gamma_BM_xmin * ( (*Ex3D)(0,j,k)   - (*Ex3D    )(0,j,k-1)+ (*Ex3D)(1,j,k)       - (*Ex3D    )(1,j,k-1)) ;
            }// k  ---end compute By
        }//j  ---end compute By
        
        // for Bz^(d,d,p)
        for (unsigned int j=0 ; j<ny_d ; j++) {
            for (unsigned int k=0 ; k<nz_p ; k++) {
                                
                (*Bz3D)(0,j,k) =  (*Bz3D_old)(1,j,k)
                +                 Alpha_BM_xmin    * ( (*Bz3D)(1,j,k  ) - (*Bz3D_old)(0,j,k) ) 
                -                 cb*Beta_BM_xmin  * ( (*Bx3D)(0,j,k+1) - (*Bx3D    )(0,j,k)  + (*Bx3D_old)(0,j,k+1) - (*Bx3D_old)(0,j,k)  )
                +                 ce*Gamma_BM_xmin * ( (*Ex3D)(0,j,k)   - (*Ex3D    )(0,j-1,k)+ (*Ex3D)(1,j,k)       - (*Ex3D    )(1,j-1,k)) ;
                
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
                
                (*By3D)(nx_d-1,j,k) = Alpha_BM_xmax   * (*Ez3D)(nx_p-1,j,k)
                +                   Beta_BM_xmax    *( (*By3D)(nx_d-2,j,k) -(*By_val)(j,k))
                +                   Gamma_BM_xmax   * byE
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
                
                (*Bz3D)(nx_d-1,j,k) = -Alpha_BM_xmax * (*Ey3D)(nx_p-1,j,k)
                +                    Beta_BM_xmax  *( (*Bz3D)(nx_d-2,j,k) -(*Bz_val)(j,k))
                +                    Gamma_BM_xmax * bzE
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
                (*Bx3D)(i,0,k) = - Alpha_BM_ymin   * (*Ez3D)(i,0,k)
                +              Beta_BM_ymin    *( (*Bx3D)(i,1,k)-(*Bx_val)(i,k))
                +              (*Bx_val)(i,k);
            }// k  ---end compute Bx
        }//i  ---end compute Bx
        
        // for Bz^(d,d,p)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int k=0 ; k<nz_p ; k++) {
                (*Bz3D)(i,0,k) = Alpha_BM_ymin   * (*Ex3D)(i,0,k)
                +              Beta_BM_ymin    *( (*Bz3D)(i,1,k)-(*Bz_val)(i,k))
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
                
                (*Bx3D)(i,ny_d-1,k) = -Alpha_BM_ymax * (*Ez3D)(i,ny_p-1,k)
                +                    Beta_BM_ymax  *( (*Bx3D)(i,ny_d-2,k) -(*Bx_val)(i,k))
                +                    (*Bx_val)(i,k);
                
            }//k  ---end compute Bz
        }//j  ---end compute Bz
        
        // for Bz^(d,d,p)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int k=0 ; k<nz_p ; k++) {
                
                (*Bz3D)(i,ny_d-1,k) = Alpha_BM_ymax   * (*Ex3D)(i,ny_p-1,k)
                +                   Beta_BM_ymax    *( (*Bz3D)(i,ny_d-2,k) -(*Bz_val)(i,k))
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
                
                (*Bx3D)(i,j,0) = Alpha_BM_zmin   * (*Ey3D)(i,j,0)
                +              Beta_BM_zmin    *( (*Bx3D)(i,j,1)-(*Bx_val)(i,j))
                +              (*Bx_val)(i,j);
            }// j  ---end compute Bx
        }//i  ---end compute Bx
        
        // for By^(d,p,d)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                
                (*By3D)(i,j,0) = - Alpha_BM_zmin   * (*Ex3D)(i,j,0)
                +              Beta_BM_zmin    *( (*By3D)(i,j,1)-(*By_val)(i,j))
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
                
                (*Bx3D)(i,j,nz_d-1) = Alpha_BM_zmax   * (*Ey3D)(i,j,nz_p-1)
                +                   Beta_BM_zmax    *( (*Bx3D)(i,j,nz_d-2) -(*Bx_val)(i,j))
                +                   (*Bx_val)(i,j);
                
            }//j  ---end compute Bx
        }//i  ---end compute Bx
        
        
        // for By^(d,p,d)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                
                (*By3D)(i,j,nz_d-1) = -Alpha_BM_zmax * (*Ex3D)(i,j,nz_p-1)
                +                    Beta_BM_zmax  *( (*By3D)(i,j,nz_d-2) -(*By_val)(i,j))
                +                    (*By_val)(i,j);
                
            }//j  ---end compute By
        }//i  ---end compute By
        
    }
}
