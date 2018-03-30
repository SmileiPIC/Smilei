#include "Interpolator1D3Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "Patch.h"

using namespace std;

Interpolator1D3Order::Interpolator1D3Order(Params &params, Patch* patch) : Interpolator1D(params, patch)
{
    dx_inv_   = 1.0/params.cell_length[0];
    
    dble_1ov6 = 1.0/6.0;
    dble_2ov3 = 2.0/3.0;
}

/***********************************************************************
    Interpolate the field fx defined on the primal grid
    with size nstp_x and space step stp_x_inv at the position
    xj and return the value fxj
***********************************************************************/
void Interpolator1D3Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc)
{
    double xjn, xi2, xi3;
    
    //!\todo Julien, can you check that this is indeed the centered B-field which is passed to the pusher?
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);
    
    
    // Calculate the normalized positions
    // ----------------------------------
    xjn    = particles.position(0, ipart)*dx_inv_;
    
    // Dual Grid : Ex, By, Bz
    // ----------------------
    id_   = (int)(xjn+0.50);        // position of the 2nd node
    xi  = xjn - (double)id_ + 0.5;    // normalized distance to the central node
    xi2 = xi*xi;
    xi3 = xi*xi2;
    
    // 3rd order interpolation on 4 nodes
    coeffd_[0] = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
    coeffd_[1]  = dble_2ov3 - xi2 + 0.5*xi3;
    coeffd_[2]  = dble_1ov6 + 0.5*(xi+xi2-xi3);
    coeffd_[3]  = xi3*dble_1ov6;
    
    id_ -= index_domain_begin;
    
    *(ELoc+0*nparts) = compute(coeffd_, Ex1D,   id_);  
    *(BLoc+1*nparts) = compute(coeffd_, By1D_m, id_);  
    *(BLoc+2*nparts) = compute(coeffd_, Bz1D_m, id_);  
    
    // Primal Grid : Ey, Ez, Bx
    // ------------------------
    ip_ = (int)xjn;            // index of the 2nd node
    xi  = xjn - (double)ip_;   //normalized distance to the 2nd node
    xi2 = xi*xi;
    xi3 = xi2*xi;
    
    // 3rd order interpolation on 4 nodes
    coeffp_[0]  = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
    coeffp_[1]  = dble_2ov3 - xi2 + 0.5*xi3;
    coeffp_[2]  = dble_1ov6 + 0.5*(xi+xi2-xi3);
    coeffp_[3]  = xi3*dble_1ov6;
    
    ip_ -= index_domain_begin;
    
    *(ELoc+1*nparts) = compute(coeffp_, Ey1D,   ip_);  
    *(ELoc+2*nparts) = compute(coeffp_, Ez1D,   ip_);  
    *(BLoc+0*nparts) = compute(coeffp_, Bx1D_m, ip_);
    
    
}//END Interpolator1D3Order


void Interpolator1D3Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc)
{
    
    int ipart = *istart;

    double *ELoc = &(smpi->dynamics_Epart[ithread][ipart]);
    double *BLoc = &(smpi->dynamics_Bpart[ithread][ipart]);

    // Interpolate E, B
    // Compute coefficient for ipart position    (*this)(EMfields, particles, ipart, ELoc, BLoc);
    double xjn, xi2, xi3;
    
    //!\todo Julien, can you check that this is indeed the centered B-field which is passed to the pusher?
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);
    
    
    // Calculate the normalized positions
    // ----------------------------------
    xjn    = particles.position(0, ipart)*dx_inv_;
    
    // Dual Grid : Ex, By, Bz
    // ----------------------
    id_   = (int)(xjn+0.50);        // position of the 2nd node
    xi  = xjn - (double)id_ + 0.5;    // normalized distance to the central node
    xi2 = xi*xi;
    xi3 = xi*xi2;
    
    // 3rd order interpolation on 4 nodes
    coeffd_[0] = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
    coeffd_[1]  = dble_2ov3 - xi2 + 0.5*xi3;
    coeffd_[2]  = dble_1ov6 + 0.5*(xi+xi2-xi3);
    coeffd_[3]  = xi3*dble_1ov6;
    
    id_ -= index_domain_begin;
    
    int nparts( particles.size() );
    
    *(ELoc+0*nparts) = compute(coeffd_, Ex1D,   id_);  
    *(BLoc+1*nparts) = compute(coeffd_, By1D_m, id_);  
    *(BLoc+2*nparts) = compute(coeffd_, Bz1D_m, id_);  
    
    // Primal Grid : Ey, Ez, Bx
    // ------------------------
    ip_ = (int)xjn;            // index of the 2nd node
    xi  = xjn - (double)ip_;   //normalized distance to the 2nd node
    xi2 = xi*xi;
    xi3 = xi2*xi;
    
    // 3rd order interpolation on 4 nodes
    coeffp_[0]  = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
    coeffp_[1]  = dble_2ov3 - xi2 + 0.5*xi3;
    coeffp_[2]  = dble_1ov6 + 0.5*(xi+xi2-xi3);
    coeffp_[3]  = xi3*dble_1ov6;
    
    ip_ -= index_domain_begin;
    
    *(ELoc+1*nparts) = compute(coeffp_, Ey1D,   ip_);  
    *(ELoc+2*nparts) = compute(coeffp_, Ez1D,   ip_);  
    *(BLoc+0*nparts) = compute(coeffp_, Bx1D_m, ip_);
    
    //!\todo Julien, can you check that this is indeed the centered B-field which is passed to the pusher?
    Field1D* Jx1D     = static_cast<Field1D*>(EMfields->Jx_);
    Field1D* Jy1D     = static_cast<Field1D*>(EMfields->Jy_);
    Field1D* Jz1D     = static_cast<Field1D*>(EMfields->Jz_);
    Field1D* Rho1D    = static_cast<Field1D*>(EMfields->rho_);
    
    // Primal Grid : Jy, Jz, Rho
    // ------------------------
    (*JLoc).y = compute(coeffp_, Jy1D,  ip_);  
    (*JLoc).z = compute(coeffp_, Jz1D,  ip_);  
    (*RhoLoc) = compute(coeffp_, Rho1D, ip_);    
    
    // Dual Grid : Jx
    // ----------------------
    (*JLoc).x = compute(coeffd_, Jx1D,  id_);  
    
}
void Interpolator1D3Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref)
{
    std::vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    
    //Loop on bin particles
    int npart_tot = particles.size();
    for (int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle
        (*this)(EMfields, particles, ipart, npart_tot, &(*Epart)[ipart], &(*Bpart)[ipart]);
        //Buffering of iol and delta
        (*iold)[ipart] = ip_;
        (*delta)[ipart] = xi;
    }
    
}

// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator1D3Order::operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> * selection)
{
    if( selection ) {
        
        int nsel_tot = selection->size();
        for (int isel=0 ; isel<nsel_tot; isel++ ) {
            (*this)(EMfields, particles, (*selection)[isel], offset, buffer+isel, buffer+isel+3*offset);
        }
        
    } else {
        
        int npart_tot = particles.size();
        for (int ipart=0 ; ipart<npart_tot; ipart++ ) {
            (*this)(EMfields, particles, ipart, offset, buffer+ipart, buffer+ipart+3*offset);
        }
        
    }
}


