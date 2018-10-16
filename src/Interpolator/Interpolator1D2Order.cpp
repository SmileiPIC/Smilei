#include "Interpolator1D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"

using namespace std;

Interpolator1D2Order::Interpolator1D2Order(Params &params, Patch* patch) : Interpolator1D(params, patch) {
    dx_inv_ = 1.0/params.cell_length[0];
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc)
{
    // Static cast of the electromagnetic fields
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);
    
    // Particle position (in units of the spatial-step)
    double xjn = particles.position(0, ipart)*dx_inv_;
    // Calculate coeffs
    coeffs(xjn);
    
    // Interpolate the fields from the Dual grid : Ex, By, Bz
    *(ELoc+0*nparts) = compute(coeffd_, Ex1D,   id_);
    *(BLoc+1*nparts) = compute(coeffd_, By1D_m, id_);
    *(BLoc+2*nparts) = compute(coeffd_, Bz1D_m, id_);
    
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    *(ELoc+1*nparts) = compute(coeffp_, Ey1D,   ip_);
    *(ELoc+2*nparts) = compute(coeffp_, Ez1D,   ip_);
    *(BLoc+0*nparts) = compute(coeffp_, Bx1D_m, ip_);
    
}//END Interpolator1D2Order

void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc)
{
    int ipart = *istart;

    double *ELoc = &(smpi->dynamics_Epart[ithread][ipart]);
    double *BLoc = &(smpi->dynamics_Bpart[ithread][ipart]);
    
    // Static cast of the electromagnetic fields
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);
    Field1D* Jx1D     = static_cast<Field1D*>(EMfields->Jx_);
    Field1D* Jy1D     = static_cast<Field1D*>(EMfields->Jy_);
    Field1D* Jz1D     = static_cast<Field1D*>(EMfields->Jz_);
    Field1D* Rho1D    = static_cast<Field1D*>(EMfields->rho_);
    
    // Particle position (in units of the spatial-step)
    double xjn = particles.position(0, ipart)*dx_inv_;
    // Calculate coeffs
    coeffs(xjn);
    
    int nparts( particles.size() );
    
    // Interpolate the fields from the Dual grid : Ex, By, Bz
    *(ELoc+0*nparts) = compute(coeffd_, Ex1D,   id_);
    *(BLoc+1*nparts) = compute(coeffd_, By1D_m, id_);
    *(BLoc+2*nparts) = compute(coeffd_, Bz1D_m, id_);
    
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    *(ELoc+1*nparts) = compute(coeffp_, Ey1D,   ip_);
    *(ELoc+2*nparts) = compute(coeffp_, Ez1D,   ip_);
    *(BLoc+0*nparts) = compute(coeffp_, Bx1D_m, ip_);
    
    // Interpolate the fields from the Primal grid : Jy, Jz, Rho
    (*JLoc).y = compute(coeffp_, Jy1D,  ip_);
    (*JLoc).z = compute(coeffp_, Jz1D,  ip_);
    (*RhoLoc) = compute(coeffp_, Rho1D, ip_);
    
    // Interpolate the fields from the Dual grid : Jx
    (*JLoc).x = compute(coeffd_, Jx1D,  id_);
    
}

// Interpolator on another field than the basic ones
void Interpolator1D2Order::operator() (Field* field, Particles &particles, int *istart, int *iend, double* FieldLoc)
{
    Field1D* F = static_cast<Field1D*>(field);
    double * coeff = field->isDual(0) ? coeffd_ : coeffp_;
    int * i = field->isDual(0) ? &id_ : &ip_;
    
    for (int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xjn = particles.position(0, ipart)*dx_inv_;
        coeffs(xjn);
        (*FieldLoc) = compute(coeff, F, *i);
    }
}

void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref)
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
        (*delta)[ipart] = xjmxi;
    }
    
}

// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> * selection)
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
