#include "Interpolator2D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D2Order::Interpolator2D2Order(Params &params, Patch *patch) : Interpolator2D(params, patch)
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D2Order::fields(ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc)
{
    // Static cast of the electromagnetic fields
    Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_m);
    Field2D* By2D = static_cast<Field2D*>(EMfields->By_m);
    Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_m);
    
    
    // Normalized particle position
    double xpn = particles.position(0, ipart)*dx_inv_;
    double ypn = particles.position(1, ipart)*dy_inv_;
    // Calculate coeffs
    coeffs(xpn, ypn);
    
    // Interpolation of Ex^(d,p)
    *(ELoc+0*nparts) = compute( &coeffxd_[1], &coeffyp_[1], Ex2D, id_, jp_);
    // Interpolation of Ey^(p,d)
    *(ELoc+1*nparts) = compute( &coeffxp_[1], &coeffyd_[1], Ey2D, ip_, jd_);
    // Interpolation of Ez^(p,p)
    *(ELoc+2*nparts) = compute( &coeffxp_[1], &coeffyp_[1], Ez2D, ip_, jp_);
    // Interpolation of Bx^(p,d)
    *(BLoc+0*nparts) = compute( &coeffxp_[1], &coeffyd_[1], Bx2D, ip_, jd_);
    // Interpolation of By^(d,p)
    *(BLoc+1*nparts) = compute( &coeffxd_[1], &coeffyp_[1], By2D, id_, jp_);
    // Interpolation of Bz^(d,d)
    *(BLoc+2*nparts) = compute( &coeffxd_[1], &coeffyd_[1], Bz2D, id_, jd_);
} // END Interpolator2D2Order

void Interpolator2D2Order::fieldsAndCurrents(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc)
{
    int ipart = *istart;

    double *ELoc = &(smpi->dynamics_Epart[ithread][ipart]);
    double *BLoc = &(smpi->dynamics_Bpart[ithread][ipart]);

    // Interpolate E, B
    // Compute coefficient for ipart position
    // Static cast of the electromagnetic fields
    Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_m);
    Field2D* By2D = static_cast<Field2D*>(EMfields->By_m);
    Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_m);
    Field2D* Jx2D = static_cast<Field2D*>(EMfields->Jx_);
    Field2D* Jy2D = static_cast<Field2D*>(EMfields->Jy_);
    Field2D* Jz2D = static_cast<Field2D*>(EMfields->Jz_);
    Field2D* Rho2D= static_cast<Field2D*>(EMfields->rho_);
    
    // Normalized particle position
    double xpn = particles.position(0, ipart)*dx_inv_;
    double ypn = particles.position(1, ipart)*dy_inv_;
    // Calculate coeffs
    coeffs(xpn, ypn);
    
    int nparts( particles.size() );
    
    // Interpolation of Ex^(d,p)
    *(ELoc+0*nparts) = compute( &coeffxd_[1], &coeffyp_[1], Ex2D, id_, jp_);
    // Interpolation of Ey^(p,d)
    *(ELoc+1*nparts) = compute( &coeffxp_[1], &coeffyd_[1], Ey2D, ip_, jd_);
    // Interpolation of Ez^(p,p)
    *(ELoc+2*nparts) = compute( &coeffxp_[1], &coeffyp_[1], Ez2D, ip_, jp_);
    // Interpolation of Bx^(p,d)
    *(BLoc+0*nparts) = compute( &coeffxp_[1], &coeffyd_[1], Bx2D, ip_, jd_);
    // Interpolation of By^(d,p)
    *(BLoc+1*nparts) = compute( &coeffxd_[1], &coeffyp_[1], By2D, id_, jp_);
    // Interpolation of Bz^(d,d)
    *(BLoc+2*nparts) = compute( &coeffxd_[1], &coeffyd_[1], Bz2D, id_, jd_);
    // Interpolation of Jx^(d,p)
    (*JLoc).x = compute( &coeffxd_[1], &coeffyp_[1], Jx2D, id_, jp_);
    // Interpolation of Jy^(p,d)
    (*JLoc).y = compute( &coeffxp_[1], &coeffyd_[1], Jy2D, ip_, jd_);
    // Interpolation of Jz^(p,p)
    (*JLoc).z = compute( &coeffxp_[1], &coeffyp_[1], Jz2D, ip_, jp_);
    // Interpolation of Rho^(p,p)
    (*RhoLoc) = compute( &coeffxp_[1], &coeffyp_[1], Rho2D, ip_, jp_);
}

// Interpolator on another field than the basic ones
void Interpolator2D2Order::oneField(Field* field, Particles &particles, int *istart, int *iend, double* FieldLoc)
{
    Field2D* F = static_cast<Field2D*>(field);
    double * coeffx = field->isDual(0) ? &coeffxd_[1] : &coeffxp_[1];
    double * coeffy = field->isDual(1) ? &coeffyd_[1] : &coeffyp_[1];
    int * i = field->isDual(0) ? &id_ : &ip_;
    int * j = field->isDual(1) ? &jd_ : &jp_;
    
    for (int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position(0, ipart)*dx_inv_;
        double ypn = particles.position(1, ipart)*dy_inv_;
        coeffs(xpn, ypn);
        FieldLoc[ipart] = compute(coeffx, coeffy, F, *i, *j);
    }
}

void Interpolator2D2Order::fields_batch(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref)
{
    std::vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    
    //Loop on bin particles
    int nparts( particles.size() );
    for (int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle
        fields(EMfields, particles, ipart, nparts, &(*Epart)[ipart], &(*Bpart)[ipart]);
        //Buffering of iol and delta
        (*iold)[ipart+0*nparts]  = ip_;
        (*iold)[ipart+1*nparts]  = jp_;
        (*delta)[ipart+0*nparts] = deltax;
        (*delta)[ipart+1*nparts] = deltay;
    }

}

// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator2D2Order::fields_selection(ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> * selection)
{
    if( selection ) {
        
        int nsel_tot = selection->size();
        for (int isel=0 ; isel<nsel_tot; isel++ ) {
            fields(EMfields, particles, (*selection)[isel], offset, buffer+isel, buffer+isel+3*offset);
        }
        
    } else {
        
        int npart_tot = particles.size();
        for (int ipart=0 ; ipart<npart_tot; ipart++ ) {
            fields(EMfields, particles, ipart, offset, buffer+ipart, buffer+ipart+3*offset);
        }
        
    }
}
