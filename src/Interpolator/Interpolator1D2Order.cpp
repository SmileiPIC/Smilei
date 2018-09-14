#include "Interpolator1D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "LaserEnvelope.h"


using namespace std;

Interpolator1D2Order::Interpolator1D2Order(Params &params, Patch* patch) : Interpolator1D(params, patch) {
    dx_inv_ = 1.0/params.cell_length[0];
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc)
{
    
    // Variable declaration
    double xjn, xjmxi2;
    
    // Static cast of the electromagnetic fields
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);
    
    
    // Particle position (in units of the spatial-step)
    xjn    = particles.position(0, ipart)*dx_inv_;
    
    // --------------------------------------------------------
    // Interpolate the fields from the Dual grid : Ex, By, Bz
    // --------------------------------------------------------
    id_      = round(xjn+0.5);        // index of the central point
    xjmxi  = xjn - (double)id_ +0.5;  // normalized distance to the central node
    xjmxi2 = xjmxi*xjmxi;            // square of the normalized distance to the central node
    
    // 2nd order interpolation on 3 nodes
    coeffd_[0] = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffd_[1] = (0.75-xjmxi2);
    coeffd_[2] = 0.5 * (xjmxi2+xjmxi+0.25);
    
    id_ -= index_domain_begin;
    
    *(ELoc+0*nparts) = compute(coeffd_, Ex1D,   id_);  
    *(BLoc+1*nparts) = compute(coeffd_, By1D_m, id_);  
    *(BLoc+2*nparts) = compute(coeffd_, Bz1D_m, id_);  
    
    // --------------------------------------------------------
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    // --------------------------------------------------------
    ip_      = round(xjn);      // index of the central point
    xjmxi  = xjn -(double)ip_;  // normalized distance to the central node
    xjmxi2 = pow(xjmxi,2);      // square of the normalized distance to the central node
    
    // 2nd order interpolation on 3 nodes
    coeffp_[0] = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffp_[1] = (0.75-xjmxi2);
    coeffp_[2] = 0.5 * (xjmxi2+xjmxi+0.25);
    
    ip_ -= index_domain_begin;
    
    *(ELoc+1*nparts) = compute(coeffp_, Ey1D,   ip_);  
    *(ELoc+2*nparts) = compute(coeffp_, Ez1D,   ip_);  
    *(BLoc+0*nparts) = compute(coeffp_, Bx1D_m, ip_);  
    
}//END Interpolator1D2Order

void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc)
{
    int ipart = *istart;

    double *ELoc = &(smpi->dynamics_Epart[ithread][ipart]);
    double *BLoc = &(smpi->dynamics_Bpart[ithread][ipart]);

    // Interpolate E, B
    // Compute coefficient for ipart position
    // Variable declaration
    double xjn, xjmxi2;
    
    // Static cast of the electromagnetic fields
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);
    
    
    // Particle position (in units of the spatial-step)
    xjn    = particles.position(0, ipart)*dx_inv_;
    
    // --------------------------------------------------------
    // Interpolate the fields from the Dual grid : Ex, By, Bz
    // --------------------------------------------------------
    id_      = round(xjn+0.5);        // index of the central point
    xjmxi  = xjn - (double)id_ +0.5;  // normalized distance to the central node
    xjmxi2 = xjmxi*xjmxi;            // square of the normalized distance to the central node
    
    // 2nd order interpolation on 3 nodes
    coeffd_[0] = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffd_[1] = (0.75-xjmxi2);
    coeffd_[2] = 0.5 * (xjmxi2+xjmxi+0.25);
    
    id_ -= index_domain_begin;
    
    int nparts( particles.size() );
    
    *(ELoc+0*nparts) = compute(coeffd_, Ex1D,   id_);  
    *(BLoc+1*nparts) = compute(coeffd_, By1D_m, id_);  
    *(BLoc+2*nparts) = compute(coeffd_, Bz1D_m, id_);  
    
    // --------------------------------------------------------
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    // --------------------------------------------------------
    ip_      = round(xjn);      // index of the central point
    xjmxi  = xjn -(double)ip_;  // normalized distance to the central node
    xjmxi2 = pow(xjmxi,2);      // square of the normalized distance to the central node
    
    // 2nd order interpolation on 3 nodes
    coeffp_[0] = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffp_[1] = (0.75-xjmxi2);
    coeffp_[2] = 0.5 * (xjmxi2+xjmxi+0.25);
    
    ip_ -= index_domain_begin;
    
    *(ELoc+1*nparts) = compute(coeffp_, Ey1D,   ip_);  
    *(ELoc+2*nparts) = compute(coeffp_, Ez1D,   ip_);  
    *(BLoc+0*nparts) = compute(coeffp_, Bx1D_m, ip_);  
    
    // Static cast of the electromagnetic fields
    Field1D* Jx1D     = static_cast<Field1D*>(EMfields->Jx_);
    Field1D* Jy1D     = static_cast<Field1D*>(EMfields->Jy_);
    Field1D* Jz1D     = static_cast<Field1D*>(EMfields->Jz_);
    Field1D* Rho1D    = static_cast<Field1D*>(EMfields->rho_);
    
    // --------------------------------------------------------
    // Interpolate the fields from the Primal grid : Jy, Jz, Rho
    // --------------------------------------------------------
    (*JLoc).y = compute(coeffp_, Jy1D,  ip_);  
    (*JLoc).z = compute(coeffp_, Jz1D,  ip_);  
    (*RhoLoc) = compute(coeffp_, Rho1D, ip_);    
    
    // --------------------------------------------------------
    // Interpolate the fields from the Dual grid : Jx
    // --------------------------------------------------------
    (*JLoc).x = compute(coeffd_, Jx1D,  id_);  
    
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

void Interpolator1D2Order::interpolate_em_fields_and_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the envelope fields
    Field1D* Phi1D = static_cast<Field1D*>(EMfields->envelope->Phi_);
    Field1D* GradPhix1D = static_cast<Field1D*>(EMfields->envelope->GradPhix_);
    Field1D* GradPhiy1D = static_cast<Field1D*>(EMfields->envelope->GradPhiy_);
    Field1D* GradPhiz1D = static_cast<Field1D*>(EMfields->envelope->GradPhiz_);

    std::vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<double> *PHIpart        = &(smpi->dynamics_PHIpart[ithread]);
    std::vector<double> *GradPHIpart    = &(smpi->dynamics_GradPHIpart[ithread]);

    std::vector<int>    *iold  = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);

    //Loop on bin particles
    int nparts( particles.size() );
    for (int ipart=*istart ; ipart<*iend; ipart++ ) {

        (*this)(EMfields, particles, ipart, nparts, &(*Epart)[ipart], &(*Bpart)[ipart]);


        // -------------------------
        // Interpolation of Phi^(p,p)
        // -------------------------
        (*PHIpart)[ipart] = compute( &coeffp_[1], Phi1D, ip_);

        // -------------------------
        // Interpolation of GradPhix^(p,p)
        // -------------------------
        (*GradPHIpart)[ipart+0*nparts] = compute( &coeffp_[1], GradPhix1D, ip_);

        // -------------------------
        // Interpolation of GradPhiy^(p,p)
        // -------------------------
        (*GradPHIpart)[ipart+1*nparts] = compute( &coeffp_[1], GradPhiy1D, ip_);
    
        // -------------------------
        // Interpolation of GradPhiz^(p,p)
        // -------------------------
        (*GradPHIpart)[ipart+2*nparts] = compute( &coeffp_[1], GradPhiz1D, ip_);


        //Buffering of iold and delta
        (*iold)[ipart+0*nparts]  = ip_;
        (*delta)[ipart+0*nparts] = deltax;

    }


} // END Interpolator1D2Order


void Interpolator1D2Order::interpolate_envelope_and_old_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the electromagnetic fields
    Field1D* Phi1D = static_cast<Field1D*>(EMfields->envelope->Phi_);
    Field1D* Phiold1D = static_cast<Field1D*>(EMfields->envelope->Phiold_);
    Field1D* GradPhix1D = static_cast<Field1D*>(EMfields->envelope->GradPhix_);
    Field1D* GradPhiy1D = static_cast<Field1D*>(EMfields->envelope->GradPhiy_);
    Field1D* GradPhiz1D = static_cast<Field1D*>(EMfields->envelope->GradPhiz_);
    Field1D* GradPhixold1D = static_cast<Field1D*>(EMfields->envelope->GradPhixold_);
    Field1D* GradPhiyold1D = static_cast<Field1D*>(EMfields->envelope->GradPhiyold_);
    Field1D* GradPhizold1D = static_cast<Field1D*>(EMfields->envelope->GradPhizold_);
    
    std::vector<double> *PHIpart        = &(smpi->dynamics_PHIpart[ithread]);
    std::vector<double> *GradPHIpart    = &(smpi->dynamics_GradPHIpart[ithread]);
    std::vector<double> *PHIoldpart     = &(smpi->dynamics_PHIoldpart[ithread]);
    std::vector<double> *GradPHIoldpart = &(smpi->dynamics_GradPHIoldpart[ithread]);
    
    std::vector<int>    *iold  = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    
    //Loop on bin particles
    int nparts( particles.size() );
    for (int ipart=*istart ; ipart<*iend; ipart++ ) {
    
        // Normalized particle position
        double xpn = particles.position(0, ipart)*dx_inv_;
    
        // Indexes of the central nodes
        ip_ = round(xpn);
        
        // Declaration and calculation of the coefficient for interpolation
        double deltax,delta2;
    
    
        deltax   = xpn - (double)ip_;
        delta2  = deltax*deltax;
        coeffp_[0] = 0.5 * (delta2-deltax+0.25);
        coeffp_[1] = 0.75 - delta2;
        coeffp_[2] = 0.5 * (delta2+deltax+0.25);
    
    
    
        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - index_domain_begin;
        
    
        // -------------------------
        // Interpolation of Phi^(p)
        // -------------------------
        (*PHIpart)[ipart] = compute( &coeffp_[1], Phi1D, ip_);
    
        // -------------------------
        // Interpolation of Phiold^(p)
        // -------------------------
        (*PHIoldpart)[ipart] = compute( &coeffp_[1], Phiold1D, ip_);
    
        // -------------------------
        // Interpolation of GradPhix^(p)
        // -------------------------
        (*GradPHIpart)[ipart+0*nparts] = compute( &coeffp_[1], GradPhix1D, ip_);
    
        // -------------------------
        // Interpolation of GradPhixold^(p)
        // -------------------------
        (*GradPHIoldpart)[ipart+0*nparts] = compute( &coeffp_[1], GradPhixold1D, ip_);
    
        // -------------------------
        // Interpolation of GradPhiy^(p)
        // -------------------------
        (*GradPHIpart)[ipart+1*nparts] = compute( &coeffp_[1], GradPhiy1D, ip_);
    
        // -------------------------
        // Interpolation of GradPhiyold^(p)
        // -------------------------
        (*GradPHIoldpart)[ipart+1*nparts] = compute( &coeffp_[1], GradPhiyold1D, ip_);
    
        // -------------------------
        // Interpolation of GradPhiz^(p)
        // -------------------------
        (*GradPHIpart)[ipart+2*nparts] = compute( &coeffp_[1], GradPhiz1D, ip_);
    
        // -------------------------
        // Interpolation of GradPhizold^(p)
        // -------------------------
        (*GradPHIoldpart)[ipart+2*nparts] = compute( &coeffp_[1], GradPhizold1D, ip_);
    
        //Buffering of iold and delta
        (*iold)[ipart+0*nparts]  = ip_;
        (*delta)[ipart+0*nparts] = deltax;
        
    
    }

} // END Interpolator1D2Order


void Interpolator1D2Order::interpolate_envelope_and_susceptibility(ElectroMagn* EMfields, Particles &particles, int ipart, double* Env_A_abs_Loc, double* Env_Chi_Loc, double* Env_E_abs_Loc)
{
    // Static cast of the electromagnetic fields
    Field1D* Env_A_abs_1D = static_cast<Field1D*>(EMfields->Env_A_abs_);
    Field1D* Env_Chi_1D = static_cast<Field1D*>(EMfields->Env_Chi_);
    Field1D* Env_E_abs_1D = static_cast<Field1D*>(EMfields->Env_E_abs_);
    
    // Normalized particle position
    double xpn = particles.position(0, ipart)*dx_inv_;
    
    // Indexes of the central nodes
    ip_ = round(xpn);

    // Declaration and calculation of the coefficient for interpolation
    double deltax,delta2;
    
    
    deltax   = xpn - (double)ip_;
    delta2  = deltax*deltax;
    coeffp_[0] = 0.5 * (delta2-deltax+0.25);
    coeffp_[1] = 0.75 - delta2;
    coeffp_[2] = 0.5 * (delta2+deltax+0.25);
    
  
    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - index_domain_begin;
    
    // -------------------------
    // Interpolation of Env_A_abs_^(p)
    // -------------------------
    *(Env_A_abs_Loc) = compute( &coeffp_[1], Env_A_abs_1D, ip_);
    
    // -------------------------
    // Interpolation of Env_Chi_^(p)
    // -------------------------
    *(Env_Chi_Loc) = compute( &coeffp_[1], Env_Chi_1D, ip_);
    
    // -------------------------
    // Interpolation of Env_E_abs_^(p)
    // -------------------------
    *(Env_E_abs_Loc) = compute( &coeffp_[1], Env_E_abs_1D, ip_);

} // END Interpolator1D2Order

