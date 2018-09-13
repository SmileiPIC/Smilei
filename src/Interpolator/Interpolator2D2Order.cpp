#include "Interpolator2D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "LaserEnvelope.h"

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
void Interpolator2D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc)
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
    
    
    // Indexes of the central nodes
    ip_ = round(xpn);
    id_ = round(xpn+0.5);
    jp_ = round(ypn);
    jd_ = round(ypn+0.5);
    
    
    // Declaration and calculation of the coefficient for interpolation
    double delta2;
    
    deltax   = xpn - (double)id_ + 0.5;
    delta2  = deltax*deltax;
    coeffxd_[0] = 0.5 * (delta2-deltax+0.25);
    coeffxd_[1] = 0.75 - delta2;
    coeffxd_[2] = 0.5 * (delta2+deltax+0.25);
    
    deltax   = xpn - (double)ip_;
    delta2  = deltax*deltax;
    coeffxp_[0] = 0.5 * (delta2-deltax+0.25);
    coeffxp_[1] = 0.75 - delta2;
    coeffxp_[2] = 0.5 * (delta2+deltax+0.25);
    
    deltay   = ypn - (double)jd_ + 0.5;
    delta2  = deltay*deltay;
    coeffyd_[0] = 0.5 * (delta2-deltay+0.25);
    coeffyd_[1] = 0.75 - delta2;
    coeffyd_[2] = 0.5 * (delta2+deltay+0.25);
    
    deltay   = ypn - (double)jp_;
    delta2  = deltay*deltay;
    coeffyp_[0] = 0.5 * (delta2-deltay+0.25);
    coeffyp_[1] = 0.75 - delta2;
    coeffyp_[2] = 0.5 * (delta2+deltay+0.25);
    
    
    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    id_ = id_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;
    jd_ = jd_ - j_domain_begin;
    
    // -------------------------
    // Interpolation of Ex^(d,p)
    // -------------------------
    *(ELoc+0*nparts) =  compute( &coeffxd_[1], &coeffyp_[1], Ex2D, id_, jp_);
    
    // -------------------------
    // Interpolation of Ey^(p,d)
    // -------------------------
    *(ELoc+1*nparts) = compute( &coeffxp_[1], &coeffyd_[1], Ey2D, ip_, jd_);
    
    // -------------------------
    // Interpolation of Ez^(p,p)
    // -------------------------
    *(ELoc+2*nparts) = compute( &coeffxp_[1], &coeffyp_[1], Ez2D, ip_, jp_);
    
    // -------------------------
    // Interpolation of Bx^(p,d)
    // -------------------------
    *(BLoc+0*nparts) = compute( &coeffxp_[1], &coeffyd_[1], Bx2D, ip_, jd_);
    
    // -------------------------
    // Interpolation of By^(d,p)
    // -------------------------
    *(BLoc+1*nparts) = compute( &coeffxd_[1], &coeffyp_[1], By2D, id_, jp_);
    
    // -------------------------
    // Interpolation of Bz^(d,d)
    // -------------------------
     *(BLoc+2*nparts) = compute( &coeffxd_[1], &coeffyd_[1], Bz2D, id_, jd_);

} // END Interpolator2D2Order

void Interpolator2D2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc)
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
    
    
    // Normalized particle position
    double xpn = particles.position(0, ipart)*dx_inv_;
    double ypn = particles.position(1, ipart)*dy_inv_;
    
    
    // Indexes of the central nodes
    ip_ = round(xpn);
    id_ = round(xpn+0.5);
    jp_ = round(ypn);
    jd_ = round(ypn+0.5);
    
    
    // Declaration and calculation of the coefficient for interpolation
    double delta2;
    
    deltax   = xpn - (double)id_ + 0.5;
    delta2  = deltax*deltax;
    coeffxd_[0] = 0.5 * (delta2-deltax+0.25);
    coeffxd_[1] = 0.75 - delta2;
    coeffxd_[2] = 0.5 * (delta2+deltax+0.25);
    
    deltax   = xpn - (double)ip_;
    delta2  = deltax*deltax;
    coeffxp_[0] = 0.5 * (delta2-deltax+0.25);
    coeffxp_[1] = 0.75 - delta2;
    coeffxp_[2] = 0.5 * (delta2+deltax+0.25);
    
    deltay   = ypn - (double)jd_ + 0.5;
    delta2  = deltay*deltay;
    coeffyd_[0] = 0.5 * (delta2-deltay+0.25);
    coeffyd_[1] = 0.75 - delta2;
    coeffyd_[2] = 0.5 * (delta2+deltay+0.25);
    
    deltay   = ypn - (double)jp_;
    delta2  = deltay*deltay;
    coeffyp_[0] = 0.5 * (delta2-deltay+0.25);
    coeffyp_[1] = 0.75 - delta2;
    coeffyp_[2] = 0.5 * (delta2+deltay+0.25);
    
    
    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    id_ = id_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;
    jd_ = jd_ - j_domain_begin;
    
    int nparts( particles.size() );    

    // -------------------------
    // Interpolation of Ex^(d,p)
    // -------------------------
    *(ELoc+0*nparts) =  compute( &coeffxd_[1], &coeffyp_[1], Ex2D, id_, jp_);
    
    // -------------------------
    // Interpolation of Ey^(p,d)
    // -------------------------
    *(ELoc+1*nparts) = compute( &coeffxp_[1], &coeffyd_[1], Ey2D, ip_, jd_);
    
    // -------------------------
    // Interpolation of Ez^(p,p)
    // -------------------------
    *(ELoc+2*nparts) = compute( &coeffxp_[1], &coeffyp_[1], Ez2D, ip_, jp_);
    
    // -------------------------
    // Interpolation of Bx^(p,d)
    // -------------------------
    *(BLoc+0*nparts) = compute( &coeffxp_[1], &coeffyd_[1], Bx2D, ip_, jd_);
    
    // -------------------------
    // Interpolation of By^(d,p)
    // -------------------------
    *(BLoc+1*nparts) = compute( &coeffxd_[1], &coeffyp_[1], By2D, id_, jp_);
    
    // -------------------------
    // Interpolation of Bz^(d,d)
    // -------------------------
     *(BLoc+2*nparts) = compute( &coeffxd_[1], &coeffyd_[1], Bz2D, id_, jd_);
    
    // Static cast of the electromagnetic fields
    Field2D* Jx2D = static_cast<Field2D*>(EMfields->Jx_);
    Field2D* Jy2D = static_cast<Field2D*>(EMfields->Jy_);
    Field2D* Jz2D = static_cast<Field2D*>(EMfields->Jz_);
    Field2D* Rho2D= static_cast<Field2D*>(EMfields->rho_);
    
    
    // -------------------------
    // Interpolation of Jx^(d,p)
    // -------------------------
    (*JLoc).x = compute( &coeffxd_[1], &coeffyp_[1], Jx2D, id_, jp_);
    
    // -------------------------
    // Interpolation of Jy^(p,d)
    // -------------------------
    (*JLoc).y = compute( &coeffxp_[1], &coeffyd_[1], Jy2D, ip_, jd_);
    
    // -------------------------
    // Interpolation of Ez^(p,p)
    // -------------------------
    (*JLoc).z = compute( &coeffxp_[1], &coeffyp_[1], Jz2D, ip_, jp_);
    
    // -------------------------
    // Interpolation of Rho^(p,p)
    // -------------------------
    (*RhoLoc) = compute( &coeffxp_[1], &coeffyp_[1], Rho2D, ip_, jp_);

}

void Interpolator2D2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref)
{
    std::vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    
    //Loop on bin particles
    int nparts( particles.size() );
    for (int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle
        (*this)(EMfields, particles, ipart, nparts, &(*Epart)[ipart], &(*Bpart)[ipart]);
        //Buffering of iol and delta
        (*iold)[ipart+0*nparts]  = ip_;
        (*iold)[ipart+1*nparts]  = jp_;
        (*delta)[ipart+0*nparts] = deltax;
        (*delta)[ipart+1*nparts] = deltay;
    }

}

// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator2D2Order::operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> * selection)
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

void Interpolator2D2Order::interpolate_em_fields_and_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the envelope fields
    Field2D* Phi2D = static_cast<Field2D*>(EMfields->envelope->Phi_);
    Field2D* GradPhix2D = static_cast<Field2D*>(EMfields->envelope->GradPhix_);
    Field2D* GradPhiy2D = static_cast<Field2D*>(EMfields->envelope->GradPhiy_);
    Field2D* GradPhiz2D = static_cast<Field2D*>(EMfields->envelope->GradPhiz_);

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
        (*PHIpart)[ipart] = compute( &coeffxp_[1], &coeffyp_[1], Phi2D, ip_, jp_);

        // -------------------------
        // Interpolation of GradPhix^(p,p)
        // -------------------------
        (*GradPHIpart)[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhix2D, ip_, jp_);

        // -------------------------
        // Interpolation of GradPhiy^(p,p)
        // -------------------------
        (*GradPHIpart)[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiy2D, ip_, jp_);
    
        // -------------------------
        // Interpolation of GradPhiz^(p,p,p)
        // -------------------------
        (*GradPHIpart)[ipart+2*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiz2D, ip_, jp_);


        //Buffering of iold and delta
        (*iold)[ipart+0*nparts]  = ip_;
        (*iold)[ipart+1*nparts]  = jp_;
        (*delta)[ipart+0*nparts] = deltax;
        (*delta)[ipart+1*nparts] = deltay;

    }


} // END Interpolator2D2Order


void Interpolator2D2Order::interpolate_envelope_and_old_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the electromagnetic fields
    Field2D* Phi2D = static_cast<Field2D*>(EMfields->envelope->Phi_);
    Field2D* Phiold2D = static_cast<Field2D*>(EMfields->envelope->Phiold_);
    Field2D* GradPhix2D = static_cast<Field2D*>(EMfields->envelope->GradPhix_);
    Field2D* GradPhiy2D = static_cast<Field2D*>(EMfields->envelope->GradPhiy_);
    Field2D* GradPhiz2D = static_cast<Field2D*>(EMfields->envelope->GradPhiz_);
    Field2D* GradPhixold2D = static_cast<Field2D*>(EMfields->envelope->GradPhixold_);
    Field2D* GradPhiyold2D = static_cast<Field2D*>(EMfields->envelope->GradPhiyold_);
    Field2D* GradPhizold2D = static_cast<Field2D*>(EMfields->envelope->GradPhizold_);
    
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
        double ypn = particles.position(1, ipart)*dy_inv_;
    
    
        // Indexes of the central nodes
        ip_ = round(xpn);
        id_ = round(xpn+0.5);
        jp_ = round(ypn);
        jd_ = round(ypn+0.5);
        
    
        // Declaration and calculation of the coefficient for interpolation
        double delta2;
    
        deltax   = xpn - (double)id_ + 0.5;
        delta2  = deltax*deltax;
        coeffxd_[0] = 0.5 * (delta2-deltax+0.25);
        coeffxd_[1] = 0.75 - delta2;
        coeffxd_[2] = 0.5 * (delta2+deltax+0.25);
    
        deltax   = xpn - (double)ip_;
        delta2  = deltax*deltax;
        coeffxp_[0] = 0.5 * (delta2-deltax+0.25);
        coeffxp_[1] = 0.75 - delta2;
        coeffxp_[2] = 0.5 * (delta2+deltax+0.25);
    
        deltay   = ypn - (double)jd_ + 0.5;
        delta2  = deltay*deltay;
        coeffyd_[0] = 0.5 * (delta2-deltay+0.25);
        coeffyd_[1] = 0.75 - delta2;
        coeffyd_[2] = 0.5 * (delta2+deltay+0.25);
    
        deltay   = ypn - (double)jp_;
        delta2  = deltay*deltay;
        coeffyp_[0] = 0.5 * (delta2-deltay+0.25);
        coeffyp_[1] = 0.75 - delta2;
        coeffyp_[2] = 0.5 * (delta2+deltay+0.25);
    
    
    
        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        id_ = id_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
        jd_ = jd_ - j_domain_begin;
        
    
        // -------------------------
        // Interpolation of Phi^(p,p,p)
        // -------------------------
        (*PHIpart)[ipart] = compute( &coeffxp_[1], &coeffyp_[1], Phi2D, ip_, jp_);
    
        // -------------------------
        // Interpolation of Phiold^(p,p,p)
        // -------------------------
        (*PHIoldpart)[ipart] = compute( &coeffxp_[1], &coeffyp_[1], Phiold2D, ip_, jp_);
    
        // -------------------------
        // Interpolation of GradPhix^(p,p,p)
        // -------------------------
        (*GradPHIpart)[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhix2D, ip_, jp_);
    
        // -------------------------
        // Interpolation of GradPhixold^(p,p,p)
        // -------------------------
        (*GradPHIoldpart)[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhixold2D, ip_, jp_);
    
        // -------------------------
        // Interpolation of GradPhiy^(p,p,p)
        // -------------------------
        (*GradPHIpart)[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiy2D, ip_, jp_);
    
        // -------------------------
        // Interpolation of GradPhiyold^(p,p,p)
        // -------------------------
        (*GradPHIoldpart)[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiyold2D, ip_, jp_);
    
        // -------------------------
        // Interpolation of GradPhiz^(p,p,p)
        // -------------------------
        (*GradPHIpart)[ipart+2*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhiz2D, ip_, jp_);
    
        // -------------------------
        // Interpolation of GradPhizold^(p,p,p)
        // -------------------------
        (*GradPHIoldpart)[ipart+2*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhizold2D, ip_, jp_);
    
        //Buffering of iold and delta
        (*iold)[ipart+0*nparts]  = ip_;
        (*iold)[ipart+1*nparts]  = jp_;
        (*delta)[ipart+0*nparts] = deltax;
        (*delta)[ipart+1*nparts] = deltay;
        
    
    }

} // END Interpolator2D2Order


void Interpolator2D2Order::interpolate_envelope_and_susceptibility(ElectroMagn* EMfields, Particles &particles, int ipart, double* Env_A_abs_Loc, double* Env_Chi_Loc, double* Env_E_abs_Loc)
{
    // // Static cast of the electromagnetic fields
    // Field3D* Env_A_abs_3D = static_cast<Field3D*>(EMfields->Env_A_abs_);
    // Field3D* Env_Chi_3D = static_cast<Field3D*>(EMfields->Env_Chi_);
    // Field3D* Env_E_abs_3D = static_cast<Field3D*>(EMfields->Env_E_abs_);
    // 
    // // Normalized particle position
    // double xpn = particles.position(0, ipart)*dx_inv_;
    // double ypn = particles.position(1, ipart)*dy_inv_;
    // double zpn = particles.position(2, ipart)*dz_inv_;
    // 
    // 
    // // Indexes of the central nodes
    // ip_ = round(xpn);
    // jp_ = round(ypn);
    // kp_ = round(zpn);
    // 
    // 
    // // Declaration and calculation of the coefficient for interpolation
    // double delta2;
    // 
    // 
    // deltax   = xpn - (double)ip_;
    // delta2  = deltax*deltax;
    // coeffxp_[0] = 0.5 * (delta2-deltax+0.25);
    // coeffxp_[1] = 0.75 - delta2;
    // coeffxp_[2] = 0.5 * (delta2+deltax+0.25);
    // 
    // deltay   = ypn - (double)jp_;
    // delta2  = deltay*deltay;
    // coeffyp_[0] = 0.5 * (delta2-deltay+0.25);
    // coeffyp_[1] = 0.75 - delta2;
    // coeffyp_[2] = 0.5 * (delta2+deltay+0.25);
    // 
    // deltaz   = zpn - (double)kp_;
    // delta2  = deltaz*deltaz;
    // coeffzp_[0] = 0.5 * (delta2-deltaz+0.25);
    // coeffzp_[1] = 0.75 - delta2;
    // coeffzp_[2] = 0.5 * (delta2+deltaz+0.25);
    // 
    // 
    // //!\todo CHECK if this is correct for both primal & dual grids !!!
    // // First index for summation
    // ip_ = ip_ - i_domain_begin;
    // jp_ = jp_ - j_domain_begin;
    // kp_ = kp_ - k_domain_begin;
    // 
    // // -------------------------
    // // Interpolation of Env_A_abs_^(p,p,p)
    // // -------------------------
    // *(Env_A_abs_Loc) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_A_abs_3D, ip_, jp_, kp_);
    // 
    // // -------------------------
    // // Interpolation of Env_Chi_^(p,p,p)
    // // -------------------------
    // *(Env_Chi_Loc) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_Chi_3D, ip_, jp_, kp_);
    // 
    // // -------------------------
    // // Interpolation of Env_E_abs_^(p,p,p)
    // // -------------------------
    // *(Env_E_abs_Loc) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_E_abs_3D, ip_, jp_, kp_);

} // END Interpolator2D2Order

