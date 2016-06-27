#include "Interpolator3D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3D2Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3D2Order::Interpolator3D2Order(Params &params, Patch *patch) : Interpolator3D(params, patch)
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];
    dz_inv_ = 1.0/params.cell_length[2];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
{
    // Static cast of the electromagnetic fields
    Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
    Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
    Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
    Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_m);
    Field3D* By3D = static_cast<Field3D*>(EMfields->By_m);
    Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_m);


    // Normalized particle position
    double xpn = particles.position(0, ipart)*dx_inv_;
    double ypn = particles.position(1, ipart)*dy_inv_;
    double zpn = particles.position(2, ipart)*dz_inv_;


    // Indexes of the central nodes
    ip_ = round(xpn);
    id_ = round(xpn+0.5);
    jp_ = round(ypn);
    jd_ = round(ypn+0.5);
    kp_ = round(zpn);
    kd_ = round(zpn+0.5);


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

    deltaz   = zpn - (double)kd_ + 0.5;
    delta2  = deltaz*deltaz;
    coeffzd_[0] = 0.5 * (delta2-deltaz+0.25);
    coeffzd_[1] = 0.75 - delta2;
    coeffzd_[2] = 0.5 * (delta2+deltaz+0.25);

    deltaz   = zpn - (double)kp_;
    delta2  = deltaz*deltaz;
    coeffzp_[0] = 0.5 * (delta2-deltaz+0.25);
    coeffzp_[1] = 0.75 - delta2;
    coeffzp_[2] = 0.5 * (delta2+deltaz+0.25);


    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    id_ = id_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;
    jd_ = jd_ - j_domain_begin;
    kp_ = kp_ - k_domain_begin;
    kd_ = kd_ - k_domain_begin;


    // -------------------------
    // Interpolation of Ex^(d,p,p)
    // -------------------------
    (*ELoc).x = compute( &coeffxd_[1], &coeffyp_[1], &coeffzp_[1], Ex3D, id_, jp_, kp_);

    // -------------------------
    // Interpolation of Ey^(p,d,p)
    // -------------------------
    (*ELoc).y = compute( &coeffxp_[1], &coeffyd_[1], &coeffzp_[1], Ey3D, ip_, jd_, kp_);

    // -------------------------
    // Interpolation of Ez^(p,p,d)
    // -------------------------
    (*ELoc).z = compute( &coeffxp_[1], &coeffyp_[1], &coeffzd_[1], Ez3D, ip_, jp_, kd_);

    // -------------------------
    // Interpolation of Bx^(p,d,d)
    // -------------------------
    (*BLoc).x = compute( &coeffxp_[1], &coeffyd_[1], &coeffzd_[1], Bx3D, ip_, jd_, kd_);

    // -------------------------
    // Interpolation of By^(d,p,d)
    // -------------------------
    (*BLoc).y = compute( &coeffxd_[1], &coeffyp_[1], &coeffzd_[1], By3D, id_, jp_, kd_);

    // -------------------------
    // Interpolation of Bz^(d,d,p)
    // -------------------------
    (*BLoc).z = compute( &coeffxd_[1], &coeffyd_[1], &coeffzp_[1], Bz3D, id_, jd_, kp_);

} // END Interpolator3D2Order

void Interpolator3D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
{
    // Interpolate E, B
    // Compute coefficient for ipart position
    (*this)(EMfields, particles, ipart, ELoc, BLoc);

    // Static cast of the electromagnetic fields
    Field3D* Jx3D = static_cast<Field3D*>(EMfields->Jx_);
    Field3D* Jy3D = static_cast<Field3D*>(EMfields->Jy_);
    Field3D* Jz3D = static_cast<Field3D*>(EMfields->Jz_);
    Field3D* Rho3D= static_cast<Field3D*>(EMfields->rho_);
    
    
    // -------------------------
    // Interpolation of Jx^(d,p,p)
    // -------------------------
    (*JLoc).x = compute( &coeffxd_[1], &coeffyp_[1], &coeffzp_[1], Jx3D, id_, jp_, kp_);
    
    // -------------------------
    // Interpolation of Jy^(p,d,p)
    // -------------------------
    (*JLoc).y = compute( &coeffxp_[1], &coeffyd_[1], &coeffzp_[1], Jx3D, ip_, jd_, kp_);
    
    // -------------------------
    // Interpolation of Jz^(p,p,d)
    // -------------------------
    (*JLoc).z = compute( &coeffxp_[1], &coeffyp_[1], &coeffzd_[1], Jz3D, ip_, jp_, kd_);
    
    // -------------------------
    // Interpolation of Rho^(p,p,p)
    // -------------------------
    (*RhoLoc) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Rho3D, ip_, jp_, kp_);

}

void Interpolator3D2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread)
{
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);

    //Loop on bin particles
    for (unsigned int ipart=istart ; ipart<iend; ipart++ ) {
        //Interpolation on current particle
        (*this)(EMfields, particles, ipart, &(*Epart)[ipart], &(*Bpart)[ipart]);
        //Buffering of iol and delta
        (*iold)[ipart*2+0] = ip_;
        (*iold)[ipart*2+1] = jp_;
        (*iold)[ipart*2+2] = kp_;
        (*delta)[ipart*2+0] = deltax;
        (*delta)[ipart*2+1] = deltay;
        (*delta)[ipart*2+2] = deltaz;
    }

}
