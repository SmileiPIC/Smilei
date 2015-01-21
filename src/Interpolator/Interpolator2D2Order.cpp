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
Interpolator2D2Order::Interpolator2D2Order(PicParams &params, SmileiMPI *smpi) : Interpolator2D(params, smpi)
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
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
    double delta, delta2;

    delta   = xpn - (double)ip_;
    delta2  = delta*delta;
    coeffxp_[0] = 0.5 * (delta2-delta+0.25);
    coeffxp_[1] = 0.75 - delta2;
    coeffxp_[2] = 0.5 * (delta2+delta+0.25);

    delta   = xpn - (double)id_ + 0.5;
    delta2  = delta*delta;
    coeffxd_[0] = 0.5 * (delta2-delta+0.25);
    coeffxd_[1] = 0.75 - delta2;
    coeffxd_[2] = 0.5 * (delta2+delta+0.25);

    delta   = ypn - (double)jp_;
    delta2  = delta*delta;
    coeffyp_[0] = 0.5 * (delta2-delta+0.25);
    coeffyp_[1] = 0.75 - delta2;
    coeffyp_[2] = 0.5 * (delta2+delta+0.25);

    delta   = ypn - (double)jd_ + 0.5;
    delta2  = delta*delta;
    coeffyd_[0] = 0.5 * (delta2-delta+0.25);
    coeffyd_[1] = 0.75 - delta2;
    coeffyd_[2] = 0.5 * (delta2+delta+0.25);

    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - 1 - i_domain_begin;
    id_ = id_ - 1 - i_domain_begin;
    jp_ = jp_ - 1 - j_domain_begin;
    jd_ = jd_ - 1 - j_domain_begin;


    // -------------------------
    // Interpolation of Ex^(d,p)
    // -------------------------
    (*ELoc).x =  compute( coeffxd_, coeffyp_, Ex2D, id_, jp_);

    // -------------------------
    // Interpolation of Ey^(p,d)
    // -------------------------
    (*ELoc).y = compute( coeffxp_, coeffyd_, Ey2D, ip_, jd_);

    // -------------------------
    // Interpolation of Ez^(p,p)
    // -------------------------
    (*ELoc).z = compute( coeffxp_, coeffyp_, Ez2D, ip_, jp_);

    // -------------------------
    // Interpolation of Bx^(p,d)
    // -------------------------
    (*BLoc).x = compute( coeffxp_, coeffyd_, Bx2D, ip_, jd_);

    // -------------------------
    // Interpolation of By^(d,p)
    // -------------------------
    (*BLoc).y = compute( coeffxd_, coeffyp_, By2D, id_, jp_);

    // -------------------------
    // Interpolation of Bz^(d,d)
    // -------------------------
    (*BLoc).z = compute( coeffxd_, coeffyd_, Bz2D, id_, jd_);

} // END Interpolator2D2Order

void Interpolator2D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
{
    // Interpolate E, B
    // Compute coefficient for ipart position
    (*this)(EMfields, particles, ipart, ELoc, BLoc);

    // Static cast of the electromagnetic fields
    Field2D* Jx2D = static_cast<Field2D*>(EMfields->Jx_);
//    Field2D* Jy2D = static_cast<Field2D*>(EMfields->Jy_);
    Field2D* Jz2D = static_cast<Field2D*>(EMfields->Jz_);
    Field2D* Rho2D= static_cast<Field2D*>(EMfields->rho_);
    
    
    // -------------------------
    // Interpolation of Jx^(d,p)
    // -------------------------
    (*JLoc).x = compute( coeffxd_, coeffyp_, Jx2D, id_, jp_);
    
    // -------------------------
    // Interpolation of Jy^(p,d)
    // -------------------------
    (*JLoc).y = compute( coeffxp_, coeffyd_, Jx2D, ip_, jd_);
    
    // -------------------------
    // Interpolation of Ez^(p,p)
    // -------------------------
    (*JLoc).z = compute( coeffxp_, coeffyp_, Jz2D, ip_, jp_);
    
    // -------------------------
    // Interpolation of Rho^(p,p)
    // -------------------------
    (*RhoLoc) = compute( coeffxp_, coeffyp_, Rho2D, ip_, jp_);

}
