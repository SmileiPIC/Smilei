#include "InterpolatorRZ2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "ElectroMagn3DRZ.h"
#include "cField2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for InterpolatorRZ2Order
// ---------------------------------------------------------------------------------------------------------------------
InterpolatorRZ2Order::InterpolatorRZ2Order(Params &params, Patch *patch) : InterpolatorRZ(params, patch)
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void InterpolatorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
{
    #ifdef _TODO_RZ
    // Loop on modes ?
    #endif

    int imode = 0;
    // Static cast of the electromagnetic fields
    cField2D* ExRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->El_[imode];
    cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Er_[imode];
    cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Et_[imode];
    cField2D* BxRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bl_m[imode];
    cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Br_m[imode];
    cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bt_m[imode];


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
    (*ELoc).x = compute( &coeffxd_[1], &coeffyp_[1], ExRZ, id_, jp_);

    // -------------------------
    // Interpolation of Er^(p,d)
    // -------------------------
    (*ELoc).y = compute( &coeffxp_[1], &coeffyd_[1], ErRZ, ip_, jd_);

    // -------------------------
    // Interpolation of Et^(p,p)
    // -------------------------
    (*ELoc).z = compute( &coeffxp_[1], &coeffyp_[1], EtRZ, ip_, jp_);

    // -------------------------
    // Interpolation of Bx^(p,d)
    // -------------------------
    (*BLoc).x = compute( &coeffxp_[1], &coeffyd_[1], BxRZ, ip_, jd_);

    // -------------------------
    // Interpolation of Br^(d,p)
    // -------------------------
    (*BLoc).y = compute( &coeffxd_[1], &coeffyp_[1], BrRZ, id_, jp_);

    // -------------------------
    // Interpolation of Bt^(d,d)
    // -------------------------
    (*BLoc).z = compute( &coeffxd_[1], &coeffyd_[1], BtRZ, id_, jd_);

} // END InterpolatorRZ2Order

void InterpolatorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
{
    #ifdef _TODO_RZ
    // Loop on modes ?
    #endif

    // Interpolate E, B
    // Compute coefficient for ipart position
    (*this)(EMfields, particles, ipart, ELoc, BLoc);

    int imode = 0;
    // Static cast of the electromagnetic fields
    cField2D* JxRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jl_[imode]; 
    cField2D* JyRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jr_[imode]; 
    cField2D* JzRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jt_[imode]; 
    cField2D* RhoRZ= (static_cast<ElectroMagn3DRZ*>(EMfields))->rho_RZ_[imode];
    
    
    // -------------------------
    // Interpolation of Jx^(d,p,p)
    // -------------------------
    (*JLoc).x = compute( &coeffxd_[1], &coeffyp_[1], JxRZ, id_, jp_);
    
    // -------------------------
    // Interpolation of Jy^(p,d,p)
    // -------------------------
    (*JLoc).y = compute( &coeffxp_[1], &coeffyd_[1], JyRZ, ip_, jd_);
    
    // -------------------------
    // Interpolation of Jz^(p,p,d)
    // -------------------------
    (*JLoc).z = compute( &coeffxp_[1], &coeffyp_[1], JzRZ, ip_, jp_);
    
    // -------------------------
    // Interpolation of Rho^(p,p,p)
    // -------------------------
    (*RhoLoc) = compute( &coeffxp_[1], &coeffyp_[1], RhoRZ, ip_, jp_);

}

void InterpolatorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread)
{
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);

    //Loop on bin particles
    for (int ipart=istart ; ipart<iend; ipart++ ) {
        //Interpolation on current particle
        (*this)(EMfields, particles, ipart, &(*Epart)[ipart], &(*Bpart)[ipart]);
        //Buffering of iol and delta
        (*iold)[ipart*3]    = ip_;
        (*iold)[ipart*3+1]  = jp_;
        (*delta)[ipart*3]   = deltax;
        (*delta)[ipart*3+1] = deltay;
    }

}
