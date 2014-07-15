#include "Interpolator1D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "SmileiMPI_Cart1D.h"

using namespace std;

Interpolator1D2Order::Interpolator1D2Order(PicParams *params, SmileiMPI* smpi) : Interpolator1D(params, smpi)
{
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

    dx_inv_ = 1.0/params->cell_length[0];

    index_domain_begin = smpi1D->getCellStartingGlobalIndex(0);

}

Interpolator1D2Order::~Interpolator1D2Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
{

    // Variable declaration
    int i;
    double xjn, xjmxi, xjmxi2;
    double coeffInf, coeffCur, coeffSup;

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
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    // --------------------------------------------------------
    i      = round(xjn);      // index of the central point
    xjmxi  = xjn -(double)i;  // normalized distance to the central node
    xjmxi2 = pow(xjmxi,2);    // square of the normalized distance to the central node

    // 2nd order interpolation on 3 nodes
    coeffInf = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffCur = (0.75-xjmxi2);
    coeffSup = 0.5 * (xjmxi2+xjmxi+0.25);

    i -= index_domain_begin;

    (*ELoc).y =  coeffInf * (*Ey1D)(i-1)   + coeffCur * (*Ey1D)(i)   + coeffSup * (*Ey1D)(i+1);
    (*ELoc).z =  coeffInf * (*Ez1D)(i-1)   + coeffCur * (*Ez1D)(i)   + coeffSup * (*Ez1D)(i+1);
    (*BLoc).x =  coeffInf * (*Bx1D_m)(i-1) + coeffCur * (*Bx1D_m)(i) + coeffSup * (*Bx1D_m)(i+1);


    // --------------------------------------------------------
    // Interpolate the fields from the Dual grid : Ex, By, Bz
    // --------------------------------------------------------
    i      = round(xjn+0.5);        // index of the central point
    xjmxi  = xjn - (double)i +0.5;  // normalized distance to the central node
    xjmxi2 = pow(xjmxi,2);          // square of the normalized distance to the central node

    // 2nd order interpolation on 3 nodes
    coeffInf = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffCur = (0.75-xjmxi2);
    coeffSup = 0.5 * (xjmxi2+xjmxi+0.25);

    i -= index_domain_begin;

    (*ELoc).x =  coeffInf * (*Ex1D)(i-1)   + coeffCur * (*Ex1D)(i)   + coeffSup * (*Ex1D)(i+1);
    (*BLoc).y =  coeffInf * (*By1D_m)(i-1) + coeffCur * (*By1D_m)(i) + coeffSup * (*By1D_m)(i+1);
    (*BLoc).z =  coeffInf * (*Bz1D_m)(i-1) + coeffCur * (*Bz1D_m)(i) + coeffSup * (*Bz1D_m)(i+1);

}//END Interpolator1D2Order

void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc){
    (*this)(EMfields, particles, ipart, ELoc, BLoc);
    
    // Variable declaration
    int i;
    double xjn, xjmxi, xjmxi2;
    double coeffInf, coeffCur, coeffSup;
    
    // Static cast of the electromagnetic fields
    Field1D* Jx1D     = static_cast<Field1D*>(EMfields->Jx_);
    Field1D* Jy1D     = static_cast<Field1D*>(EMfields->Jy_);
    Field1D* Jz1D     = static_cast<Field1D*>(EMfields->Jz_);
    Field1D* Rho1D     = static_cast<Field1D*>(EMfields->rho_);
    
    
    // Particle position (in units of the spatial-step)
    xjn    = particles.position(0, ipart)*dx_inv_;
    
    
    // --------------------------------------------------------
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    // --------------------------------------------------------
    i      = round(xjn);      // index of the central point
    xjmxi  = xjn -(double)i;  // normalized distance to the central node
    xjmxi2 = pow(xjmxi,2);    // square of the normalized distance to the central node
    
    // 2nd order interpolation on 3 nodes
    coeffInf = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffCur = (0.75-xjmxi2);
    coeffSup = 0.5 * (xjmxi2+xjmxi+0.25);
    
    i -= index_domain_begin;

    (*JLoc).y =  coeffInf * (*Jy1D)(i-1)   + coeffCur * (*Jy1D)(i)   + coeffSup * (*Jy1D)(i+1);
    (*JLoc).z =  coeffInf * (*Jz1D)(i-1)   + coeffCur * (*Jz1D)(i)   + coeffSup * (*Jz1D)(i+1);

    
    // --------------------------------------------------------
    // Interpolate the fields from the Dual grid : Ex, By, Bz
    // --------------------------------------------------------
    i      = round(xjn+0.5);        // index of the central point
    xjmxi  = xjn - (double)i +0.5;  // normalized distance to the central node
    xjmxi2 = pow(xjmxi,2);          // square of the normalized distance to the central node
    
    // 2nd order interpolation on 3 nodes
    coeffInf = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffCur = (0.75-xjmxi2);
    coeffSup = 0.5 * (xjmxi2+xjmxi+0.25);
    
    i -= index_domain_begin;
        
    (*JLoc).x =  coeffInf * (*Jx1D)(i-1)   + coeffCur * (*Jx1D)(i)   + coeffSup * (*Jx1D)(i+1);
    (*RhoLoc) =  coeffInf * (*Rho1D)(i-1)   + coeffCur * (*Rho1D)(i)   + coeffSup * (*Rho1D)(i+1);
    
}
