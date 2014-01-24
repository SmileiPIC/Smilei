#include "Interpolator2D4Order.h"
#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particle.h"
#include "SmileiMPI_Cart2D.h"

#include <iostream>
#include <cmath>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D4Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D4Order::Interpolator2D4Order(PicParams *params, SmileiMPI* smpi) : Interpolator2D(params, smpi)
{
	SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);

	dx_inv_ = 1.0/params->cell_length[0];
    dy_inv_ = 1.0/params->cell_length[1];
    
    
	i_domain_begin = smpi2D->getCellStartingGlobalIndex(0);
    j_domain_begin = smpi2D->getCellStartingGlobalIndex(1);
    
    //double defined for use in coefficients
    dble_1_ov_384 = 1.0/384.0;
    dble_1_ov_48 = 1.0/48.0;
    dble_1_ov_16 = 1.0/16.0;
    dble_1_ov_12 = 1.0/12.0;
    dble_1_ov_24 = 1.0/24.0;
    dble_19_ov_96 = 19.0/96.0;
    dble_11_ov_24 = 11.0/24.0;
    dble_1_ov_4 = 1.0/4.0;
    dble_1_ov_6 = 1.0/6.0;
    dble_115_ov_192 = 115.0/192.0;
    dble_5_ov_8 = 5.0/8.0;


}

Interpolator2D4Order::~Interpolator2D4Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D4Order::operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc)
{
    // Static cast of the electromagnetic fields
	Field2D* Ex2D = static_cast<Field2D*>(champs->Ex_);
	Field2D* Ey2D = static_cast<Field2D*>(champs->Ey_);
	Field2D* Ez2D = static_cast<Field2D*>(champs->Ez_);
	Field2D* Bx2D = static_cast<Field2D*>(champs->Bx_m);
	Field2D* By2D = static_cast<Field2D*>(champs->By_m);
	Field2D* Bz2D = static_cast<Field2D*>(champs->Bz_m);
    
    
    // Normalized particle position
    double xpn = part->position(0)*dx_inv_;
    double ypn = part->position(1)*dy_inv_;
    
    
    // Indexes of the central nodes
    int ic_p = round(xpn);
    int ic_d = round(xpn+0.5);
    int jc_p = round(ypn);
    int jc_d = round(ypn+0.5);
    
    
    // Declaration and calculation of the coefficient for interpolation
    double delta, delta2, delta3, delta4;
    
    std::vector<double> Cx_p(5);
    delta   = xpn - (double)ic_p;
    delta2  = delta*delta;
    delta3  = delta2*delta;
    delta4  = delta3*delta;
    Cx_p[0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_12 * delta4;
    Cx_p[1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
	Cx_p[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
    Cx_p[3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Cx_p[4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_12 * delta4;
    
    std::vector<double> Cx_d(5);
    delta   = xpn - (double)ic_d + 0.5;
    delta2  = delta*delta;
    delta3  = delta2*delta;
    delta4  = delta3*delta;
    Cx_d[0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_12 * delta4;
    Cx_d[1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
	Cx_d[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
    Cx_d[3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Cx_d[4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_12 * delta4;
    
    std::vector<double> Cy_p(5);
    delta   = ypn - (double)jc_p;
    delta2  = delta*delta;
    delta3  = delta2*delta;
    delta4  = delta3*delta;
    Cy_p[0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_12 * delta4;
    Cy_p[1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
	Cy_p[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
    Cy_p[3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Cy_p[4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_12 * delta4;
    
    
    std::vector<double> Cy_d(5);
    delta   = ypn - (double)jc_d + 0.5;
    delta2  = delta*delta;
    delta2  = delta*delta;
    delta2  = delta*delta;
    delta3  = delta2*delta;
    delta4  = delta3*delta;
    Cy_d[0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_12 * delta4;
    Cy_d[1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
	Cy_d[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
    Cy_d[3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Cy_d[4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_12 * delta4;
    
    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    int ip = ic_p - 2 - i_domain_begin;
    int id = ic_d - 2 - i_domain_begin;
    int jp = jc_p - 2 - j_domain_begin;
    int jd = jc_d - 2 - j_domain_begin;
    
    
    // -------------------------
	// Interpolation of Ex^(d,p)
    // -------------------------
    (*ELoc).x = 0.0;
    for (int iloc=0 ; iloc<5 ; iloc++) {
        for (int jloc=0 ; jloc<5 ; jloc++) {
            (*ELoc).x += Cx_d[iloc] * Cy_p[jloc] * (*Ex2D)(id+iloc,jp+jloc);
        }
    }
    
    // -------------------------
	// Interpolation of Ey^(p,d)
    // -------------------------
    (*ELoc).y = 0.0;
    for (int iloc=0 ; iloc<5 ; iloc++) {
        for (int jloc=0 ; jloc<5 ; jloc++) {
            (*ELoc).y += Cx_p[iloc] * Cy_d[jloc] * (*Ey2D)(ip+iloc,jd+jloc);
        }
    }
    
    // -------------------------
	// Interpolation of Ez^(p,p)
    // -------------------------
    (*ELoc).z = 0.0;
    for (int iloc=0 ; iloc<5 ; iloc++) {
        for (int jloc=0 ; jloc<5 ; jloc++) {
            (*ELoc).z += Cx_p[iloc] * Cy_p[jloc] * (*Ez2D)(ip+iloc,jp+jloc);
        }
    }

    // -------------------------
	// Interpolation of Bx^(p,d)
    // -------------------------
    (*BLoc).x = 0.0;
    for (int iloc=0 ; iloc<5 ; iloc++) {
        for (int jloc=0 ; jloc<5 ; jloc++) {
            (*BLoc).x += Cx_p[iloc] * Cy_d[jloc] * (*Bx2D)(ip+iloc,jd+jloc);
        }
    }
    
    // -------------------------
	// Interpolation of By^(d,p)
    // -------------------------
    (*BLoc).y = 0.0;
    for (int iloc=0 ; iloc<5 ; iloc++) {
        for (int jloc=0 ; jloc<5 ; jloc++) {
            (*BLoc).y += Cx_d[iloc] * Cy_p[jloc] * (*By2D)(id+iloc,jp+jloc);
        }
    }
    
    // -------------------------
	// Interpolation of Bz^(d,d)
    // -------------------------
    (*BLoc).z = 0.0;
    for (int iloc=0 ; iloc<5 ; iloc++) {
        for (int jloc=0 ; jloc<5 ; jloc++) {
            (*BLoc).z += Cx_d[iloc] * Cy_d[jloc] * (*Bz2D)(id+iloc,jd+jloc);
        }
    }
    
} // END Interpolator2D4Order
