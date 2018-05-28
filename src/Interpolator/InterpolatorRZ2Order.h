#ifndef INTERPOLATORRZ2ORDER_H
#define INTERPOLATORRZ2ORDER_H


#include "InterpolatorRZ.h"
#include "cField2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class InterpolatorRZ2Order : public InterpolatorRZ
{

public:
    InterpolatorRZ2Order(Params&, Patch*);
    ~InterpolatorRZ2Order() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final;

    inline double compute( double* coeffx, double* coeffy, cField2D* f, int idx, int idy) {
	double interp_res(0.);
        //unroll ?
	for (int iloc=-1 ; iloc<2 ; iloc++) {
	    for (int jloc=-1 ; jloc<2 ; jloc++) {
                #ifdef _TODO_RZ
                #endif
                interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ( real( (*f)(idx+iloc,idy+jloc) ) + imag( (*f)(idx+iloc,idy+jloc) ) );
	    }
	}
	return interp_res;
    };  

private:
    // Last prim index computed
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Last delta computed
    double deltax, deltay;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3];
    //! Number of modes;
    unsigned int nmodes;


};//END class

#endif
