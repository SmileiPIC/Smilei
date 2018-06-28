#ifndef INTERPOLATOR2D2ORDERV_H
#define INTERPOLATOR2D2ORDERV_H


#include "Interpolator2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D2OrderV : public Interpolator2D
{

public:
    Interpolator2D2OrderV(Params&, Patch*);
    ~Interpolator2D2OrderV() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double* ELoc, double* BLoc);
    // Sorting
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final;
    // Probes
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final ;

    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final{};

    inline double compute( double* coeffx, double* coeffy, Field2D* f, int idx, int idy) {
        double interp_res(0.);
        //unroll ?
        for (int iloc=-1 ; iloc<2 ; iloc++) {
            for (int jloc=-1 ; jloc<2 ; jloc++) {
                interp_res += *(coeffx+iloc) * *(coeffy+jloc) * (*f)(idx+iloc,idy+jloc);
            }
        }
        return interp_res;
    };  

private:


};//END class

#endif
