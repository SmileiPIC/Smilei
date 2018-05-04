#ifndef INTERPOLATOR3D4ORDERV_H
#define INTERPOLATOR3D4ORDERV_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D4OrderV : public Interpolator3D
{

public:
    Interpolator3D4OrderV(Params&, Patch*);
    ~Interpolator3D4OrderV() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double* ELoc, double* BLoc);
    // Sorting
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final;
    // Probes
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final ;

    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final{};

private:
    double dble_1_ov_384 ;
    double dble_1_ov_48 ;
    double dble_1_ov_16 ;
    double dble_1_ov_12 ;
    double dble_1_ov_24 ;
    double dble_19_ov_96 ;
    double dble_11_ov_24 ;
    double dble_1_ov_4 ;
    double dble_1_ov_6 ;
    double dble_115_ov_192 ;
    double dble_5_ov_8 ;


};//END class

#endif
