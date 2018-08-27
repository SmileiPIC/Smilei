#ifndef INTERPOLATORRZ_H
#define INTERPOLATORRZ_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator RZ
//  --------------------------------------------------------------------------------------------------------------------
class InterpolatorRZ : public Interpolator
{
public:
    InterpolatorRZ(Params& params, Patch *patch);

    virtual ~InterpolatorRZ() override {} ;

    virtual void operator()  (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread) override = 0  ;
    virtual void operator()  (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc) override = 0;
    virtual void operator()  (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override = 0;

protected:
    //! Inverse of the spatial-step
    double dl_inv_;
    double dr_inv_;
    int i_domain_begin;
    int j_domain_begin;
    double dr;
};

#endif
