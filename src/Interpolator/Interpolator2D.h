#ifndef INTERPOLATOR2D_H
#define INTERPOLATOR2D_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 2D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D : public Interpolator
{
public:
    Interpolator2D(Params& params, Patch *patch);

    virtual ~Interpolator2D() override {} ;

    virtual void fieldsAndCurrents(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override = 0;
    virtual void fields_batch     (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override = 0  ;
    virtual void fields_selection (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override = 0;
    virtual void oneField         (Field* field, Particles &particles, int *istart, int *iend, double* FieldLoc) override = 0;

protected:
    //! Inverse of the spatial-step
    double dx_inv_;
    double dy_inv_;
    int i_domain_begin;
    int j_domain_begin;
    double D_inv[2];
};
#endif
