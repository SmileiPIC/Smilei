#ifndef INTERPOLATOR2D_H
#define INTERPOLATOR2D_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 2D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D : public Interpolator
{
public:
    Interpolator2D( Params &params, Patch *patch );

    virtual ~Interpolator2D() override {} ;

    virtual void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override = 0;
    virtual void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override = 0  ;
    virtual void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override = 0;
    virtual void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override = 0;

protected:

    //! Inverse of the spatial-step
    double d_inv_[2];

    int i_domain_begin;
    int j_domain_begin;

};
#endif
