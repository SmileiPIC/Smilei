#ifndef INTERPOLATOR3D_H
#define INTERPOLATOR3D_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 3D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D : public Interpolator
{
public:

    //! Creator for Interpolator3D
    Interpolator3D( Params &params, Patch *patch );

    virtual ~Interpolator3D() override {} ;

    //! Interpolation of all fields and currents for a single particles located at istart.
    virtual void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override = 0;

    virtual void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override = 0  ;

    virtual void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override = 0;

    virtual void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override = 0;

protected:
    //! Inverse of the spatial-step
    double dx_inv_;
    double dy_inv_;
    double dz_inv_;
    int i_domain_begin;
    int j_domain_begin;
    int k_domain_begin;
    double d_inv_[3];
};

#endif
