#ifndef INTERPOLATORAM_H
#define INTERPOLATORAM_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator AM
//  --------------------------------------------------------------------------------------------------------------------
class InterpolatorAM : public Interpolator
{
public:
    InterpolatorAM( Params &params, Patch *patch );
    
    virtual ~InterpolatorAM() override {} ;
    
    virtual void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override = 0 ;
    virtual void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override = 0  ;
    virtual void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override = 0;
    virtual void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override = 0;
    
protected:
    //! Inverse of the spatial-step
    int i_domain_begin_;
    int j_domain_begin_;
    double D_inv_[2];
};

#endif
