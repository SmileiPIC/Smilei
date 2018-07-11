#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Field.h"
#include "Params.h"

class Params;
class Patch;
class ElectroMagn;
class Particles;


//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator
{
public:
    Interpolator(Params& params, Patch* patch);
    virtual ~Interpolator() {};
    
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) = 0;
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) = 0;
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) = 0;

    virtual void interpolate_em_fields_and_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };

    virtual void interpolate_envelope_and_old_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };


private:

};//END class

#endif
