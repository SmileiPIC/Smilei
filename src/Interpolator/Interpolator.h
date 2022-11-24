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
    Interpolator() {};
    virtual ~Interpolator() {};
    
    virtual void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) = 0;
    virtual void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) = 0;
    virtual void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) = 0;
    virtual void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) =0;
    
    virtual void fieldsAndEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int , int = 0 )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };
    
    virtual void timeCenteredEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int , int = 0 )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };
    
    virtual void envelopeAndSusceptibility( ElectroMagn *, Particles &, int , double *, double *, double *, double * )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };
 
    virtual void envelopeFieldForIonization( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int , int = 0 )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };
    
private:

};//END class

#endif
