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
    Interpolator( Params &params, Patch *patch );
    virtual ~Interpolator() {};
    
    virtual void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) = 0;
    virtual void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) = 0;
    virtual void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) = 0;
    virtual void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) =0;
    
    virtual void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };
    
    virtual void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };
    
    virtual void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };
 
    virtual void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };
    
private:

};//END class

#endif
