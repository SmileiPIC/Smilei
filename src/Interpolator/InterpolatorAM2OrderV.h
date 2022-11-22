#ifndef INTERPOLATORAM2ORDERV_H
#define INTERPOLATORAM2ORDERV_H


#include "InterpolatorAM.h"
#include "cField2D.h"
#include "Field2D.h"
#include "Pragma.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for AM spectral simulations
//  --------------------------------------------------------------------------------------------------------------------
class InterpolatorAM2OrderV final : public InterpolatorAM
{

public:
    InterpolatorAM2OrderV( Params &, Patch * );
    ~InterpolatorAM2OrderV() override final {};
    
    void fieldsAndCurrents( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, LocalFields *, double * ) override final {};
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell, int ipart_ref = 0 ) override final ;
    void fieldsSelection( ElectroMagn *, Particles &, double *, int, std::vector<unsigned int> * ) override final {};
    void oneField( Field **, Particles &, int *, int *, double *, double * = NULL, double * = NULL, double * = NULL ) override final {};
    
    void fieldsAndEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int = 0 ) override final {};
    void timeCenteredEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int = 0 ) override final {};
    void envelopeAndSusceptibility( ElectroMagn *, Particles &, int , double *, double *, double *, double * ) override final {};
    void envelopeFieldForIonization( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int , int = 0 ) override final {};    


private:
    
    //! Number of modes;
    unsigned int nmodes_;
    
    
};//END class

#endif
