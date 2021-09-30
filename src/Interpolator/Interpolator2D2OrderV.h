#ifndef INTERPOLATOR2D2ORDERV_H
#define INTERPOLATOR2D2ORDERV_H


#include "Interpolator2D.h"
#include "Interpolator2D2Order.h"
#include "Field2D.h"
#include "Pragma.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 2d3v simulations
//! This class inherits from Interpolator2D2Order for some scalar versions
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D2OrderV final : public Interpolator2D2Order
{

public:

    //! Creator for Interpolator2D2Order
    Interpolator2D2OrderV( Params &, Patch * );

    //! Destructor for Interpolator2D2Order
    ~Interpolator2D2OrderV() override final {};

    // Scalar basic all field interpolation
    // inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, double *ELoc, double *BLoc );

    //! Interpolation of all fields and currents for a single particles located at istart.
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;

    //! Wrapper called by the particle dynamics section
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    // void fieldsSelection(   ElectroMagn *EMfields,
    //                         Particles &particles,
    //                         double *buffer,
    //                         int offset,
    //                         std::vector<unsigned int> *selection );

    //! Interpolator on another field than the basic ones
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;

    //! Interpolator specific to the envelope model
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    //! Interpolator specific to the envelope model
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    //! Interpolator specific to the envelope model
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;

private:


};//END class

#endif
