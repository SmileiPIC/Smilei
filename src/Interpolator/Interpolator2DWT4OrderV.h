#ifndef INTERPOLATOR2DWT4ORDERV_H
#define INTERPOLATOR2DWT4ORDERV_H


#include "Interpolator2DWT4Order.h"
#include "Field2D.h"
#include "Pragma.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for vectorized 4th order WT interpolator for 2Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2DWT4OrderV final : public Interpolator2DWT4Order
{

public:

    //! Creator for Interpolator2DWT4OrderV
    Interpolator2DWT4OrderV( Params &, Patch * );

    //! Destructor for Interpolator2DWT4OrderV
    ~Interpolator2DWT4OrderV() override final {};

    //inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );

    //! Interpolation of all fields and currents for a single particles located at istart.
    // void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;

    //void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final {};

    //! Wrapper called by the particle dynamics section
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell, int ipart_ref = 0 ) override final ;

    //! Interpolator on another field than the basic ones
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;

    //! Interpolator specific to the envelope model
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    //! Interpolator specific to the envelope model
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    //! Interpolator specific to the envelope model
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;

private:

    double dble_1_ov_6 ;
    double dble_1_ov_24 ;
    double dble_11_ov_24 ;
    double dble_19_ov_96 ;
    double dble_115_ov_192 ;
    double dt_ov_D[2] ;
    double dt2_ov_D2[2] ;
    double D_ov_96dt[2] ;

};//END class

#endif
