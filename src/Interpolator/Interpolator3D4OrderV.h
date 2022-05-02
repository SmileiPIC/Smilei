#ifndef INTERPOLATOR3D4ORDERV_H
#define INTERPOLATOR3D4ORDERV_H


#include "Interpolator3D4Order.h"
#include "Field3D.h"
#include "Pragma.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class for vectorized 4th order interpolator for 3d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D4OrderV final : public Interpolator3D4Order
{

public:

    //! Creator for Interpolator3D4Order
    Interpolator3D4OrderV( Params &, Patch * );

    //! Destructor for Interpolator3D4Order
    ~Interpolator3D4OrderV() override final {};

    //inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, double *ELoc, double *BLoc );
    // inline void fields( ElectroMagn *EMfields, Particles &particles, int ipart, double *ELoc, double *BLoc );

    //! Interpolation of all fields and currents for a single particles located at istart.
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;

    //! Wrapper called by the particle dynamics section
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override final;

    // void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final {};

    //! Interpolator on another field than the basic ones
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;

    //! Interpolator specific to the envelope model
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    //! Interpolator specific to the envelope model
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    //! Interpolator specific to the envelope model
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;

#if defined __INTEL_COMPILER

    __declspec(noinline) void interp_Bx( int* idxO, int np_computed, double *coeffxp, double *coeffyd, double *coeffzd, int *dualy, int* dualz, Field3D *Bx3D, double *Bpart );
    __declspec(noinline) void interp_By( int* idxO, int np_computed, double *coeffxd, double *coeffyp, double *coeffzd, int *dualx, int* dualz, Field3D *By3D, double *Bpart );
    __declspec(noinline) void interp_Bz( int* idxO, int np_computed, double *coeffxd, double *coeffyd, double *coeffzp, int *dualx, int* dualy, Field3D *Bz3D, double *Bpart );

#else

    void interp_Bx( int* idxO, int np_computed, double *coeffxp, double *coeffyd, double *coeffzd, int *dualy, int* dualz, Field3D *Bx3D, double *Bpart );
    void interp_By( int* idxO, int np_computed, double *coeffxd, double *coeffyp, double *coeffzd, int *dualx, int* dualz, Field3D *By3D, double *Bpart );
    void interp_Bz( int* idxO, int np_computed, double *coeffxd, double *coeffyd, double *coeffzp, int *dualx, int* dualy, Field3D *Bz3D, double *Bpart );

#endif

private:

    double dble_1_ov_384 ;
    double dble_1_ov_48 ;
    double dble_1_ov_16 ;
    double dble_1_ov_12 ;
    double dble_1_ov_24 ;
    double dble_19_ov_96 ;
    double dble_11_ov_24 ;
    double dble_1_ov_4 ;
    double dble_1_ov_6 ;
    double dble_115_ov_192 ;
    double dble_5_ov_8 ;


};//END class

#endif
