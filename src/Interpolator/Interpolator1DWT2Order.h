#ifndef INTERPOLATOR1DWT2ORDER_H
#define INTERPOLATOR1DWT2ORDER_H


#include "Interpolator1D.h"

#include "Field1D.h"
//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order WT interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1DWT2Order final : public Interpolator1D
{

public:
    Interpolator1DWT2Order( Params &, Patch * );
    ~Interpolator1DWT2Order() override final {};

    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override final;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;

    inline double __attribute__((always_inline)) compute( double *coeff, Field1D *f, int idx )
    {
        double interp_res =  coeff[0] * ( *f )( idx-1 )   + coeff[1] * ( *f )( idx )   + coeff[2] * ( *f )( idx+1 );
        return interp_res;
    };

    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;
    void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

private:
    inline void __attribute__((always_inline)) coeffs( double xjn )
    {
        double var1;

        // Dual
        id_      = round( xjn+0.5 );      // index of the central point
        xjmxi  = xjn - ( double )id_ +0.5; // normalized distance to the central node
        var1 = xjmxi*xjmxi;            // square of the normalized distance to the central node

        // 2nd order interpolation on 3 nodes
        coeffd_[0] = 0.5 * ( var1-xjmxi+0.25 );
        coeffd_[1] = ( 0.75-var1 );
        coeffd_[2] = 0.5 * ( var1+xjmxi+0.25 );

        id_ -= i_domain_begin_;

        // Primal
        ip_      = round( xjn );    // index of the central point
        xjmxi  = xjn -( double )ip_; // normalized distance to the central node
        var1 = xjmxi * xjmxi;   // square of the normalized distance to the central node

        // 2nd order interpolation on 3 nodes
        coeffp_[0] = 0.5 * ( var1-xjmxi+0.25 );
        coeffp_[1] = ( 0.75-var1 );
        coeffp_[2] = 0.5 * ( var1+xjmxi+0.25 );

        // 2nd order WT interpolation on 3 nodes
        var1 = dx_ov_8dt * ( fabs( dt_ov_dx-xjmxi ) * ( dt_ov_dx-xjmxi ) + fabs( dt_ov_dx+xjmxi ) * ( dt_ov_dx+xjmxi ) );
        coeffpt_[0] = var1 - 0.5 * xjmxi;
        coeffpt_[1] = 1.0 - 2.0 * var1;
        coeffpt_[2] = var1 + 0.5 * xjmxi;
        
        ip_ -= i_domain_begin_;
    }

    // Coefficients for WT
    double dt_ov_dx;
    double dx_ov_8dt;
    
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Last delta computed
    double xjmxi, deltax;
    // Interpolation coefficient on Prim grid
    double coeffp_[3];
    // WT interpolation coefficient on Prim grid
    double coeffpt_[3];
    // Interpolation coefficient on Dual grid
    double coeffd_[3];


};//END class

#endif
