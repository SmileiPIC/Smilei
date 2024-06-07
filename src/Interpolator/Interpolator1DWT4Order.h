#ifndef INTERPOLATOR1DWT4ORDER_H
#define INTERPOLATOR1DWT4ORDER_H


#include "Interpolator1D.h"
#include "Field1D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 4th order WT interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1DWT4Order final : public Interpolator1D
{

public:
    Interpolator1DWT4Order( Params &, Patch * );
    ~Interpolator1DWT4Order() override final {};

    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override final;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;

    inline double __attribute__((always_inline)) compute( double *coeff, Field1D *f, int idx )
    {
        double interp_res =  coeff[0] * ( *f )( idx-2 )   + coeff[1] * ( *f )( idx-1 )   + coeff[2] * ( *f )( idx ) + coeff[3] * ( *f )( idx+1 ) + coeff[4] * ( *f )( idx+2 );
        return interp_res;
    };

    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;

private:
    inline void __attribute__((always_inline)) coeffs( double xjn )
    {
        double var1, var2, var3, var4, var5;

        // Dual
        id_      = round( xjn+0.5 ); // index of the central point
        xjmxi  = xjn -( double )id_+0.5; // normalized distance to the central node
        var1 = dble_11_ov_24 * xjmxi;
        var2 = xjmxi * xjmxi;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - xjmxi;
        var4 = var4 * var4;
        var5 = 0.5 + xjmxi;
        var5 = var5 * var5;

        // coefficients for the 4th order interpolation on 5 nodes
        coeffd_[0] = dble_1_ov_24 * var4 * var4;
        coeffd_[1] = dble_19_ov_96 - var1 + var3 * ( 1.5+xjmxi-var2 );
        coeffd_[2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeffd_[3] = dble_19_ov_96 + var1 + var3 * ( 1.5-xjmxi -var2 );
        coeffd_[4] = dble_1_ov_24 * var5 * var5;

        id_ -= i_domain_begin_;

        // Primal
        ip_      = round( xjn );    // index of the central point
        xjmxi  = xjn -( double )ip_; // normalized distance to the central node
        var1 = dble_11_ov_24 * xjmxi;
        var2 = xjmxi * xjmxi;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - xjmxi;
        var4 = var4 * var4;
        var5 = 0.5 + xjmxi;
        var5 = var5 * var5;

        // coefficients for the 4th order interpolation on 5 nodes
        coeffp_[0] = dble_1_ov_24 * var4 * var4;
        coeffp_[1] = dble_19_ov_96 - var1 + var3 * ( 1.5+xjmxi-var2 );
        coeffp_[2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeffp_[3] = dble_19_ov_96 + var1 + var3 * ( 1.5-xjmxi-var2 );
        coeffp_[4] = dble_1_ov_24 * var5 * var5;

        // Primal
        var1 = dt_ov_dx - xjmxi;
        var1 = var1 * var1;
        var1 = dx_ov_96dt * var1 * var1;
        var2 = dt_ov_dx + xjmxi;
        var2 = var2 * var2;
        var2 = dx_ov_96dt * var2 * var2;
        var3 = copysign( var1, dt_ov_dx-xjmxi ) + copysign( var2, dt_ov_dx+xjmxi);
        var4 =  dble_1_ov_24 * ((( 3.0+xjmxi ) * xjmxi  - ( 3.0-dt2_ov_dx2 )) * xjmxi + (1.0+dt2_ov_dx2 ));
        var5 =  dble_1_ov_24 * ((( 3.0-xjmxi ) * xjmxi  + ( 3.0-dt2_ov_dx2 )) * xjmxi + (1.0+dt2_ov_dx2 ));

        // coefficients for the 4th order WT interpolation on 5 nodes
        coeffpt_[0] = var3 + var1 - var2;
        coeffpt_[1] = 4.0 * ( var4-var3 );
        coeffpt_[2] = 1.0 + 6.0 * var3 - 4.0 * ( var4+var5 );
        coeffpt_[3] = 4.0 * ( var5-var3 );
        coeffpt_[4] = var3 + var2 - var1;
        
        
        ip_ -= i_domain_begin_;
    }

    double dble_1_ov_6 ;
    double dble_1_ov_24 ;
    double dble_11_ov_24 ;
    double dble_19_ov_96 ;
    double dble_115_ov_192 ;
    double dt_ov_dx ;
    double dt2_ov_dx2 ;
    double dx_ov_96dt ;
    
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Last delta computed
    double xjmxi;
    // Interpolation coefficient on Prim grid
    double coeffp_[5];
    // WT interpolation coefficient on Prim grid
    double coeffpt_[5];
    // Interpolation coefficient on Dual grid
    double coeffd_[5];


};//END class

#endif
