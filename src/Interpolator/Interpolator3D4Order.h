#ifndef INTERPOLATOR3D4ORDER_H
#define INTERPOLATOR3D4ORDER_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 3Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D4Order : public Interpolator3D
{

public:

    //! Creator for Interpolator3D4Order
    Interpolator3D4Order( Params &, Patch * );
    ~Interpolator3D4Order() override {};
    
    //! 4th Order Interpolation of the fields at a the particle position
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override ;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override ;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;

    //! Interpolator on another field than the basic ones
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override;

    //! Interpolator specific to the envelope model
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override ;

    //! Interpolator specific to the envelope model
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override ;

    //! Interpolator specific to the envelope model
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override ;

    //! Computation of a field from provided coefficients
    inline double __attribute__((always_inline)) compute( double *coeffx, double *coeffy, double *coeffz, Field3D *f, int idx, int idy, int idz )
    {
        double interp_res( 0. );
        //unroll ?
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                    interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * *( coeffz+kloc ) * ( *f )( idx+iloc, idy+jloc, idz+kloc );
                }
            }
        }
        return interp_res;
    };

private:

    inline void __attribute__((always_inline)) coeffs( double xpn, double ypn, double zpn, int* idx_p, int* idx_d,
                        double *coeffxp, double *coeffyp, double *coeffzp,
                        double *coeffxd, double *coeffyd, double *coeffzd, double* delta_p )
    {
        // Indexes of the central nodes
        idx_p[0] = round( xpn );
        idx_p[1] = round( ypn );
        idx_p[2] = round( zpn );
        
        // Declaration and calculation of the coefficient for interpolation
        double delta, delta2, delta3, delta4;
        
        delta_p[0] = xpn - ( double )idx_p[0];
        delta2     = delta_p[0]*delta_p[0];
        delta3     = delta2*delta_p[0];
        delta4     = delta3*delta_p[0];
        coeffxp[0] = dble_1_ov_384   - dble_1_ov_48  * delta_p[0]  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffxp[1] = dble_19_ov_96   - dble_11_ov_24 * delta_p[0]  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxp[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffxp[3] = dble_19_ov_96   + dble_11_ov_24 * delta_p[0]  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxp[4] = dble_1_ov_384   + dble_1_ov_48  * delta_p[0]  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        delta_p[1] = ypn - ( double )idx_p[1];
        delta2     = delta_p[1]*delta_p[1];
        delta3     = delta2*delta_p[1];
        delta4     = delta3*delta_p[1];
        coeffyp[0] = dble_1_ov_384   - dble_1_ov_48  * delta_p[1]  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffyp[1] = dble_19_ov_96   - dble_11_ov_24 * delta_p[1]  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyp[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffyp[3] = dble_19_ov_96   + dble_11_ov_24 * delta_p[1]  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyp[4] = dble_1_ov_384   + dble_1_ov_48  * delta_p[1]  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        delta_p[2] = zpn - ( double )idx_p[2];
        delta2     = delta_p[2]*delta_p[2];
        delta3     = delta2*delta_p[2];
        delta4     = delta3*delta_p[2];
        coeffzp[0] = dble_1_ov_384   - dble_1_ov_48  * delta_p[2]  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffzp[1] = dble_19_ov_96   - dble_11_ov_24 * delta_p[2]  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzp[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffzp[3] = dble_19_ov_96   + dble_11_ov_24 * delta_p[2]  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzp[4] = dble_1_ov_384   + dble_1_ov_48  * delta_p[2]  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        // First index for summation
        idx_p[0]   = idx_p[0] - i_domain_begin;
        idx_p[1]   = idx_p[1] - j_domain_begin;
        idx_p[2]   = idx_p[2] - k_domain_begin;

        if(idx_d){
            idx_d[0] = round( xpn+0.5 );
            idx_d[1] = round( ypn+0.5 );
            idx_d[2] = round( zpn+0.5 );


            delta    = xpn - ( double )idx_d[0] + 0.5;
            delta2     = delta*delta;
            delta3     = delta2*delta;
            delta4     = delta3*delta;
            coeffxd[0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeffxd[1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeffxd[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeffxd[3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeffxd[4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

            delta    = ypn - ( double )idx_d[1] + 0.5;
            delta2     = delta*delta;
            delta3     = delta2*delta;
            delta4     = delta3*delta;
            coeffyd[0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeffyd[1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeffyd[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeffyd[3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeffyd[4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

            delta    = zpn - ( double )idx_d[2] + 0.5;
            delta2     = delta*delta;
            delta3     = delta2*delta;
            delta4     = delta3*delta;
            coeffzd[0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeffzd[1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeffzd[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeffzd[3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeffzd[4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

            idx_d[0]   = idx_d[0] - i_domain_begin;
            idx_d[1]   = idx_d[1] - j_domain_begin;
            idx_d[2]   = idx_d[2] - k_domain_begin;
        }

    }
    
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
