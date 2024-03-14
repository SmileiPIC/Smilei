#ifndef INTERPOLATOR2D4ORDER_H
#define INTERPOLATOR2D4ORDER_H


#include "Interpolator2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 4th order interpolator for 2Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D4Order : public Interpolator2D
{

public:

    //! Creator for Interpolator2D2Order
    Interpolator2D4Order( Params &, Patch * );
    ~Interpolator2D4Order() override {};
    
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override ;

    //! Wrapper called by the particle dynamics section
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override ;

    //! Interpolator specific to tracked particles. A selection of particles may be provided
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;

    //! Interpolator on another field than the basic ones
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override;

    //! Computation of a field from provided coefficients
    inline double __attribute__((always_inline)) compute( double *coeffx, double *coeffy, Field2D *f, int idx, int idy )
    {
        double interp_res( 0. );
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( *f )( idx+iloc, idy+jloc );
            }
        }
        return interp_res;
    };

    //! Interpolator specific to the envelope model
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override ;

    //! Interpolator specific to the envelope model
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override ;

    //! Interpolator specific to the envelope model
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override ;

private:
    inline void __attribute__((always_inline)) coeffs( double xpn, double ypn, int* idx_p, int* idx_d,
                        double *coeffxp, double *coeffyp,
                        double *coeffxd, double *coeffyd, double* delta_p )
    {
        // Indexes of the central nodes
        idx_p[0] = round( xpn );
        idx_d[0] = round( xpn+0.5 );
        idx_p[1] = round( ypn );
        idx_d[1] = round( ypn+0.5 );
        
        // Declaration and calculation of the coefficient for interpolation
        double delta_x, delta_y, delta2, delta3, delta4;
        
        delta_x     = xpn - ( double )idx_d[0] + 0.5;
        delta2     = delta_x*delta_x;
        delta3     = delta2*delta_x;
        delta4     = delta3*delta_x;
        coeffxd[0] = dble_1_ov_384   - dble_1_ov_48  * delta_x  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffxd[1] = dble_19_ov_96   - dble_11_ov_24 * delta_x  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxd[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffxd[3] = dble_19_ov_96   + dble_11_ov_24 * delta_x  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxd[4] = dble_1_ov_384   + dble_1_ov_48  * delta_x  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        delta_p[0] = xpn - ( double )idx_p[0];
        delta2     = delta_p[0]*delta_p[0];
        delta3     = delta2*delta_p[0];
        delta4     = delta3*delta_p[0];
        coeffxp[0] = dble_1_ov_384   - dble_1_ov_48  * delta_p[0]  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffxp[1] = dble_19_ov_96   - dble_11_ov_24 * delta_p[0]  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxp[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffxp[3] = dble_19_ov_96   + dble_11_ov_24 * delta_p[0]  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxp[4] = dble_1_ov_384   + dble_1_ov_48  * delta_p[0]  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        delta_y     = ypn - ( double )idx_d[1] + 0.5;
        delta2     = delta_y*delta_y;
        delta3     = delta2*delta_y;
        delta4     = delta3*delta_y;
        coeffyd[0] = dble_1_ov_384   - dble_1_ov_48  * delta_y  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffyd[1] = dble_19_ov_96   - dble_11_ov_24 * delta_y  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyd[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffyd[3] = dble_19_ov_96   + dble_11_ov_24 * delta_y  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyd[4] = dble_1_ov_384   + dble_1_ov_48  * delta_y  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        delta_p[1] = ypn - ( double )idx_p[1];
        delta2     = delta_p[1]*delta_p[1];
        delta3     = delta2*delta_p[1];
        delta4     = delta3*delta_p[1];
        coeffyp[0] = dble_1_ov_384   - dble_1_ov_48  * delta_p[1]  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffyp[1] = dble_19_ov_96   - dble_11_ov_24 * delta_p[1]  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyp[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffyp[3] = dble_19_ov_96   + dble_11_ov_24 * delta_p[1]  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyp[4] = dble_1_ov_384   + dble_1_ov_48  * delta_p[1]  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        // First index for summation
        idx_p[0]   = idx_p[0] - i_domain_begin;
        idx_d[0]   = idx_d[0] - i_domain_begin;
        idx_p[1]   = idx_p[1] - j_domain_begin;
        idx_d[1]   = idx_d[1] - j_domain_begin;

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
