#ifndef INTERPOLATOR1D4ORDER_H
#define INTERPOLATOR1D4ORDER_H


#include "Interpolator1D.h"
#include "Field1D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 4th order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D4Order final : public Interpolator1D
{

public:
    Interpolator1D4Order( Params &, Patch * );
    ~Interpolator1D4Order() override final {};
    
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
    inline void coeffs( double xpn, int* idx_p, int* idx_d,
                        double *coeffxp, double *coeffxd, double* delta_p )
    {
        double delta, delta2, delta3, delta4 ;
     
                
        // Primal
        idx_p[0]   = round( xpn );            // index of the central point
        delta_p[0] = xpn -( double )idx_p[0]; // normalized distance to the central node
        delta2     = delta_p[0]*delta_p[0];   // square of the normalized distance to the central node
        delta3     = delta2*delta_p[0];       // cube of the normalized distance to the central node
        delta4     = delta3*delta_p[0];       // 4th power of the normalized distance to the central node
        
        // coefficients for the 4th order interpolation on 5 nodes
        coeffxp[0] = dble_1_ov_384   - dble_1_ov_48  * delta_p[0]  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffxp[1] = dble_19_ov_96   - dble_11_ov_24 * delta_p[0]  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxp[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffxp[3] = dble_19_ov_96   + dble_11_ov_24 * delta_p[0]  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxp[4] = dble_1_ov_384   + dble_1_ov_48  * delta_p[0]  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        idx_p[0]  -= i_domain_begin_;

        if(idx_d){   
            // Dual
            idx_d[0]   = round( xpn+0.5 );       // index of the central point
            delta      = xpn -( double )idx_d[0]+0.5; // normalized distance to the central node
            delta2     = delta*delta;            // square of the normalized distance to the central node
            delta3     = delta2*delta;           // cube of the normalized distance to the central node
            delta4     = delta3*delta;           // 4th power of the normalized distance to the central node
            
            // coefficients for the 4th order interpolation on 5 nodes
            coeffxd[0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeffxd[1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeffxd[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeffxd[3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeffxd[4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            
            idx_d[0]  -= i_domain_begin_;
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
