#ifndef INTERPOLATOR1D2ORDER_H
#define INTERPOLATOR1D2ORDER_H


#include "Interpolator1D.h"

#include "Field1D.h"
#include "gpu.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D2Order final : public Interpolator1D
{

public:
    Interpolator1D2Order( Params &, Patch * );
    ~Interpolator1D2Order() override {}; //final
    
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    inline void __attribute__((always_inline)) fieldsForTasks( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc, int *iold, double *delta );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override final;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;

    inline double __attribute__((always_inline)) 
    compute( double *coeff, Field1D *f, int idx )
    {
        double interp_res =  coeff[0] * ( *f )( idx-1 )   + coeff[1] * ( *f )( idx )   + coeff[2] * ( *f )( idx+1 );
        return interp_res;
    }

    SMILEI_ACCELERATOR_DECLARE_ROUTINE
    static inline double __attribute__((always_inline))
    compute( const double *__restrict__ coeff,
             const double *__restrict__ f, 
             int idx )
    {
        double interp_res = coeff[0] * f[idx-1] + coeff[1] * f[idx] + coeff[2] * f[idx+1];
        return interp_res;
    }
    SMILEI_ACCELERATOR_DECLARE_ROUTINE_END

    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;
    void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

private:
    inline void __attribute__((always_inline)) coeffs( double xjn )
    {
        double xjmxi2;

        // Dual
        id_    = std::round( xjn + 0.5 );      // index of the central point
        xjmxi  = xjn - static_cast<double>(id_) + 0.5; // normalized distance to the central node
        xjmxi2 = xjmxi*xjmxi;            // square of the normalized distance to the central node

        // 2nd order interpolation on 3 nodes
        coeffd_[0] = 0.5 * ( xjmxi2-xjmxi + 0.25 );
        coeffd_[1] = ( 0.75 - xjmxi2 );
        coeffd_[2] = 0.5 * ( xjmxi2+xjmxi + 0.25 );

        id_ -= i_domain_begin_;

        // Primal
        ip_      = std::round( xjn );    // index of the central point
        xjmxi  = xjn - static_cast<double>(ip_); // normalized distance to the central node
        xjmxi2 = xjmxi * xjmxi;   // square of the normalized distance to the central node

        // 2nd order interpolation on 3 nodes
        coeffp_[0] = 0.5 * ( xjmxi2 - xjmxi + 0.25 );
        coeffp_[1] = ( 0.75 - xjmxi2 );
        coeffp_[2] = 0.5 * ( xjmxi2 + xjmxi + 0.25 );

        ip_ -= i_domain_begin_;
    }
    
    // 2nd order interpolation on 3 nodes
    SMILEI_ACCELERATOR_DECLARE_ROUTINE
    inline void __attribute__( ( always_inline ) )
    coeffs( double xpn, int* idx_p, int* idx_d,
            double *coeffxp, double *coeffxd, double* delta_p ) const
    {
        double delta, delta2;
        
        // index of the central point
        idx_p[0]   = std::round( xpn );
        idx_d[0]   = std::round( xpn + 0.5 );

        delta      = xpn - static_cast<double>( idx_d[0] ) + 0.5; // normalized distance to the central node
        delta2     = delta * delta;                   // square of the normalized distance to the central node
        
        coeffxd[0] = 0.5 * ( delta2 - delta + 0.25 );
        coeffxd[1] = ( 0.75 - delta2 );
        coeffxd[2] = 0.5 * ( delta2 + delta + 0.25 );

        delta      = xpn - static_cast<double>( idx_p[0] );
        delta2     = delta * delta; // pow( delta_p[0], 2 );   // square of the normalized distance to the central node

        delta_p[0] = delta;   // normalized distance to the central node	
        coeffxp[0] = 0.5 * ( delta2 - delta_p[0] + 0.25 );
        coeffxp[1] = ( 0.75 - delta2 );
        coeffxp[2] = 0.5 * ( delta2 + delta_p[0] + 0.25 );
        
        idx_p[0] = idx_p[0] - i_domain_begin_;
        idx_d[0] = idx_d[0] - i_domain_begin_;
        
    }    
    SMILEI_ACCELERATOR_DECLARE_ROUTINE_END
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Last delta computed
    double xjmxi;
    // Interpolation coefficient on Prim grid
    double coeffp_[3];
    // Interpolation coefficient on Dual grid
    double coeffd_[3];


};//END class

#endif
