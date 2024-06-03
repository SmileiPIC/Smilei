#ifndef INTERPOLATOR3D2ORDER_H
#define INTERPOLATOR3D2ORDER_H

#include "Field3D.h"
#include "Interpolator3D.h"
#include "gpu.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 3Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D2Order : public Interpolator3D
{

public:

    //! Creator for Interpolator3D2Order
    Interpolator3D2Order( Params &, Patch * );

    //! Destructor for Interpolator3D2Order
    ~Interpolator3D2Order() override {};

    //! 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );

    //! Interpolation of all fields and currents for a single particles located at istart.
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override ;

    //! Wrapper called by the particle dynamics section
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override ;

    //! Interpolator specific to tracked particles. A selection of particles may be provided
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;

    //! Interpolator on another field than the basic ones
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override ;

    //! Computation of a field from provided coefficients
    inline double __attribute__((always_inline)) compute( double *coeffx, double *coeffy, double *coeffz, const Field3D *const f, int idx, int idy, int idz )
    {
        double interp_res( 0. );
        //unroll ?
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                    interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * *( coeffz+kloc ) * ( *f )( idx+iloc, idy+jloc, idz+kloc );
                }
            }
        }
        return interp_res;
    }

    SMILEI_ACCELERATOR_DECLARE_ROUTINE
    //! Computation of a field from provided coefficients
    static inline double __attribute__((always_inline)) compute( 
        const double *const __restrict__ coeffx, 
        const double *const __restrict__ coeffy, 
        const double *const __restrict__ coeffz, 
        const double *const __restrict__ f, 
        int idx, 
        int idy, 
        int idz, 
        int /*nx*/, 
        int ny, 
        int nz )
    {
        double interp_res( 0. );
        //unroll ?
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                    interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * *( coeffz+kloc ) * f[ (idx+iloc)*ny*nz + (idy+jloc)*nz + (idz+kloc) ];
                }
            }
        }
        return interp_res;
    }
    SMILEI_ACCELERATOR_DECLARE_ROUTINE_END

    //! Interpolator specific to the envelope model
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override;

    //! Interpolator specific to the envelope model
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override;

    //! Interpolator specific to the envelope model
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override;

    //! Interpolator specific to the envelope model
    void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override;

private:

    SMILEI_ACCELERATOR_DECLARE_ROUTINE
    //! Compuation of coefficients for interpolation using particle normalized positions xpn, ypn, zpn
    inline void __attribute__((always_inline)) coeffs( double xpn, double ypn, double zpn, int* idx_p, int* idx_d,
                        double *coeffxp, double *coeffyp, double *coeffzp,
                        double *coeffxd, double *coeffyd, double *coeffzd, double* delta_p ) const
    {
        // Indexes of the central nodes
        idx_p[0] = std::round( xpn );
        idx_p[1] = std::round( ypn );
        idx_p[2] = std::round( zpn );

        // Declaration and calculation of the coefficient for interpolation
        double delta, delta2;

        delta_p[0]   = xpn - ( double )idx_p[0];
        delta2  = delta_p[0]*delta_p[0];
        coeffxp[0] = 0.5 * ( delta2-delta_p[0]+0.25 );
        coeffxp[1] = 0.75 - delta2;
        coeffxp[2] = 0.5 * ( delta2+delta_p[0]+0.25 );

        delta_p[1]   = ypn - ( double )idx_p[1];
        delta2  = delta_p[1]*delta_p[1];
        coeffyp[0] = 0.5 * ( delta2-delta_p[1]+0.25 );
        coeffyp[1] = 0.75 - delta2;
        coeffyp[2] = 0.5 * ( delta2+delta_p[1]+0.25 );

        delta_p[2]   = zpn - ( double )idx_p[2];
        delta2  = delta_p[2]*delta_p[2];
        coeffzp[0] = 0.5 * ( delta2-delta_p[2]+0.25 );
        coeffzp[1] = 0.75 - delta2;
        coeffzp[2] = 0.5 * ( delta2+delta_p[2]+0.25 );

        // First index for summation
        idx_p[0] = idx_p[0] - i_domain_begin;
        idx_p[1] = idx_p[1] - j_domain_begin;
        idx_p[2] = idx_p[2] - k_domain_begin;

        if(idx_d){
            idx_d[0] = std::round( xpn+0.5 );
            idx_d[1] = std::round( ypn+0.5 );
            idx_d[2] = std::round( zpn+0.5 );

            delta   = xpn - ( double )idx_d[0] + 0.5;
            delta2  = delta*delta;
            coeffxd[0] = 0.5 * ( delta2-delta+0.25 );
            coeffxd[1] = 0.75 - delta2;
            coeffxd[2] = 0.5 * ( delta2+delta+0.25 );

            delta   = ypn - ( double )idx_d[1] + 0.5;
            delta2  = delta*delta;
            coeffyd[0] = 0.5 * ( delta2-delta+0.25 );
            coeffyd[1] = 0.75 - delta2;
            coeffyd[2] = 0.5 * ( delta2+delta+0.25 );

            delta   = zpn - ( double )idx_d[2] + 0.5;
            delta2  = delta*delta;
            coeffzd[0] = 0.5 * ( delta2-delta+0.25 );
            coeffzd[1] = 0.75 - delta2;
            coeffzd[2] = 0.5 * ( delta2+delta+0.25 );

            idx_d[0] = idx_d[0] - i_domain_begin;
            idx_d[1] = idx_d[1] - j_domain_begin;
            idx_d[2] = idx_d[2] - k_domain_begin;
        }

    }
    SMILEI_ACCELERATOR_DECLARE_ROUTINE_END

};//END class

#endif
