#ifndef INTERPOLATOR3D2ORDER_H
#define INTERPOLATOR3D2ORDER_H


#include "Interpolator3D.h"
#include "Field3D.h"


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
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override ;

    //! Interpolator specific to tracked particles. A selection of particles may be provided
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;

    //! Interpolator on another field than the basic ones
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override ;

    //! Computation of a field from provided coefficients
    inline double __attribute__((always_inline)) compute( double *coeffx, double *coeffy, double *coeffz, Field3D *f, int idx, int idy, int idz )
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
    };

    //! Computation of a field from provided coefficients
    inline double __attribute__((always_inline)) compute( double *coeffx, double *coeffy, double *coeffz, double *f, int idx, int idy, int idz, int nx, int ny, int nz )
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
    };

    //! Interpolator specific to the envelope model
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override ;

    //! Interpolator specific to the envelope model
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override ;

    //! Interpolator specific to the envelope model
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override ;

    //! Interpolator specific to the envelope model
    void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

private:

    //! Compuation of coefficients for interpolation using particle normalized positions xpn, ypn, zpn
    inline void __attribute__((always_inline)) coeffs( double xpn, double ypn, double zpn )
    {
        // Indexes of the central nodes
        ip_ = round( xpn );
        id_ = round( xpn+0.5 );
        jp_ = round( ypn );
        jd_ = round( ypn+0.5 );
        kp_ = round( zpn );
        kd_ = round( zpn+0.5 );

        // Declaration and calculation of the coefficient for interpolation
        double delta2;

        deltax   = xpn - ( double )id_ + 0.5;
        delta2  = deltax*deltax;
        coeffxd_[0] = 0.5 * ( delta2-deltax+0.25 );
        coeffxd_[1] = 0.75 - delta2;
        coeffxd_[2] = 0.5 * ( delta2+deltax+0.25 );

        deltax   = xpn - ( double )ip_;
        delta2  = deltax*deltax;
        coeffxp_[0] = 0.5 * ( delta2-deltax+0.25 );
        coeffxp_[1] = 0.75 - delta2;
        coeffxp_[2] = 0.5 * ( delta2+deltax+0.25 );

        deltay   = ypn - ( double )jd_ + 0.5;
        delta2  = deltay*deltay;
        coeffyd_[0] = 0.5 * ( delta2-deltay+0.25 );
        coeffyd_[1] = 0.75 - delta2;
        coeffyd_[2] = 0.5 * ( delta2+deltay+0.25 );

        deltay   = ypn - ( double )jp_;
        delta2  = deltay*deltay;
        coeffyp_[0] = 0.5 * ( delta2-deltay+0.25 );
        coeffyp_[1] = 0.75 - delta2;
        coeffyp_[2] = 0.5 * ( delta2+deltay+0.25 );

        deltaz   = zpn - ( double )kd_ + 0.5;
        delta2  = deltaz*deltaz;
        coeffzd_[0] = 0.5 * ( delta2-deltaz+0.25 );
        coeffzd_[1] = 0.75 - delta2;
        coeffzd_[2] = 0.5 * ( delta2+deltaz+0.25 );

        deltaz   = zpn - ( double )kp_;
        delta2  = deltaz*deltaz;
        coeffzp_[0] = 0.5 * ( delta2-deltaz+0.25 );
        coeffzp_[1] = 0.75 - delta2;
        coeffzp_[2] = 0.5 * ( delta2+deltaz+0.25 );

        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        id_ = id_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
        jd_ = jd_ - j_domain_begin;
        kp_ = kp_ - k_domain_begin;
        kd_ = kd_ - k_domain_begin;
    }

    //! Compuation of coefficients for interpolation using particle normalized positions xpn, ypn, zpn
    inline void __attribute__((always_inline)) coeffs( double xpn, double ypn, double zpn, int* idx_p, int* idx_d,
                        double *coeffxp, double *coeffyp, double *coeffzp,
                        double *coeffxd, double *coeffyd, double *coeffzd, double* delta_p )
    {
        // Indexes of the central nodes
        idx_p[0] = round( xpn );
        idx_d[0] = round( xpn+0.5 );
        idx_p[1] = round( ypn );
        idx_d[1] = round( ypn+0.5 );
        idx_p[2] = round( zpn );
        idx_d[2] = round( zpn+0.5 );

        // Declaration and calculation of the coefficient for interpolation
        double delta, delta2;

        delta   = xpn - ( double )idx_d[0] + 0.5;
        delta2  = delta*delta;
        coeffxd[0] = 0.5 * ( delta2-delta+0.25 );
        coeffxd[1] = 0.75 - delta2;
        coeffxd[2] = 0.5 * ( delta2+delta+0.25 );

        delta_p[0]   = xpn - ( double )idx_p[0];
        delta2  = delta_p[0]*delta_p[0];
        coeffxp[0] = 0.5 * ( delta2-delta_p[0]+0.25 );
        coeffxp[1] = 0.75 - delta2;
        coeffxp[2] = 0.5 * ( delta2+delta_p[0]+0.25 );

        delta   = ypn - ( double )idx_d[1] + 0.5;
        delta2  = delta*delta;
        coeffyd[0] = 0.5 * ( delta2-delta+0.25 );
        coeffyd[1] = 0.75 - delta2;
        coeffyd[2] = 0.5 * ( delta2+delta+0.25 );

        delta_p[1]   = ypn - ( double )idx_p[1];
        delta2  = delta_p[1]*delta_p[1];
        coeffyp[0] = 0.5 * ( delta2-delta_p[1]+0.25 );
        coeffyp[1] = 0.75 - delta2;
        coeffyp[2] = 0.5 * ( delta2+delta_p[1]+0.25 );

        delta   = zpn - ( double )idx_d[2] + 0.5;
        delta2  = delta*delta;
        coeffzd[0] = 0.5 * ( delta2-delta+0.25 );
        coeffzd[1] = 0.75 - delta2;
        coeffzd[2] = 0.5 * ( delta2+delta+0.25 );

        delta_p[2]   = zpn - ( double )idx_p[2];
        delta2  = delta_p[2]*delta_p[2];
        coeffzp[0] = 0.5 * ( delta2-delta_p[2]+0.25 );
        coeffzp[1] = 0.75 - delta2;
        coeffzp[2] = 0.5 * ( delta2+delta_p[2]+0.25 );

        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        idx_p[0] = idx_p[0] - i_domain_begin;
        idx_d[0] = idx_d[0] - i_domain_begin;
        idx_p[1] = idx_p[1] - j_domain_begin;
        idx_d[1] = idx_d[1] - j_domain_begin;
        idx_p[2] = idx_p[2] - k_domain_begin;
        idx_d[2] = idx_d[2] - k_domain_begin;
    }

    // Last prim index computed
    int ip_, jp_, kp_;
    // Last dual index computed
    int id_, jd_, kd_;
    // Last delta computed
    double deltax, deltay, deltaz;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3], coeffzp_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3], coeffzd_[3];


};//END class

#endif
