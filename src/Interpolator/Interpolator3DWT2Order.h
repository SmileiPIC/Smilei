#ifndef INTERPOLATOR3DWT2ORDER_H
#define INTERPOLATOR3DWT2ORDER_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order WT interpolator for 3Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3DWT2Order : public Interpolator3D
{

public:

    //! Creator for Interpolator3DWT2Order
    Interpolator3DWT2Order( Params &, Patch * );

    //! Destructor for Interpolator3DWT2Order
    ~Interpolator3DWT2Order() override {};

    //! 2nd Order WT Interpolation of the fields at a the particle position (3 nodes are used)
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
        double var1;

        deltax   = xpn - ( double )id_ + 0.5;
        var1  = deltax*deltax;
        coeffxd_[0] = 0.5 * ( var1-deltax+0.25 );
        coeffxd_[1] = 0.75 - var1;
        coeffxd_[2] = 0.5 * ( var1+deltax+0.25 );

        deltax   = xpn - ( double )ip_;
        var1  = deltax*deltax;
        coeffxp_[0] = 0.5 * ( var1-deltax+0.25 );
        coeffxp_[1] = 0.75 - var1;
        coeffxp_[2] = 0.5 * ( var1+deltax+0.25 );

        var1 = D_ov_8dt[0] * ( fabs( dt_ov_D[0]-deltax ) * ( dt_ov_D[0]-deltax ) + fabs( dt_ov_D[0]+deltax ) * ( dt_ov_D[0]+deltax ) );
        coeffxpt_[0] = var1 - 0.5 * deltax;
        coeffxpt_[1] = 1.0 - 2.0 * var1;
        coeffxpt_[2] = var1 + 0.5 * deltax;
        
        deltay   = ypn - ( double )jd_ + 0.5;
        var1  = deltay*deltay;
        coeffyd_[0] = 0.5 * ( var1-deltay+0.25 );
        coeffyd_[1] = 0.75 - var1;
        coeffyd_[2] = 0.5 * ( var1+deltay+0.25 );

        deltay   = ypn - ( double )jp_;
        var1  = deltay*deltay;
        coeffyp_[0] = 0.5 * ( var1-deltay+0.25 );
        coeffyp_[1] = 0.75 - var1;
        coeffyp_[2] = 0.5 * ( var1+deltay+0.25 );

        var1 = D_ov_8dt[1] * ( fabs( dt_ov_D[1]-deltay ) * ( dt_ov_D[1]-deltay ) + fabs( dt_ov_D[1]+deltay ) * ( dt_ov_D[1]+deltay ) );
        coeffypt_[0] = var1 - 0.5 * deltay;
        coeffypt_[1] = 1.0 - 2.0 * var1;
        coeffypt_[2] = var1 + 0.5 * deltay;
        
        deltaz   = zpn - ( double )kd_ + 0.5;
        var1  = deltaz*deltaz;
        coeffzd_[0] = 0.5 * ( var1-deltaz+0.25 );
        coeffzd_[1] = 0.75 - var1;
        coeffzd_[2] = 0.5 * ( var1+deltaz+0.25 );

        deltaz   = zpn - ( double )kp_;
        var1  = deltaz*deltaz;
        coeffzp_[0] = 0.5 * ( var1-deltaz+0.25 );
        coeffzp_[1] = 0.75 - var1;
        coeffzp_[2] = 0.5 * ( var1+deltaz+0.25 );

        var1 = D_ov_8dt[2] * ( fabs( dt_ov_D[2]-deltaz ) * ( dt_ov_D[2]-deltaz ) + fabs( dt_ov_D[2]+deltaz ) * ( dt_ov_D[2]+deltaz ) );
        coeffzpt_[0] = var1 - 0.5 * deltaz;
        coeffzpt_[1] = 1.0 - 2.0 * var1;
        coeffzpt_[2] = var1 + 0.5 * deltaz;
        
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
                        double *coeffxpt, double *coeffypt, double *coeffzpt,
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
        double delta, var1;

        delta   = xpn - ( double )idx_d[0] + 0.5;
        var1  = delta*delta;
        coeffxd[0] = 0.5 * ( var1-delta+0.25 );
        coeffxd[1] = 0.75 - var1;
        coeffxd[2] = 0.5 * ( var1+delta+0.25 );

        delta_p[0]   = xpn - ( double )idx_p[0];
        var1 = D_ov_8dt[0] * ( fabs( dt_ov_D[0]-delta_p[0] ) * ( dt_ov_D[0]-delta_p[0] ) + fabs( dt_ov_D[0]+delta_p[0] ) * ( dt_ov_D[0]+delta_p[0] ) );
        coeffxpt[0] = var1 - 0.5 * delta_p[0];
        coeffxpt[1] = 1.0 - 2.0 * var1;
        coeffxpt[2] = var1 + 0.5 * delta_p[0];

        delta   = ypn - ( double )idx_d[1] + 0.5;
        var1  = delta*delta;
        coeffyd[0] = 0.5 * ( var1-delta+0.25 );
        coeffyd[1] = 0.75 - var1;
        coeffyd[2] = 0.5 * ( var1+delta+0.25 );

        delta_p[1]   = ypn - ( double )idx_p[1];
        var1 = D_ov_8dt[1] * ( fabs( dt_ov_D[1]-delta_p[1] ) * ( dt_ov_D[1]-delta_p[1] ) + fabs( dt_ov_D[1]+delta_p[1] ) * ( dt_ov_D[1]+delta_p[1] ) );
        coeffypt[0] = var1 - 0.5 * delta_p[1];
        coeffypt[1] = 1.0 - 2.0 * var1;
        coeffypt[2] = var1 + 0.5 * delta_p[1];

        delta   = zpn - ( double )idx_d[2] + 0.5;
        var1  = delta*delta;
        coeffzd[0] = 0.5 * ( var1-delta+0.25 );
        coeffzd[1] = 0.75 - var1;
        coeffzd[2] = 0.5 * ( var1+delta+0.25 );

        delta_p[2]   = zpn - ( double )idx_p[2];
        var1 = D_ov_8dt[2] * ( fabs( dt_ov_D[2]-delta_p[2] ) * ( dt_ov_D[2]-delta_p[2] ) + fabs( dt_ov_D[2]+delta_p[2] ) * ( dt_ov_D[2]+delta_p[2] ) );
        coeffzpt[0] = var1 - 0.5 * delta_p[2];
        coeffzpt[1] = 1.0 - 2.0 * var1;
        coeffzpt[2] = var1 + 0.5 * delta_p[2];

        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        idx_p[0] = idx_p[0] - i_domain_begin;
        idx_d[0] = idx_d[0] - i_domain_begin;
        idx_p[1] = idx_p[1] - j_domain_begin;
        idx_d[1] = idx_d[1] - j_domain_begin;
        idx_p[2] = idx_p[2] - k_domain_begin;
        idx_d[2] = idx_d[2] - k_domain_begin;
    }

    // Coefficients for WT
    double dt_ov_D[3];
    double D_ov_8dt[3];
    
    // Last prim index computed
    int ip_, jp_, kp_;
    // Last dual index computed
    int id_, jd_, kd_;
    // Last delta computed
    double deltax, deltay, deltaz;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3], coeffzp_[3];
    // WT interpolation coefficient on Prim grid
    double coeffxpt_[3], coeffypt_[3], coeffzpt_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3], coeffzd_[3];


};//END class

#endif
