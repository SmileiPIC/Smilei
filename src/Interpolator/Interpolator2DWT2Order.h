#ifndef INTERPOLATOR2DWT2ORDER_H
#define INTERPOLATOR2DWT2ORDER_H


#include "Interpolator2D.h"
#include "Field2D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order WT interpolator for 2Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2DWT2Order : public Interpolator2D
{

public:

    //! Creator for Interpolator2DWT2Order
    Interpolator2DWT2Order( Params &, Patch * );

    //! Destructor for Interpolator2DWT2Order
    ~Interpolator2DWT2Order() override {};

    //! 2nd Order WT Interpolation of the fields at a the particle position (3 nodes are used)
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );

    //! Interpolation of all fields and currents for a single particles located at istart.
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override ;

    //! Wrapper called by the particle dynamics section
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override ;

    //! Interpolator specific to tracked particles. A selection of particles may be provided
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override;

    //! Interpolator on another field than the basic ones
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override;

    //! Computation of a field from provided coefficients
    inline double __attribute__((always_inline)) compute( double *coeffx, double *coeffy, Field2D *f, int idx, int idy )
    {
        double interp_res( 0. );
        //unroll ?
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
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

    //! Interpolator specific to the envelope model
    void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override;

private:

    //! Compuation of coefficients for interpolation using particle normalized positions xpn and ypn
    inline void __attribute__((always_inline)) coeffs( double xpn, double ypn )
    {
        // Indexes of the central nodes
        ip_ = round( xpn );
        id_ = round( xpn+0.5 );
        jp_ = round( ypn );
        jd_ = round( ypn+0.5 );

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
        
        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        id_ = id_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
        jd_ = jd_ - j_domain_begin;
    }

    // Coefficients for WT
    double D_ov_8dt [2];
    double dt_ov_D [2];
    
    // Last prim index computed
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Last delta computed
    double deltax, deltay;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3];
    // WT Interpolation coefficient on Prim grid
    double coeffxpt_[3], coeffypt_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3];


};//END class

#endif
