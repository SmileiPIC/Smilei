#ifndef INTERPOLATOR2DWT4ORDER_H
#define INTERPOLATOR2DWT4ORDER_H

#include <cmath>

#include "Interpolator2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 4th order WT interpolator for 2Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2DWT4Order : public Interpolator2D
{

public:

    //! Creator for Interpolator2DWT2Order
    Interpolator2DWT4Order( Params &, Patch * );

    //! Destructor for Interpolator2DWT2Order
    ~Interpolator2DWT4Order() override {};

    //! 4th Order WT Interpolation of the fields at a the particle position
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );

    //! Interpolation of all fields and currents for a single particles located at istart.
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;

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

    //! Compuation of coefficients for interpolation using particle normalized positions xpn and ypn
    inline void __attribute__((always_inline)) coeffs( double xpn, double ypn )
    {
        // Indexes of the central nodes
        ip_ = round( xpn );
        id_ = round( xpn+0.5 );
        jp_ = round( ypn );
        jd_ = round( ypn+0.5 );

        // Declaration and calculation of the coefficient for interpolation
        double var1, var2, var3, var4, var5;

        deltax   = xpn - ( double )id_ + 0.5;
        var1 = dble_11_ov_24 * deltax;
        var2 = deltax * deltax;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - deltax;
        var4 = var4 * var4;
        var5 = 0.5 + deltax;
        var5 = var5 * var5;
        coeffxd_[0] = dble_1_ov_24 * var4 * var4;
        coeffxd_[1] = dble_19_ov_96 - var1 + var3 * ( 1.5+deltax-var2 );
        coeffxd_[2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeffxd_[3] = dble_19_ov_96 + var1 + var3 * ( 1.5-deltax -var2 );
        coeffxd_[4] = dble_1_ov_24 * var5 * var5;
        
        deltax   = xpn - ( double )ip_;
        var1 = dble_11_ov_24 * deltax;
        var2 = deltax * deltax;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - deltax;
        var4 = var4 * var4;
        var5 = 0.5 + deltax;
        var5 = var5 * var5;
        coeffxp_[0] = dble_1_ov_24 * var4 * var4;
        coeffxp_[1] = dble_19_ov_96 - var1 + var3 * ( 1.5+deltax-var2 );
        coeffxp_[2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeffxp_[3] = dble_19_ov_96 + var1 + var3 * ( 1.5-deltax-var2 );
        coeffxp_[4] = dble_1_ov_24 * var5 * var5;

        // std::cerr << " xpn: " << xpn
        //           << " delta: " << deltax
        //           << " coeff: "<< coeffxp_[0]
        //           << std::endl;

        var1 = dt_ov_D[0] - deltax;
        var1 = var1 * var1;
        var1 = D_ov_96dt[0] * var1 * var1;
        var2 = dt_ov_D[0] + deltax;
        var2 = var2 * var2;
        var2 = D_ov_96dt[0] * var2 * var2;
        var3 = copysign( var1, dt_ov_D[0]-deltax ) + copysign( var2, dt_ov_D[0]+deltax);
        var4 =  dble_1_ov_24 * ((( 3.0+deltax ) * deltax  - ( 3.0-dt2_ov_D2[0] )) * deltax + (1.0+dt2_ov_D2[0] ));
        var5 =  dble_1_ov_24 * ((( 3.0-deltax ) * deltax  + ( 3.0-dt2_ov_D2[0] )) * deltax + (1.0+dt2_ov_D2[0] ));
        coeffxpt_[0] = var3 + var1 - var2;
        coeffxpt_[1] = 4.0 * ( var4-var3 );
        coeffxpt_[2] = 1.0 + 6.0 * var3 - 4.0 * ( var4+var5 );
        coeffxpt_[3] = 4.0 * ( var5-var3 );
        coeffxpt_[4] = var3 + var2 - var1;
        
        deltay   = ypn - ( double )jd_ + 0.5;
        var1 = dble_11_ov_24 * deltay;
        var2 = deltay * deltay;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - deltay;
        var4 = var4 * var4;
        var5 = 0.5 + deltay;
        var5 = var5 * var5;
        coeffyd_[0] = dble_1_ov_24 * var4 * var4;
        coeffyd_[1] = dble_19_ov_96 - var1 + var3 * ( 1.5+deltay-var2 );
        coeffyd_[2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeffyd_[3] = dble_19_ov_96 + var1 + var3 * ( 1.5-deltay -var2 );
        coeffyd_[4] = dble_1_ov_24 * var5 * var5;
        
        deltay   = ypn - ( double )jp_;
        var1 = dble_11_ov_24 * deltay;
        var2 = deltay * deltay;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - deltay;
        var4 = var4 * var4;
        var5 = 0.5 + deltay;
        var5 = var5 * var5;
        coeffyp_[0] = dble_1_ov_24 * var4 * var4;
        coeffyp_[1] = dble_19_ov_96 - var1 + var3 * ( 1.5+deltay-var2 );
        coeffyp_[2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeffyp_[3] = dble_19_ov_96 + var1 + var3 * ( 1.5-deltay-var2 );
        coeffyp_[4] = dble_1_ov_24 * var5 * var5;

        var1 = dt_ov_D[1] - deltay;
        var1 = var1 * var1;
        var1 = D_ov_96dt[1] * var1 * var1;
        var2 = dt_ov_D[1] + deltay;
        var2 = var2 * var2;
        var2 = D_ov_96dt[1] * var2 * var2;
        var3 = copysign( var1, dt_ov_D[1]-deltay ) + copysign( var2, dt_ov_D[1]+deltay);
        var4 =  dble_1_ov_24 * ((( 3.0+deltay ) * deltay  - ( 3.0-dt2_ov_D2[1] )) * deltay + (1.0+dt2_ov_D2[1] ));
        var5 =  dble_1_ov_24 * ((( 3.0-deltay ) * deltay  + ( 3.0-dt2_ov_D2[1] )) * deltay + (1.0+dt2_ov_D2[1] ));
        coeffypt_[0] = var3 + var1 - var2;
        coeffypt_[1] = 4.0 * ( var4-var3 );
        coeffypt_[2] = 1.0 + 6.0 * var3 - 4.0 * ( var4+var5 );
        coeffypt_[3] = 4.0 * ( var5-var3 );
        coeffypt_[4] = var3 + var2 - var1;
        
        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        id_ = id_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
        jd_ = jd_ - j_domain_begin;

        // std::cerr << " coefxp: " << coeffxp_[0]
        //           << " " << coeffxp_[1]
        //           << " " << coeffxp_[2]
        //           << " " << coeffxp_[3]
        //           << " " << coeffxp_[4] << std::endl;
        //
        // std::cerr << " coefyp: " << coeffyp_[0]
        //           << " " << coeffyp_[1]
        //           << " " << coeffyp_[2]
        //           << " " << coeffyp_[3]
        //           << " " << coeffyp_[4] << std::endl;
        //
        // std::cerr << " coefxd: " << coeffxd_[0]
        //           << " " << coeffxd_[1]
        //           << " " << coeffxd_[2]
        //           << " " << coeffxd_[3]
        //           << " " << coeffxd_[4] << std::endl;
        //
        // std::cerr << " coefyd: " << coeffyd_[0]
        //           << " " << coeffyd_[1]
        //           << " " << coeffyd_[2]
        //           << " " << coeffyd_[3]
        //           << " " << coeffyd_[4] << std::endl;
        //
        // std::cerr << " xpn: " << xpn
        //           << " ypn: " << ypn
        //           << " ip: " << ip_
        //           << " id: " << id_
        //           << " jp: " << jp_
        //           << " jd: " << jd_ << std::endl;

    }


    double dble_1_ov_6 ;
    double dble_1_ov_24 ;
    double dble_11_ov_24 ;
    double dble_19_ov_96 ;
    double dble_115_ov_192 ;
    double dt_ov_D [2] ;
    double dt2_ov_D2 [2] ;
    double D_ov_96dt [2] ;

    // Last prim index computed
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Last delta computed
    double deltax, deltay;
    // Interpolation coefficient on Prim grid
    double coeffxp_[5], coeffyp_[5];
    // WT interpolation coefficient on Prim grid
    double coeffxpt_[5], coeffypt_[5];
    // Interpolation coefficient on Dual grid
    double coeffxd_[5], coeffyd_[5];

};//END class

#endif
