#ifndef INTERPOLATOR3DWT4ORDER_H
#define INTERPOLATOR3DWT4ORDER_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 4th order WT interpolator for 3Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3DWT4Order : public Interpolator3D
{

public:

    //! Creator for Interpolator3D4Order
    Interpolator3DWT4Order( Params &, Patch * );

    //! Destructor for Interpolator3D4Order
    ~Interpolator3DWT4Order() override {};

    //! 4th Order WT Interpolation of the fields at a the particle position
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    //inline void fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );

    //! Interpolation of all fields and currents for a single particles located at istart.
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override ;

    //! Wrapper called by the particle dynamics section
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override ;

    //! Interpolator specific to tracked particles. A selection of particles may be provided
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
    //inline double compute( double *coeffx, double *coeffy, double *coeffz, Field3D *f, int idx, int idy, int idz )
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

    //! Compuation of coefficients for interpolation using particle normalized positions xpn, ypn, zpn
    inline void __attribute__((always_inline)) coeffs( double xpn, double ypn, double zpn )
    //inline void coeffs( double xpn, double ypn, double zpn )
    {
        // Indexes of the central nodes
        ip_ = round( xpn );
        id_ = round( xpn+0.5 );
        jp_ = round( ypn );
        jd_ = round( ypn+0.5 );
        kp_ = round( zpn );
        kd_ = round( zpn+0.5 );

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
        
        deltaz   = zpn - ( double )kd_ + 0.5;
        var1 = dble_11_ov_24 * deltaz;
        var2 = deltaz * deltaz;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - deltaz;
        var4 = var4 * var4;
        var5 = 0.5 + deltaz;
        var5 = var5 * var5;
        coeffzd_[0] = dble_1_ov_24 * var4 * var4;
        coeffzd_[1] = dble_19_ov_96 - var1 + var3 * ( 1.5+deltaz-var2 );
        coeffzd_[2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeffzd_[3] = dble_19_ov_96 + var1 + var3 * ( 1.5-deltaz -var2 );
        coeffzd_[4] = dble_1_ov_24 * var5 * var5;
        
        deltaz   = zpn - ( double )kp_;
        var1 = dble_11_ov_24 * deltaz;
        var2 = deltaz * deltaz;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - deltaz;
        var4 = var4 * var4;
        var5 = 0.5 + deltaz;
        var5 = var5 * var5;
        coeffzp_[0] = dble_1_ov_24 * var4 * var4;
        coeffzp_[1] = dble_19_ov_96 - var1 + var3 * ( 1.5+deltaz-var2 );
        coeffzp_[2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeffzp_[3] = dble_19_ov_96 + var1 + var3 * ( 1.5-deltaz-var2 );
        coeffzp_[4] = dble_1_ov_24 * var5 * var5;

        var1 = dt_ov_D[2] - deltaz;
        var1 = var1 * var1;
        var1 = D_ov_96dt[2] * var1 * var1;
        var2 = dt_ov_D[2] + deltaz;
        var2 = var2 * var2;
        var2 = D_ov_96dt[2] * var2 * var2;
        var3 = copysign( var1, dt_ov_D[2]-deltaz ) + copysign( var2, dt_ov_D[2]+deltaz);
        var4 =  dble_1_ov_24 * ((( 3.0+deltaz ) * deltaz  - ( 3.0-dt2_ov_D2[2] )) * deltaz + (1.0+dt2_ov_D2[2] ));
        var5 =  dble_1_ov_24 * ((( 3.0-deltaz ) * deltaz  + ( 3.0-dt2_ov_D2[2] )) * deltaz + (1.0+dt2_ov_D2[2] ));
        coeffzpt_[0] = var3 + var1 - var2;
        coeffzpt_[1] = 4.0 * ( var4-var3 );
        coeffzpt_[2] = 1.0 + 6.0 * var3 - 4.0 * ( var4+var5 );
        coeffzpt_[3] = 4.0 * ( var5-var3 );
        coeffzpt_[4] = var3 + var2 - var1;
        
        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        id_ = id_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
        jd_ = jd_ - j_domain_begin;
        kp_ = kp_ - k_domain_begin;
        kd_ = kd_ - k_domain_begin;
    };

    double dble_1_ov_6 ;
    double dble_1_ov_24 ;
    double dble_11_ov_24 ;
    double dble_19_ov_96 ;
    double dble_115_ov_192 ;
    double dt_ov_D[3] ;
    double dt2_ov_D2[3] ;
    double D_ov_96dt[3] ;
    
    // Last prim index computed
    int ip_, jp_, kp_;
    // Last dual index computed
    int id_, jd_, kd_;
    // Last delta computed
    double deltax, deltay, deltaz;
    // Interpolation coefficient on Prim grid
    double coeffxp_[5], coeffyp_[5], coeffzp_[5];
    // WT interpolation coefficient on Prim grid
    double coeffxpt_[5], coeffypt_[5], coeffzpt_[5];
    // Interpolation coefficient on Dual grid
    double coeffxd_[5], coeffyd_[5], coeffzd_[5];


};//END class

#endif
