#ifndef INTERPOLATOR3D4ORDER_H
#define INTERPOLATOR3D4ORDER_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D4Order final : public Interpolator3D
{

public:
    Interpolator3D4Order( Params &, Patch * );
    ~Interpolator3D4Order() override final {};
    
    inline void fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    inline void fieldsForTasks( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc, int *iold, double *delta );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final ;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;
    
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;
    
    inline double compute( double *coeffx, double *coeffy, double *coeffz, Field3D *f, int idx, int idy, int idz )
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
    inline void coeffs( double xpn, double ypn, double zpn )
    {
        // Indexes of the central nodes
        ip_ = round( xpn );
        id_ = round( xpn+0.5 );
        jp_ = round( ypn );
        jd_ = round( ypn+0.5 );
        kp_ = round( zpn );
        kd_ = round( zpn+0.5 );
        
        // Declaration and calculation of the coefficient for interpolation
        double delta2, delta3, delta4;
        
        deltax   = xpn - ( double )id_ + 0.5;
        delta2  = deltax*deltax;
        delta3  = delta2*deltax;
        delta4  = delta3*deltax;
        coeffxd_[0] = dble_1_ov_384   - dble_1_ov_48  * deltax  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffxd_[1] = dble_19_ov_96   - dble_11_ov_24 * deltax  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxd_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffxd_[3] = dble_19_ov_96   + dble_11_ov_24 * deltax  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxd_[4] = dble_1_ov_384   + dble_1_ov_48  * deltax  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        deltax   = xpn - ( double )ip_;
        delta2  = deltax*deltax;
        delta3  = delta2*deltax;
        delta4  = delta3*deltax;
        coeffxp_[0] = dble_1_ov_384   - dble_1_ov_48  * deltax  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffxp_[1] = dble_19_ov_96   - dble_11_ov_24 * deltax  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxp_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffxp_[3] = dble_19_ov_96   + dble_11_ov_24 * deltax  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffxp_[4] = dble_1_ov_384   + dble_1_ov_48  * deltax  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        deltay   = ypn - ( double )jd_ + 0.5;
        delta2  = deltay*deltay;
        delta3  = delta2*deltay;
        delta4  = delta3*deltay;
        coeffyd_[0] = dble_1_ov_384   - dble_1_ov_48  * deltay  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffyd_[1] = dble_19_ov_96   - dble_11_ov_24 * deltay  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyd_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffyd_[3] = dble_19_ov_96   + dble_11_ov_24 * deltay  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyd_[4] = dble_1_ov_384   + dble_1_ov_48  * deltay  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        deltay   = ypn - ( double )jp_;
        delta2  = deltay*deltay;
        delta3  = delta2*deltay;
        delta4  = delta3*deltay;
        coeffyp_[0] = dble_1_ov_384   - dble_1_ov_48  * deltay  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffyp_[1] = dble_19_ov_96   - dble_11_ov_24 * deltay  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyp_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffyp_[3] = dble_19_ov_96   + dble_11_ov_24 * deltay  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffyp_[4] = dble_1_ov_384   + dble_1_ov_48  * deltay  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        deltaz   = zpn - ( double )kd_ + 0.5;
        delta2  = deltaz*deltaz;
        delta3  = delta2*deltaz;
        delta4  = delta3*deltaz;
        coeffzd_[0] = dble_1_ov_384   - dble_1_ov_48  * deltaz  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffzd_[1] = dble_19_ov_96   - dble_11_ov_24 * deltaz  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzd_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffzd_[3] = dble_19_ov_96   + dble_11_ov_24 * deltaz  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzd_[4] = dble_1_ov_384   + dble_1_ov_48  * deltaz  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        deltaz   = zpn - ( double )kp_;
        delta2  = deltaz*deltaz;
        delta3  = delta2*deltaz;
        delta4  = delta3*deltaz;
        coeffzp_[0] = dble_1_ov_384   - dble_1_ov_48  * deltaz  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffzp_[1] = dble_19_ov_96   - dble_11_ov_24 * deltaz  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzp_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffzp_[3] = dble_19_ov_96   + dble_11_ov_24 * deltaz  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzp_[4] = dble_1_ov_384   + dble_1_ov_48  * deltaz  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        id_ = id_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
        jd_ = jd_ - j_domain_begin;
        kp_ = kp_ - k_domain_begin;
        kd_ = kd_ - k_domain_begin;
    };

    inline void coeffs( double xpn, double ypn, double zpn, int* idx_p, int* idx_d,
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
        double delta_x, delta_y, delta_z, delta2, delta3, delta4;
        
        delta_x    = xpn - ( double )idx_d[0] + 0.5;
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
        
        delta_y    = ypn - ( double )idx_d[1] + 0.5;
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
        
        delta_z    = zpn - ( double )idx_d[2] + 0.5;
        delta2     = delta_z*delta_z;
        delta3     = delta2*delta_z;
        delta4     = delta3*delta_z;
        coeffzd[0] = dble_1_ov_384   - dble_1_ov_48  * delta_z  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffzd[1] = dble_19_ov_96   - dble_11_ov_24 * delta_z  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzd[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffzd[3] = dble_19_ov_96   + dble_11_ov_24 * delta_z  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzd[4] = dble_1_ov_384   + dble_1_ov_48  * delta_z  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        delta_p[2] = zpn - ( double )idx_p[2];
        delta2     = delta_p[2]*delta_p[2];
        delta3     = delta2*delta_p[2];
        delta4     = delta3*delta_p[2];
        coeffzp[0] = dble_1_ov_384   - dble_1_ov_48  * delta_p[2]  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        coeffzp[1] = dble_19_ov_96   - dble_11_ov_24 * delta_p[2]  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzp[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
        coeffzp[3] = dble_19_ov_96   + dble_11_ov_24 * delta_p[2]  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
        coeffzp[4] = dble_1_ov_384   + dble_1_ov_48  * delta_p[2]  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
        
        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        idx_p[0]   = idx_p[0] - i_domain_begin;
        idx_d[0]   = idx_d[0] - i_domain_begin;
        idx_p[1]   = idx_p[1] - j_domain_begin;
        idx_d[1]   = idx_d[1] - j_domain_begin;
        idx_p[2]   = idx_p[2] - k_domain_begin;
        idx_d[2]   = idx_d[2] - k_domain_begin;

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
    
    // Last prim index computed
    int ip_, jp_, kp_;
    // Last dual index computed
    int id_, jd_, kd_;
    // Last delta computed
    double deltax, deltay, deltaz;
    // Interpolation coefficient on Prim grid
    double coeffxp_[5], coeffyp_[5], coeffzp_[5];
    // Interpolation coefficient on Dual grid
    double coeffxd_[5], coeffyd_[5], coeffzd_[5];
    
    
};//END class

#endif
