#ifndef INTERPOLATORAM2ORDER_H
#define INTERPOLATORAM2ORDER_H


#include "InterpolatorAM.h"
#include "cField2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for AM spectral simulations
//  --------------------------------------------------------------------------------------------------------------------
class InterpolatorAM2Order final : public InterpolatorAM
{

public:
    InterpolatorAM2Order( Params &, Patch * );
    ~InterpolatorAM2Order() override final {};
    
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    inline void __attribute__((always_inline)) fieldsForTasks( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc, int *iold, double *delta, std::complex<double> *eitheta_old  );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override final ;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;

    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void fieldsAndEnvelopeForTasks( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void timeCenteredEnvelopeForTasks( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;
    void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void envelopeFieldForIonizationTasks( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    inline std::complex<double> __attribute__((always_inline)) compute( double *coeffx, double *coeffy, cField2D *f, int idx, int idy )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                //std::cout << "Er standard idx = " << idx << " idy = " << idy << " iloc = " << iloc << " jloc = " << jloc << " " <<  *( coeffx+iloc ) << " " <<  *( coeffy+jloc ) << " " <<  ( *f )( idx+iloc, idy+jloc )  << std::endl;
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;

                //std::cout<<"f "<<std::fixed << std::setprecision(3)<<(*f)(idx+iloc,idy+jloc)<<std::endl;
            }
        }
        //std::cout<<"interp res "<< interp_res <<std::endl;
        return interp_res;
    };

    inline std::complex<double> __attribute__((always_inline)) compute_0_T( double *coeffx, double *coeffy, cField2D *f, int idx, int idy )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                if( jloc+idy+j_domain_begin_==0 ) {
                    interp_res -= *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
                } else {
                    interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
                }
            }
        }
        return interp_res;
    };
    inline std::complex<double> __attribute__((always_inline)) compute_0_L( double *coeffx, double *coeffy, cField2D *f, int idx, int idy )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
            }
        }
        return interp_res;
    };
    inline std::complex<double> __attribute__((always_inline)) compute_1_T( double *coeffx, double *coeffy, cField2D *f, int idx, int idy, std::complex<double> *exptheta )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) )*( *exptheta ) ;
            }
        }
        return interp_res;
    };
    inline std::complex<double> __attribute__((always_inline)) compute_1_L( double *coeffx, double *coeffy, cField2D *f, int idx, int idy, std::complex<double> *exptheta )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                if( jloc+idy+j_domain_begin_==0 ) {
                    interp_res -= *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
                } else {
                    interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) )*( *exptheta );
                }
            }
        }
        return interp_res;
    };

    inline double __attribute__((always_inline)) compute( double *coeffx, double *coeffy, Field2D *f, int idx, int idy )
    {
        double interp_res( 0. );
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;

                //std::cout<<"f "<<std::fixed << std::setprecision(3)<<(*f)(idx+iloc,idy+jloc)<<std::endl;
            }
        }
        //std::cout<<"interp res "<< interp_res <<std::endl;
        return interp_res;
    };

private:
    inline void coeffs( double xpn, double rpn )
    {
        // Indexes of the central nodes
        ip_ = round( xpn );
        id_ = round( xpn+0.5 );
        jp_ = round( rpn );
        jd_ = round( rpn+0.5 );

        //std::cout << "xpn = " << ip_ << " rpn = " << jp_ << std::endl;

        // Declaration and calculation of the coefficient for interpolation
        double delta2;
        
        deltax_   = xpn - ( double )id_ + 0.5;
        delta2  = deltax_*deltax_;
        coeffxd_[0] = 0.5 * ( delta2-deltax_+0.25 );
        coeffxd_[1] = 0.75 - delta2;
        coeffxd_[2] = 0.5 * ( delta2+deltax_+0.25 );
        
        deltax_   = xpn - ( double )ip_;
        delta2  = deltax_*deltax_;
        coeffxp_[0] = 0.5 * ( delta2-deltax_+0.25 );
        coeffxp_[1] = 0.75 - delta2;
        coeffxp_[2] = 0.5 * ( delta2+deltax_+0.25 );
        
        deltar_   = rpn - ( double )jd_ + 0.5;
        delta2  = deltar_*deltar_;
        coeffyd_[0] = 0.5 * ( delta2-deltar_+0.25 );
        coeffyd_[1] = 0.75 - delta2;
        coeffyd_[2] = 0.5 * ( delta2+deltar_+0.25 );
        
        deltar_   = rpn - ( double )jp_;
        delta2  = deltar_*deltar_;
        coeffyp_[0] = 0.5 * ( delta2-deltar_+0.25 );
        coeffyp_[1] = 0.75 - delta2;
        coeffyp_[2] = 0.5 * ( delta2+deltar_+0.25 );
        
        //std::cout << "coeffxp = " <<  coeffxp_[0] << " " <<  coeffxp_[1] << " " <<  coeffxp_[2] << std::endl;
        //std::cout << "coeffxd = " <<  coeffxd_[0] << " " <<  coeffxd_[1] << " " <<  coeffxd_[2] << std::endl;
        //std::cout << "coeffyp = " <<  coeffyp_[0] << " " <<  coeffyp_[1] << " " <<  coeffyp_[2] << std::endl;
        //std::cout << "coeffyd = " <<  coeffyd_[0] << " " <<  coeffyd_[1] << " " <<  coeffyd_[2] << std::endl;

        // First index for summation
        ip_ = ip_ - i_domain_begin_;
        id_ = id_ - i_domain_begin_;
        jp_ = jp_ - j_domain_begin_;
        jd_ = jd_ - j_domain_begin_;
    };
    
    inline void coeffs( double xpn, double ypn, int* idx_p, int* idx_d,
                        double *coeffxp, double *coeffyp,
                        double *coeffxd, double *coeffyd, double* delta_p )
    {
        // Indexes of the central nodes
        idx_p[0] = round( xpn );
        idx_d[0] = round( xpn+0.5 );
        idx_p[1] = round( ypn );
        idx_d[1] = round( ypn+0.5 );

        // Declaration and calculation of the coefficient for interpolation
        double delta, delta2;

        delta      = xpn - ( double )idx_d[0] + 0.5;
        delta2     = delta*delta;
        coeffxd[0] = 0.5 * ( delta2-delta+0.25 );
        coeffxd[1] = 0.75 - delta2;
        coeffxd[2] = 0.5 * ( delta2+delta+0.25 );

        delta_p[0] = xpn - ( double )idx_p[0];
        delta2     = delta_p[0]*delta_p[0];
        coeffxp[0] = 0.5 * ( delta2-delta_p[0]+0.25 );
        coeffxp[1] = 0.75 - delta2;
        coeffxp[2] = 0.5 * ( delta2+delta_p[0]+0.25 );

        delta      = ypn - ( double )idx_d[1] + 0.5;
        delta2     = delta*delta;
        coeffyd[0] = 0.5 * ( delta2-delta+0.25 );
        coeffyd[1] = 0.75 - delta2;
        coeffyd[2] = 0.5 * ( delta2+delta+0.25 );

        delta_p[1] = ypn - ( double )idx_p[1];
        delta2     = delta_p[1]*delta_p[1];
        coeffyp[0] = 0.5 * ( delta2-delta_p[1]+0.25 );
        coeffyp[1] = 0.75 - delta2;
        coeffyp[2] = 0.5 * ( delta2+delta_p[1]+0.25 );

        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        idx_p[0]   = idx_p[0] - i_domain_begin_;
        idx_d[0]   = idx_d[0] - i_domain_begin_;
        idx_p[1]   = idx_p[1] - j_domain_begin_;
        idx_d[1]   = idx_d[1] - j_domain_begin_;
        
    }

    // Last prim index computed
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Last delta computed
    double deltax_, deltar_ ;
    // exp m theta
    std::complex<double> exp_m_theta_;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3];
    //! Number of modes;
    unsigned int nmodes_;
    
};//END class

#endif
