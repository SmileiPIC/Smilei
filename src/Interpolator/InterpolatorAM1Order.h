#ifndef INTERPOLATORAM1ORDER_H
#define INTERPOLATORAM1ORDER_H


#include "InterpolatorAM.h"
#include "cField2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 1st order interpolator for AM simulations
//  --------------------------------------------------------------------------------------------------------------------
class InterpolatorAM1Order : public InterpolatorAM
{

public:
    InterpolatorAM1Order( Params &, Patch * );
    ~InterpolatorAM1Order() override final {};
    
    inline void fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final ;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;
    
    
    inline std::complex<double> compute( double *coeffx, double *coeffy, cField2D *f, int idx, int idy )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=0 ; iloc<2 ; iloc++ ) {
            for( int jloc=0 ; jloc<2 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
                
                //std::cout<<"f "<<std::fixed << std::setprecision(3)<<(*f)(idx+iloc,idy+jloc)<<std::endl;
            }
        }
        //std::cout<<"interp res "<< interp_res <<std::endl;
        return interp_res;
    };
     
    inline std::complex<double> compute_0_T( double *coeffx, double *coeffy, cField2D *f, int idx, int idy )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=0 ; iloc<2 ; iloc++ ) {
            for( int jloc=0 ; jloc<2 ; jloc++ ) {
                //if( jloc+idy+j_domain_begin==0 ) {
                //    interp_res -= *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
                //} else {
                    interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
                //}
            }
        }
        return interp_res;
    };
    inline std::complex<double> compute_0_L( double *coeffx, double *coeffy, cField2D *f, int idx, int idy )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=0 ; iloc<2 ; iloc++ ) {
            for( int jloc=0 ; jloc<2 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
            }
        }
        return interp_res;
    };
    inline std::complex<double> compute_1_T( double *coeffx, double *coeffy, cField2D *f, int idx, int idy, std::complex<double> *exptheta )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=0 ; iloc<2 ; iloc++ ) {
            for( int jloc=0 ; jloc<2 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) )*( *exptheta ) ;
            }
        }
        return interp_res;
    };
    inline std::complex<double> compute_1_L( double *coeffx, double *coeffy, cField2D *f, int idx, int idy, std::complex<double> *exptheta )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=0 ; iloc<2 ; iloc++ ) {
            for( int jloc=0 ; jloc<2 ; jloc++ ) {
                //if( jloc+idy+j_domain_begin==0 ) {
                //    interp_res -= *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;
                //} else {
                    interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) )*( *exptheta );
                //}
            }
        }
        return interp_res;
    };
private:
    inline void coeffs( double xpn, double rpn )
    {
        // Indexes of the central nodes
        ip_ = floor( xpn );
        jp_ = floor( rpn );
        
        // Declaration and calculation of the coefficient for interpolation
        
        deltax   = xpn - ( double )ip_;
        coeffxp_[0] = 1. - deltax;
        coeffxp_[1] = deltax;
        
        deltar   = rpn - ( double )jp_;
        coeffyp_[0] = 1. - deltar;
        coeffyp_[1] = deltar;
        coeffyp_[2] = coeffyp_[0];
        coeffyp_[3] = coeffyp_[1];

        if (rpn < 0.){ // If particle is between 0 and dr/2 initial jp_=-1
            jp_ = 0;
            // coeffs 2-3 are used when F(-dr/2) = - F(dr/2) <==> when field mode is zero on axis
            coeffyp_[2] = coeffyp_[1] - coeffyp_[0];    
            coeffyp_[3] = 0.; // Terms are already acuumulated in coeffyp_[2] 
            // coeffs 0-1 are used when F(-dr/2) = + F(dr/2) <==> when field is constant on axis
            coeffyp_[0] = 1.; // = coeffyp_[1] + coeffyp_[0];    
            coeffyp_[1] = 0.; // Terms are already acuumulated in coeffyp_[0] 
            deltar -= 1.; // To account for the cell shift and proper recomputation of r_old in projector
        }
        
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
    };
    
    // Last prim index computed
    int ip_, jp_;
    // Last delta computed
    double deltax, deltar ;
    // exp m theta
    std::complex<double> exp_m_theta;
    // Interpolation coefficient on Prim grid
    double coeffxp_[2], coeffyp_[4] ;// coeff[2-3] to use when fields are 0 on axis.
    //! Number of modes;
    unsigned int nmodes;
    
    
};//END class

#endif
