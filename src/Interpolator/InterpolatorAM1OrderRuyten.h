#ifndef INTERPOLATORAM1ORDERRUYTEN_H
#define INTERPOLATORAM1ORDERRUYTEN_H


#include "InterpolatorAM.h"
#include "cField2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 1st order interpolator with Ruyten correction for AM simulations based on the article from W. Ruyten: Density-Conserving Shape Factors for Particle Simulations in Cylindrical and Spherical Coordinates J. Comp. Phys. 105, 224-232 (1993)
//  --------------------------------------------------------------------------------------------------------------------
class InterpolatorAM1OrderRuyten final : public InterpolatorAM
{
public:
    InterpolatorAM1OrderRuyten( Params &, Patch * );
    ~InterpolatorAM1OrderRuyten() override final {};
    
    inline void __attribute__((always_inline)) fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell = 0, int ipart_ref = 0 ) override final ;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;

    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;
    void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;

    inline std::complex<double> __attribute__((always_inline)) compute( double *coeffx, double *coeffy, cField2D *f, int idx, int idy )
    {
        std::complex<double> interp_res( 0. );
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=0 ; jloc<2 ; jloc++ ) {
                //std::cout << "Er standard idx = " << idx << " idy = " << idy << " iloc = " << iloc << " jloc = " << jloc << " " <<  *( coeffx+iloc ) << " " <<  *( coeffy+jloc ) << " " <<  ( *f )( idx+iloc, idy+jloc )  << std::endl;
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;

                //std::cout<<"f "<<std::fixed << std::setprecision(3)<<(*f)(idx+iloc,idy+jloc)<<std::endl;
            }
        }
        //std::cout<<"interp res "<< interp_res <<std::endl;
        return interp_res;
    };

    inline double __attribute__((always_inline)) compute( double *coeffx, double *coeffy, Field2D *f, int idx, int idy )
    {
        double interp_res( 0. );
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=0 ; jloc<2 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( ( *f )( idx+iloc, idy+jloc ) ) ;

                //std::cout<<"f "<<std::fixed << std::setprecision(3)<<(*f)(idx+iloc,idy+jloc)<<std::endl;
            }
        }
        //std::cout<<"interp res "<< interp_res <<std::endl;
        return interp_res;
    };

private:
    
    inline void coeffs( double xpn, double ypn, int* idx_p, int* idx_d,
                        double *coeffxp, double *coeffyp,
                        double *coeffxd, double *coeffyd, double* delta_p) const
    {
        // Indexes of the central nodes
        idx_p[0] = round( xpn );
        idx_p[1] = floor( ypn );

        // Declaration and calculation of the coefficient for interpolation
        double delta, delta2;

        delta_p[0] = xpn - ( double )idx_p[0];
        delta2     = delta_p[0]*delta_p[0];
        coeffxp[0] = 0.5 * ( delta2-delta_p[0]+0.25 );
        coeffxp[1] = 0.75 - delta2;
        coeffxp[2] = 0.5 * ( delta2+delta_p[0]+0.25 );


        delta_p[1] = ypn - ( double )idx_p[1]; // x - x_n
        double x_n = idx_p[1];
        //Ruyten scheme
        coeffyp[0] = (1. - delta_p[1])*(5*x_n + 2 - ypn)/(4.*x_n + 2.);
        coeffyp[1] = 1. - coeffyp[0]; //conservation de la charge

        // First index for summation
        idx_p[0]   = idx_p[0] - i_domain_begin_;
        idx_p[1]   = idx_p[1] - j_domain_begin_;

        if (idx_d){
            idx_d[0] = round( xpn+0.5 );
            idx_d[1] = floor( ypn+0.5 );

            delta      = xpn - ( double )idx_d[0] + 0.5;
            delta2     = delta*delta;
            coeffxd[0] = 0.5 * ( delta2-delta+0.25 );
            coeffxd[1] = 0.75 - delta2;
            coeffxd[2] = 0.5 * ( delta2+delta+0.25 );

            delta      = ypn - ( double )idx_d[1] + 0.5;
            double x_np1 = ( double )idx_d[1] + 0.5 ; 
            x_n = ( (( double )idx_d[1] - 0.5 > 0.) ? ( double )idx_d[1] - 0.5 : 0. ); 
            //Ruyten scheme
            coeffyd[0] = 0.5*(1. - delta)*(5*x_n + 2 - ypn)/(x_np1*x_np1-x_n*x_n);
            coeffyd[1] = 1. - coeffyd[0]; //conservation de la charge

            idx_d[0]   = idx_d[0] - i_domain_begin_;
            idx_d[1]   = idx_d[1] - j_domain_begin_;
        }
    }

    // exp m theta
    std::complex<double> exp_m_theta_;
    //! Number of modes;
    unsigned int nmodes_;
    
};//END class

#endif
