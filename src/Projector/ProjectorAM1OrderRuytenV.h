#ifndef ProjectorAM1OrderRuytenV_H
#define ProjectorAM1OrderRuytenV_H

#include <cmath>
#include "ProjectorAM.h"
#include <complex>
#include "dcomplex.h"
#include "Pragma.h"
#include <iomanip>

using namespace std;

class ProjectorAM1OrderRuytenV : public ProjectorAM
{
public:
    ProjectorAM1OrderRuytenV( Params &, Patch *patch );
    ~ProjectorAM1OrderRuytenV();
    
    //! Project global current densities (EMfields->Jl_/Jr_/Jt_)
    void currents(ElectroMagnAM *emAM, Particles &particles, unsigned int istart, unsigned int iend, double *invgf, int *iold, double *deltaold, std::complex<double> *array_eitheta_old, int npart_total, int ipart_ref = 0, int ispec =0 );

    //! Project global current densities (EMfields->Jl_/Jr_/Jt_/rho), diagFields timestep
    void currentsAndDensity(ElectroMagnAM *emAM, Particles &particles, unsigned int istart, unsigned int iend, double *invgf, int *iold, double *deltaold, std::complex<double> *array_eitheta_old, int npart_total, int ipart_ref = 0, int ispec = 0 );

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basicForComplex( std::complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode ) override final;

    //! Apply boundary conditions on Rho and J
    void axisBC( ElectroMagnAM *emAM, bool diag_flag ) override final;
    void apply_axisBC(std::complex<double> *rhoj,std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, unsigned int imode, bool diag_flag );

    //! Apply boundary conditions on Env_Chi
    void axisBCEnvChi( double *EnvChi ) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion ) override final;
    //
    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell,  int ipart_ref ) override final;
    
    // Project susceptibility
    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref ) override final;
    
private:

    inline void __attribute__((always_inline)) compute_distances(  double * __restrict__ position_x,
                                                                   double * __restrict__ position_y,
                                                                   double * __restrict__ position_z,
                                                                   int * __restrict__,
                                                                   int npart_total, int ipart, int istart, int ipart_ref,
                                                                   double *deltaold, std::complex<double> *array_eitheta_old, int *iold,
                                                                   double *Sl0, double *Sr0, double *DSl, double *DSr,
                                                                   double *r_bar, std::complex<double> *e_bar, std::complex<double> *e_delta_m1)
    {

        int ipo = iold[0];
        int jpo = iold[1];
        int vecSize = 8;

        // locate the particle on the primal grid at former time-step & calculate coeff. S0
        //                            L                                 //
        double delta = deltaold[istart+ipart-ipart_ref];
        double delta2 = delta*delta;
        Sl0[          ipart] = 0.5 * ( delta2-delta+0.25 );
        Sl0[  vecSize+ipart] = 0.75-delta2;
        Sl0[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );
        Sl0[3*vecSize+ipart] = 0.;
        //                            R                                 //
        delta = deltaold[istart+ipart-ipart_ref+npart_total]; // delta = x - x_n si delta > 0 et delta = x - x_n -1 si delta < 0. 
        double pos = (double)jpo + j_domain_begin_ + delta; 
        double x_n = floor(pos);
        double coeff0 = (x_n+1-pos)*(5*x_n + 2 - pos)/(4.*x_n + 2.);
    
        Sr0[          ipart] =      coeff0 * (delta < 0.     );
        Sr0[2*vecSize+ipart] = (1.-coeff0) * (delta >= 0.    );
        Sr0[  vecSize+ipart] = 1. - Sr0[ipart] - Sr0[2*vecSize+ipart];
        Sr0[3*vecSize+ipart] = 0.;


        // locate the particle on the primal grid at current time-step & calculate coeff. S1
        //                            L                                 //
        pos = position_x[istart + ipart] * dl_inv_;
        int cell = round( pos );
        int cell_shift = cell-ipo-i_domain_begin_;
        delta  = pos - ( double )cell;
        delta2 = delta*delta;
        double deltam =  0.5 * ( delta2-delta+0.25 );
        double deltap =  0.5 * ( delta2+delta+0.25 );
        delta2 = 0.75 - delta2;
        double m1 = ( cell_shift == -1 );
        double c0 = ( cell_shift ==  0 );
        double p1 = ( cell_shift ==  1 );
        DSl [          ipart] = m1 * deltam                             ;
        DSl [  vecSize+ipart] = c0 * deltam + m1 * delta2               -  Sl0[          ipart];
        DSl [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap  -  Sl0[  vecSize+ipart];
        DSl [3*vecSize+ipart] =               p1 * delta2 + c0* deltap  -  Sl0[2*vecSize+ipart];
        DSl [4*vecSize+ipart] =                             p1* deltap  ;
        
        double rp = sqrt( position_y[istart+ipart]*position_y[istart+ipart] +  position_z[istart+ipart]*position_z[istart+ipart] );
        pos = rp * dr_inv_;

        x_n = floor(pos);
        coeff0 = (x_n+1-pos)*(5*x_n + 2 - pos)/(4.*x_n + 2.);

        cell = round( pos );
        cell_shift = cell-jpo-j_domain_begin_;
        delta  = pos - ( double )cell;
        deltam =           coeff0 * (double)(delta < 0  );
        deltap =       (1.-coeff0) * (double)(delta >= 0 );
        delta2 =       1. - deltam - deltap;
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );
        DSr [          ipart] = m1 * deltam                            ;
        DSr [  vecSize+ipart] = c0 * deltam + m1 * delta2              -  Sr0[          ipart]                 ;
        DSr [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap -  Sr0[  vecSize+ipart] ;
        DSr [3*vecSize+ipart] =               p1 * delta2 + c0* deltap -  Sr0[2*vecSize+ipart] ;
        DSr [4*vecSize+ipart] =                             p1* deltap  ;

        r_bar[ipart] = ((jpo + j_domain_begin_ + deltaold[istart+ipart-ipart_ref+npart_total])*dr + rp) * 0.5; // r at t = t0 - dt/2
        std::complex<double> eitheta = ( position_y[istart+ipart] + Icpx * position_z[istart+ipart] ) / rp ; //exp(i theta)
        e_delta_m1[ipart] = std::sqrt(eitheta * (2.*std::real(array_eitheta_old[istart+ipart-ipart_ref]) - array_eitheta_old[istart+ipart-ipart_ref]));
        e_bar[ipart] = array_eitheta_old[istart+ipart-ipart_ref] * e_delta_m1[ipart];

    }

    inline void __attribute__((always_inline)) computeJl( int ipart, double *charge_weight, double *DSl, double *DSr, double *Sr0, std::complex<double> *bJ, double, double *invR_local, std::complex<double> *e_bar )
    {

        int vecSize = 8;
        double sum[5];
        std::complex<double> C_m =2.;
        double crl_p = charge_weight[ipart]*dl_ov_dt_;
        
        sum[0] = 0.;
        UNROLL_S(4)
        for( unsigned int k=1 ; k<5 ; k++ ) {
            sum[k] = sum[k-1]-DSl[( k-1 )*vecSize+ipart];
        }
        
        //mode 0
        double tmp = crl_p * ( 0.5*DSr[ipart] ) * invR_local[0] ;
        UNROLL_S(4)
        for( unsigned int i=1 ; i<5 ; i++ ) {
            bJ [( i*5 )*vecSize+ipart] += sum[i] * tmp;
        }
        UNROLL_S(4)
        for ( unsigned int j=1; j<5 ; j++ ) {
            tmp =  crl_p * ( Sr0[(j-1)*vecSize+ipart] + 0.5*DSr[j*vecSize+ipart] ) * invR_local[j];
            UNROLL_S(4)
            for( unsigned int i=1 ; i<5 ; i++ ) {
                bJ [(i*5+j )*vecSize+ipart] += sum[i] * tmp;
            }
        }
        //mode > 0
        for (unsigned int imode=1; imode<Nmode_; imode++){ 
            C_m *= e_bar[ipart];
            tmp = crl_p * ( 0.5*DSr[ipart] ) * invR_local[0] ;
            UNROLL_S(4)
            for( unsigned int i=1 ; i<5 ; i++ ) {
                bJ [200*imode + (i*5 )*vecSize+ipart] += sum[i] * tmp * C_m;
            }
            UNROLL_S(4)
            for ( unsigned int j=1; j<5 ; j++ ) {
                tmp =  crl_p * ( Sr0[(j-1)*vecSize+ipart] + 0.5*DSr[j*vecSize+ipart] ) * invR_local[j];
                UNROLL_S(4)
                for( unsigned int i=1 ; i<5 ; i++ ) {
                    bJ [200*imode + (i*5+j )*vecSize+ipart] += sum[i] * tmp * C_m;
                }
            }
        }
    }

    inline void __attribute__((always_inline)) computeJr( int ipart, double *charge_weight, double *DSl, double *DSr, double *Sl0, std::complex<double> *bJ, double one_ov_dt, double *invRd_local, std::complex<double> *e_bar, int jpo )
    {

        int vecSize = 8;
        double sum[5];
        std::complex<double> C_m =2.;
        double crr_p = charge_weight[ipart]*one_ov_dt;
        
        sum[4] = 0.;
        UNROLL_S(4)
        for( int j=3 ; j>=0 ; j-- ) {
            sum[j] = sum[j+1] * abs( jpo+j+1 + j_domain_begin_ + 0.5 )*dr * invRd_local[j+1] +  crr_p * DSr[(j+1)*vecSize+ipart] * invRd_local[j+1]*dr;
        }
        
        //mode 0
        double tmp = 0.5*DSl[ipart];
        UNROLL_S(4)
        for( unsigned int j=0 ; j<4 ; j++ ) {
            bJ [( j+1 )*vecSize+ipart] += sum[j] * tmp;
        }
        UNROLL_S(4)
        for ( unsigned int i=1; i<5 ; i++ ) {
            tmp = Sl0[(i-1)*vecSize + ipart] + 0.5*DSl[i*vecSize + ipart];
            UNROLL_S(4)
            for( unsigned int j=0 ; j<4 ; j++ ) {
                bJ [(i*5+j+1 )*vecSize + ipart] += sum[j] * tmp;
            }
        }

        //mode > 0
        for (unsigned int imode=1; imode<Nmode_; imode++){ 
            C_m *= e_bar[ipart];
            tmp = 0.5*DSl[ipart];
            UNROLL_S(4)
            for( unsigned int j=0 ; j<4 ; j++ ) {
                bJ [200*imode + ( j+1 )*vecSize+ipart] += sum[j] * tmp * C_m;
            }
            UNROLL_S(4)
            for ( unsigned int i=1; i<5 ; i++ ) {
                tmp = Sl0[(i-1)*vecSize + ipart] + 0.5*DSl[i*vecSize + ipart];
                UNROLL_S(4)
                for( unsigned int j=0 ; j<4 ; j++ ) {
                    bJ [200*imode + (i*5+j+1 )*vecSize + ipart] += sum[j] * tmp * C_m;
                }
            }
        }
    }
 
    inline void __attribute__((always_inline)) computeJt( int ipart, 
                                                                   double * __restrict__ momentum_y,
                                                                   double * __restrict__ momentum_z,
                                                                   double *charge_weight,
                                                                   double *invgf,
                                                                   double *DSl, double *DSr, double *Sl0_buff_vect, double *Sr0_buff_vect,
                                                                   std::complex<double> *bJ, double *invR_local, double *r_bar, std::complex<double> *e_bar_m1, std::complex<double> *e_delta_m1, 
                                                                   double one_ov_dt)
    {

        int vecSize = 8;
        //mode 0
        //i=0 and i=4: Sl0=Sr0=0
        //Sr0 and Sr1 need to be divided by inv_R
        // S1 = DS + S0
        double Sl1[5], Sr1[5];

        Sl1[0] = DSl[            ipart] ;
        Sl1[1] = DSl[1*vecSize + ipart]  + Sl0_buff_vect[            ipart];
        Sl1[2] = DSl[2*vecSize + ipart]  + Sl0_buff_vect[1*vecSize + ipart];
        Sl1[3] = DSl[3*vecSize + ipart]  + Sl0_buff_vect[2*vecSize + ipart];
        Sl1[4] = DSl[4*vecSize + ipart] ; 

        Sr1[0] = (DSr[            ipart]                                    ) * invR_local[0];
        Sr1[1] = (DSr[1*vecSize + ipart]  + Sr0_buff_vect[            ipart]) * invR_local[1];
        Sr1[2] = (DSr[2*vecSize + ipart]  + Sr0_buff_vect[1*vecSize + ipart]) * invR_local[2];
        Sr1[3] = (DSr[3*vecSize + ipart]  + Sr0_buff_vect[2*vecSize + ipart]) * invR_local[3];
        Sr1[4] = (DSr[4*vecSize + ipart]                                    ) * invR_local[4]; 

        Sr0_buff_vect[            ipart] *= invR_local[1];
        Sr0_buff_vect[1*vecSize + ipart] *= invR_local[2];
        Sr0_buff_vect[2*vecSize + ipart] *= invR_local[3];

        //mode 0
        std::complex<double> crt_p= charge_weight[ipart]*( momentum_z[ipart]* real(e_bar_m1[ipart]) - momentum_y[ipart]*imag(e_bar_m1[ipart]) ) * invgf[ipart];
        std::complex<double> e_delta = 0.5;
        std::complex<double> e_bar = 1.;


        for (unsigned int imode=0; imode<Nmode_; imode++){ 
            if (imode > 0){
                e_delta *= e_delta_m1[ipart];
                e_bar *= e_bar_m1[ipart];
                crt_p = charge_weight[ipart]*Icpx*e_bar * one_ov_dt * 2. * r_bar[ipart] /( double )imode ;
            }

            //j=0 case
            UNROLL_S(5)
            for( unsigned int i=0 ; i<5 ; i++ ) {
                bJ [200*imode + (i*5 )*vecSize + ipart] -= crt_p*(Sr1[0]*Sl1[i]*( e_delta-1. ) );
            }
            //i=0 case
            UNROLL_S(4)
            for( unsigned int j=1 ; j<5 ; j++ ) {
                bJ [200*imode + (j)*vecSize + ipart] -= crt_p*(Sr1[j]*Sl1[0]*( e_delta-1. ) );
            }


            UNROLL_S(4)
            for( unsigned int i=1 ; i<5 ; i++ ) {
                UNROLL_S(4)
                for ( unsigned int j=1; j<5 ; j++ ) {
                    bJ [200*imode + (i*5+j )*vecSize + ipart] -= crt_p*(Sr1[j]*Sl1[i]*( e_delta-1. ) - Sr0_buff_vect[(j-1)*vecSize + ipart]*Sl0_buff_vect[(i-1)*vecSize + ipart]*(std::conj(e_delta) - 1.));
                }
            }

            if (imode == 0) e_delta = 1. ; //Restore e_delta correct initial value.
       }
    }

    inline void __attribute__((always_inline)) computeRho( int ipart, 
                                                                   double *charge_weight,
                                                                   double *DSl, double *DSr, double *Sl0_buff_vect, double *Sr0_buff_vect,
                                                                   std::complex<double> *brho, double *invR_local, std::complex<double> *e_bar_m1) 
    {

        int vecSize = 8;
        //mode 0
        //i=0 and i=4: Sl0=Sr0=0
        //Sr0 and Sr1 need to be divided by inv_R
        // S1 = DS + S0
        double Sl1[5], Sr1[5];

        Sl1[0] = DSl[            ipart] ;
        Sl1[1] = DSl[1*vecSize + ipart]  + Sl0_buff_vect[            ipart];
        Sl1[2] = DSl[2*vecSize + ipart]  + Sl0_buff_vect[1*vecSize + ipart];
        Sl1[3] = DSl[3*vecSize + ipart]  + Sl0_buff_vect[2*vecSize + ipart];
        Sl1[4] = DSl[4*vecSize + ipart] ; 

        Sr1[0] = (DSr[            ipart]                                    ) * invR_local[0];
        Sr1[1] = (DSr[1*vecSize + ipart]  + Sr0_buff_vect[            ipart]) * invR_local[1];
        Sr1[2] = (DSr[2*vecSize + ipart]  + Sr0_buff_vect[1*vecSize + ipart]) * invR_local[2];
        Sr1[3] = (DSr[3*vecSize + ipart]  + Sr0_buff_vect[2*vecSize + ipart]) * invR_local[3];
        Sr1[4] = (DSr[4*vecSize + ipart]                                    ) * invR_local[4]; 

        Sr0_buff_vect[            ipart] *= invR_local[1];
        Sr0_buff_vect[1*vecSize + ipart] *= invR_local[2];
        Sr0_buff_vect[2*vecSize + ipart] *= invR_local[3];

        //mode 0
        std::complex<double> C_m = 1.;
        std::complex<double> e_bar = 1.;

        for (unsigned int imode=0; imode<Nmode_; imode++){ 
            if (imode > 0){
                e_bar *= e_bar_m1[ipart];
                C_m = 2. * e_bar;
            }

            UNROLL_S(5)
            for( unsigned int i=0 ; i<5 ; i++ ) {
                UNROLL_S(5)
                for ( unsigned int j=0; j<5 ; j++ ) {
                    brho [200*imode + (i*5+j )*vecSize + ipart] += C_m * charge_weight[ipart]*(Sr1[j]*Sl1[i]);
                }
            }
       }
    }

};

#endif

