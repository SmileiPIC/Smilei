#ifndef PROJECTORAM2ORDERV_H
#define PROJECTORAM2ORDERV_H

#include "ProjectorAM.h"


class ProjectorAM2OrderV : public ProjectorAM
{
public:
    ProjectorAM2OrderV( Params &, Patch *patch );
    ~ProjectorAM2OrderV();
    
    //! Project global current densities (EMfields->Jl_/Jr_/Jt_)
    void currents(ElectroMagnAM *emAM, Particles &particles, unsigned int istart, unsigned int iend, std::vector<double> *invgf, int *iold, double *deltaold, std::complex<double> *array_eitheta_old, int ipart_ref = 0 );
    ////! Project global current densities (EMfields->Jl_/Jr_/Jt_/rho), diagFields timestep
    //inline void currentsAndDensity( double *Jl, double *Jr, double *Jt, double *rho, Particles &particles, unsigned int istart, unsigned int iend, std::vector<double> *invgf, int *iold, double *deltaold, int ipart_ref );
    //
    ////! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    //void basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int bin ) override final;
    //
    ////! Project global current densities if Ionization in Species::dynamics,
    //void ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion ) override final;
    //
    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell,  int ipart_ref ) override final;
    //
    // Project susceptibility
    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref ) override final;
    
private:

    inline void compute_distances( Particles &particles, int npart_total, int ipart, int istart, int ipart_ref, double *deltaold, int *iold, double *Sl0, double *Sr0, double *DSl, double *DSr, std::complex<double> *e_bar)
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
        delta = deltaold[istart+ipart-ipart_ref+npart_total];
        delta2 = delta*delta;
        Sr0[          ipart] = 0.5 * ( delta2-delta+0.25 );
        Sr0[  vecSize+ipart] = 0.75-delta2;
        Sr0[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );
        Sr0[3*vecSize+ipart] = 0.;


        // locate the particle on the primal grid at current time-step & calculate coeff. S1
        //                            L                                 //
        double pos = particles.position( 0, istart + ipart ) * dl_inv_;
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
        
        //                            R                                 //
        double rp = sqrt( particles.position( 1, istart+ipart )*particles.position( 1, istart+ipart )+particles.position( 2, istart+ipart )*particles.position( 2, istart+ipart ) );
        pos = rp * dr_inv_;
        cell = round( pos );
        cell_shift = cell-jpo-j_domain_begin_;
        delta  = pos - ( double )cell;
        delta2 = delta*delta;
        deltam =  0.5 * ( delta2-delta+0.25 );
        deltap =  0.5 * ( delta2+delta+0.25 );
        delta2 = 0.75 - delta2;
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );
        DSr [          ipart] = m1 * deltam                            ;
        DSr [  vecSize+ipart] = c0 * deltam + m1 * delta2              -  Sr0[          ipart]                 ;
        DSr [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap -  Sr0[  vecSize+ipart] ;
        DSr [3*vecSize+ipart] =               p1 * delta2 + c0* deltap -  Sr0[2*vecSize+ipart] ;
        DSr [4*vecSize+ipart] =                             p1* deltap  ;
    }

    inline void computeJl( int ipart, double *charge_weight, double *DSl, double *DSr, double *Sr0, std::complex<double> *bJ, double dl_ov_dt, double *invR_local, std::complex<double> *C_m, std::complex<double> *e_bar )
    {

        int vecSize = 8;
        double sum[5];

        double crl_p = charge_weight[ipart]*dl_ov_dt_;
        
        sum[0] = 0.;
        for( unsigned int k=1 ; k<5 ; k++ ) {
            sum[k] = sum[k-1]-DSl[( k-1 )*vecSize+ipart];
        }
        
        //mode 0
        double tmp( crl_p * ( 0.5*DSr[ipart] ) * invR_local[0] );
        for( unsigned int i=1 ; i<5 ; i++ ) {
            bJ [( i*5 )*vecSize+ipart] += sum[i] * tmp;
        }
        for ( unsigned int j=1; j<5 ; j++ ) {
            tmp =  crl_p * ( Sr0[(j-1)*vecSize+ipart] + 0.5*DSr[j*vecSize+ipart] ) * invR_local[j];
            for( unsigned int i=1 ; i<5 ; i++ ) {
                bJ [(i*5+j )*vecSize+ipart] += sum[i] * tmp;
            }
        }
        //mode > 0
        for (unsigned int imode=1; imode<Nmode_; imode++){ 
            C_m[ipart] *= e_bar[ipart];
            double tmp( crl_p * ( 0.5*DSr[ipart] ) * invR_local[0] );
            for( unsigned int i=1 ; i<5 ; i++ ) {
                bJ [(200*imode + i*5 )*vecSize+ipart] += sum[i] * tmp * C_m[ipart];
            }
            for ( unsigned int j=1; j<5 ; j++ ) {
                tmp =  crl_p * ( Sr0[(j-1)*vecSize+ipart] + 0.5*DSr[j*vecSize+ipart] ) * invR_local[j];
                for( unsigned int i=1 ; i<5 ; i++ ) {
                    bJ [(200*imode + i*5+j )*vecSize+ipart] += sum[i] * tmp * C_m[ipart];
                }
            }
        }
    }

};

#endif

