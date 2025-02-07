#ifndef PROJECTOR3D2ORDERV_H
#define PROJECTOR3D2ORDERV_H

#include "Projector3D.h"
#include "Pragma.h"

class Projector3D2OrderV : public Projector3D
{
public:
    Projector3D2OrderV( Params &, Patch *patch );
    ~Projector3D2OrderV();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    inline void currents( double * __restrict__ Jx,
                          double * __restrict__ Jy,
                          double * __restrict__ Jz,
                          Particles &particles,
                          unsigned int istart,
                          unsigned int iend,
                          double * __restrict__ invgf,
                          int    * __restrict__ iold,
                          double * __restrict__ deltaold,
                          unsigned int buffer_size,
                          int ipart_ref = 0, int bin_shift = 0);

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    inline void currentsAndDensity( double * __restrict__ Jx,
                                    double * __restrict__ Jy,
                                    double * __restrict__ Jz,
                                    double * __restrict__ rho,
                                    Particles &particles,
                                    unsigned int istart,
                                    unsigned int iend,
                                    double * __restrict__ invgf,
                                    int    * __restrict__ iold,
                                    double * __restrict__ deltaold,
                                    unsigned int buffer_size,
                                    int ipart_ref = 0, int bin_shift = 0 );

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int bin, int bin_shift = 0 ) override final;
    
    //! Project global current densities if Ionization in SpeciesV::dynamics,
    void ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion ) override final;

    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell,  int ipart_ref ) override final;
    
    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref ) override;

private:
    double dt, dts2, dts4;

    inline void __attribute__((always_inline)) compute_distances(  double * __restrict__ position_x,
                                                                   double * __restrict__ position_y,
                                    double * __restrict__ position_z,
                                    int npart_total, int ipart, int istart, int ipart_ref,
                                    double *delta0, int *iold, double *Sx0, double *Sy0,
                                    double *Sz0, double *DSx, double *DSy, double *DSz )
    {

        int ipo = iold[0];
        int jpo = iold[1];
        int kpo = iold[2];

        int vecSize = 8;

        double delta = delta0[istart-ipart_ref+ipart];
        double delta2 = delta*delta;

        Sx0[          ipart] = 0.5 * ( delta2-delta+0.25 );
        Sx0[  vecSize+ipart] = 0.75-delta2;
        Sx0[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );
        Sx0[3*vecSize+ipart] = 0.;

        //                            Y                                 //
        delta = delta0[istart-ipart_ref+ipart+npart_total];
        delta2 = delta*delta;

        Sy0[          ipart] = 0.5 * ( delta2-delta+0.25 );
        Sy0[  vecSize+ipart] = 0.75-delta2;
        Sy0[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );
        Sy0[3*vecSize+ipart] = 0.;

        //                            Z                                 //
        delta = delta0[istart-ipart_ref+ipart+2*npart_total];
        delta2 = delta*delta;

        Sz0[          ipart] = 0.5 * ( delta2-delta+0.25 );
        Sz0[  vecSize+ipart] = 0.75-delta2;
        Sz0[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );
        Sz0[3*vecSize+ipart] = 0.;


        // locate the particle on the primal grid at current time-step & calculate coeff. S1
        //                            X                                 //
        double pos = position_x[istart+ipart] * dx_inv_;
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
        DSx [          ipart] = m1 * deltam                             ;
        DSx [  vecSize+ipart] = c0 * deltam + m1 * delta2               -  Sx0[          ipart];
        DSx [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap  -  Sx0[  vecSize+ipart];
        DSx [3*vecSize+ipart] =               p1 * delta2 + c0* deltap  -  Sx0[2*vecSize+ipart];
        DSx [4*vecSize+ipart] =                             p1* deltap  ;
        //                            Y                                 //
        pos = position_y[istart+ipart] * dy_inv_;
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
        DSy [          ipart] = m1 * deltam                            ;
        DSy [  vecSize+ipart] = c0 * deltam + m1 * delta2              -  Sy0[          ipart]                 ;
        DSy [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap -  Sy0[  vecSize+ipart] ;
        DSy [3*vecSize+ipart] =               p1 * delta2 + c0* deltap -  Sy0[2*vecSize+ipart] ;
        DSy [4*vecSize+ipart] =                             p1* deltap  ;
        //                            Z                                 //
        pos = position_z[istart+ipart] * dz_inv_;
        cell = round( pos );
        cell_shift = cell-kpo-k_domain_begin_;
        delta  = pos - ( double )cell;
        delta2 = delta*delta;
        deltam =  0.5 * ( delta2-delta+0.25 );
        deltap =  0.5 * ( delta2+delta+0.25 );
        delta2 = 0.75 - delta2;
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );
        DSz [          ipart] = m1 * deltam                                                            ;
        DSz [  vecSize+ipart] = c0 * deltam + m1 * delta2              -  Sz0[          ipart]                 ;
        DSz [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap -  Sz0[  vecSize+ipart] ;
        DSz [3*vecSize+ipart] =               p1 * delta2 + c0* deltap -  Sz0[2*vecSize+ipart] ;
        DSz [4*vecSize+ipart] =                             p1* deltap  ;

    };


    inline void __attribute__((always_inline)) compute_distances( Particles &particles, int, int ipart, int istart, int, double *, int *iold, double *Sx1, double *Sy1, double *Sz1 )
    {

        int ipo = iold[0];
        int jpo = iold[1];
        int kpo = iold[2];

        int vecSize = 8;

        // locate the particle on the primal grid at current time-step & calculate coeff. S1
        //                            X                                 //
        double pos = particles.position( 0, istart+ipart ) * dx_inv_;
        int cell = round( pos );
        int cell_shift = cell-ipo-i_domain_begin_;
        double delta  = pos - ( double )cell;
        double delta2 = delta*delta;
        double deltam =  0.5 * ( delta2-delta+0.25 );
        double deltap =  0.5 * ( delta2+delta+0.25 );
        delta2 = 0.75 - delta2;
        double m1 = ( cell_shift == -1 );
        double c0 = ( cell_shift ==  0 );
        double p1 = ( cell_shift ==  1 );
        Sx1 [          ipart] = m1 * deltam                            ;
        Sx1 [  vecSize+ipart] = c0 * deltam + m1 * delta2              ;
        Sx1 [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap ;
        Sx1 [3*vecSize+ipart] =               p1 * delta2 + c0* deltap ;
        Sx1 [4*vecSize+ipart] =                             p1* deltap ;
        //                            Y                                 //
        pos = particles.position( 1, istart+ipart ) * dy_inv_;
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
        Sy1 [          ipart] = m1 * deltam                             ;
        Sy1 [  vecSize+ipart] = c0 * deltam + m1 * delta2               ;
        Sy1 [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap  ;
        Sy1 [3*vecSize+ipart] =               p1 * delta2 + c0* deltap  ;
        Sy1 [4*vecSize+ipart] =                             p1* deltap  ;
        //                            Z                                 //
        pos = particles.position( 2, istart+ipart ) * dz_inv_;
        cell = round( pos );
        cell_shift = cell-kpo-k_domain_begin_;
        delta  = pos - ( double )cell;
        delta2 = delta*delta;
        deltam =  0.5 * ( delta2-delta+0.25 );
        deltap =  0.5 * ( delta2+delta+0.25 );
        delta2 = 0.75 - delta2;
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );
        Sz1 [          ipart] = m1 * deltam                             ;
        Sz1 [  vecSize+ipart] = c0 * deltam + m1 * delta2               ;
        Sz1 [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap  ;
        Sz1 [3*vecSize+ipart] =               p1 * delta2 + c0* deltap  ;
        Sz1 [4*vecSize+ipart] =                             p1* deltap  ;

    };

    inline void __attribute__((always_inline)) compute_distances(  double * __restrict__ position_x,
                                    double * __restrict__ position_y,
                                    double * __restrict__ position_z,
                                    int,
                                    int ipart, int istart,
                                    double *, int *iold, double *Sx1, double *Sy1, double *Sz1 )
    {

        int ipo = iold[0];
        int jpo = iold[1];
        int kpo = iold[2];

        int vecSize = 8;

        // locate the particle on the primal grid at current time-step & calculate coeff. S1
        //                            X                                 //
        double pos = position_x[istart+ipart] * dx_inv_;
        int cell = round( pos );
        int cell_shift = cell-ipo-i_domain_begin_;
        double delta  = pos - ( double )cell;
        double delta2 = delta*delta;
        double deltam =  0.5 * ( delta2-delta+0.25 );
        double deltap =  0.5 * ( delta2+delta+0.25 );
        delta2 = 0.75 - delta2;
        double m1 = ( cell_shift == -1 );
        double c0 = ( cell_shift ==  0 );
        double p1 = ( cell_shift ==  1 );
        Sx1 [          ipart] = m1 * deltam                            ;
        Sx1 [  vecSize+ipart] = c0 * deltam + m1 * delta2              ;
        Sx1 [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap ;
        Sx1 [3*vecSize+ipart] =               p1 * delta2 + c0* deltap ;
        Sx1 [4*vecSize+ipart] =                             p1* deltap ;
        //                            Y                                 //
        pos = position_y[istart+ipart] * dy_inv_;
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
        Sy1 [          ipart] = m1 * deltam                             ;
        Sy1 [  vecSize+ipart] = c0 * deltam + m1 * delta2               ;
        Sy1 [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap  ;
        Sy1 [3*vecSize+ipart] =               p1 * delta2 + c0* deltap  ;
        Sy1 [4*vecSize+ipart] =                             p1* deltap  ;
        //                            Z                                 //
        pos = position_z[istart+ipart] * dz_inv_;
        cell = round( pos );
        cell_shift = cell-kpo-k_domain_begin_;
        delta  = pos - ( double )cell;
        delta2 = delta*delta;
        deltam =  0.5 * ( delta2-delta+0.25 );
        deltap =  0.5 * ( delta2+delta+0.25 );
        delta2 = 0.75 - delta2;
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );
        Sz1 [          ipart] = m1 * deltam                             ;
        Sz1 [  vecSize+ipart] = c0 * deltam + m1 * delta2               ;
        Sz1 [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap  ;
        Sz1 [3*vecSize+ipart] =               p1 * delta2 + c0* deltap  ;
        Sz1 [4*vecSize+ipart] =                             p1* deltap  ;

    };

    inline void __attribute__((always_inline)) computeJ( int ipart, double *charge_weight,
                                                        double *DSx, double *DSy, double *DSz,
                                                        double *Sy0, double *Sz0, double *bJx,
                                                        double dxovdt, int nx, int ny, int nz )
    {
        //optrpt complains about the following loop but not unrolling it actually seems to give better result.
        double crx_p = charge_weight[ipart]*dxovdt;

        int vecSize = 8;

        double sum[5];
        sum[0] = 0.;
        UNROLL_S(4)
        for( unsigned int k=1 ; k<5 ; k++ ) {
            sum[k] = sum[k-1]-DSx[( k-1 )*vecSize+ipart];
        }

        double tmp( crx_p * ( one_third*DSy[ipart]*DSz[ipart] ) );
        UNROLL_S(4)
        for( unsigned int i=1 ; i<5 ; i++ ) {
            bJx [( ( i )*nx )*vecSize+ipart] += sum[i]*tmp;
        }

        UNROLL_S(4)
        for( unsigned int k=1 ; k<5 ; k++ ) {
            tmp = crx_p * ( 0.5*DSy[ipart]*Sz0[( k-1 )*vecSize+ipart] + one_third*DSy[ipart]*DSz[k*vecSize+ipart] );
            int index( ( k*nz )*vecSize+ipart );
            UNROLL_S(4)
            for( unsigned int i=1 ; i<5 ; i++ ) {
                bJx [ index+nx*( i )*vecSize ] += sum[i]*tmp;
            }

        }
        UNROLL_S(4)
        for( unsigned int j=1 ; j<5 ; j++ ) {
            tmp = crx_p * ( 0.5*DSz[ipart]*Sy0[( j-1 )*vecSize+ipart] + one_third*DSy[j*vecSize+ipart]*DSz[ipart] );
            int index( ( j*ny )*vecSize+ipart );
            UNROLL_S(4)
            for( unsigned int i=1 ; i<5 ; i++ ) {
                bJx [ index+nx*( i )*vecSize ] += sum[i]*tmp;
            }
        }//i
        UNROLL_S(4)
        for( int j=1 ; j<5 ; j++ ) {
            UNROLL_S(4)
            for( int k=1 ; k<5 ; k++ ) {
                tmp = crx_p * ( Sy0[( j-1 )*vecSize+ipart]*Sz0[( k-1 )*vecSize+ipart]
                                + 0.5*DSy[j*vecSize+ipart]*Sz0[( k-1 )*vecSize+ipart]
                                + 0.5*DSz[k*vecSize+ipart]*Sy0[( j-1 )*vecSize+ipart]
                                + one_third*DSy[j*vecSize+ipart]*DSz[k*vecSize+ipart] );
                int index( ( j*ny + k*nz )*vecSize+ipart );
                UNROLL_S(4)
                for( int i=1 ; i<5 ; i++ ) {
                    bJx [ index+nx*( i )*vecSize ] += sum[i]*tmp;
                }
            }
        }//i
    }

};

#endif
