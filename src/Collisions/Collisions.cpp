
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <fstream>

#include "Collisions.h"

using namespace std;


// Constructor
Collisions::Collisions(
    Params &params,
    double coulomb_log,
    double coulomb_log_factor
) :
    coulomb_log_( coulomb_log ),
    coulomb_log_factor_( coulomb_log_factor ),
    coeff1_( 4.046650232e-21*params.reference_angular_frequency_SI ), // h*omega/(2*me*c^2)
    coeff2_( 2.817940327e-15*params.reference_angular_frequency_SI/299792458. ), // re omega / c
    coeff3_( coeff2_ * coulomb_log_factor_ ),
    coeff4_( pow( 3.*coeff2_, -1./3. ) )
{
}

Collisions::Collisions() :
    coulomb_log_( -1. ),
    coulomb_log_factor_( 0 ),
    coeff1_( 0 ),
    coeff2_( 0 ),
    coeff3_( 0 ),
    coeff4_( 0 )
{
}

void Collisions::prepare()
{
    npairs_tot_  = 0.;
    smean_       = 0.;
    logLmean_    = 0.;
}

#ifdef SMILEI_ACCELERATOR_GPU_OMP
#pragma omp declare target
#endif
void Collisions::apply( Random *random, BinaryProcessData &D, size_t n )
{
    // Buffer intermediate quantities
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        D.buffer1[i] = 1./( D.m[0][i] * D.p_COM[i] );
    }
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        D.buffer2[i] = D.gamma0[i] * D.buffer1[i];
    }
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        D.buffer3[i] = D.R[i] / ( D.p_COM[i] * D.gamma_tot_COM[i] );
    }
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        D.buffer4[i] = ( double ) D.q[0][i] * ( double ) D.q[1][i];
    }
    
    // Calculate coulomb log
    // and the collision parameter s12 (similar to number of real collisions)
    
    // if auto coulomb log
    if( coulomb_log_ <= 0. ) {
    
        // Calculate the minimum impact parameter and the Debye-screening logarithm
        double bmin[SMILEI_BINARYPROCESS_BUFFERSIZE];
        double lnLD[SMILEI_BINARYPROCESS_BUFFERSIZE];
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( size_t i = 0; i<n; i++ ) {
            // Note : 0.00232282 is coeff2 / coeff1
            double b1 = 0.00232282*D.buffer4[i]*D.buffer3[i]*D.buffer2[i];
            if( b1 < 0 ) {
                b1 = -b1;
            }
            bmin[i] = coeff1_ * ( D.buffer1[i] > b1 ? D.buffer1[i] : b1 );
            lnLD[i] = D.debye > 7.839*bmin[i] ? log( D.debye / bmin[i] ) : 2.;
        }
        
        // If no Thomas-Fermi screening
        if( D.screening_group == 0 ) {
            
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( size_t i = 0; i<n; i++ ) {
                double qqqqlogL = D.buffer4[i] * D.buffer4[i] * lnLD[i];
                D.buffer5[i] = coeff3_ * qqqqlogL * D.buffer2[i] * D.buffer2[i] * D.buffer3[i] / ( D.gamma[0][i] * D.gamma[1][i] );
                logLmean_ += lnLD[i];
                smean_    += D.buffer5[i];
            }
        
        // If Thomas-Fermi screening
        } else {
            
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( size_t i = 0; i<n; i++ ) {
                double logL;
                double qqqqlogL; // q1^2 q2^2 logL
                double qqqq = D.buffer4[i] * D.buffer4[i];
                // For e-i collisions, consider the bound-electron screening
                if( D.lTF[i] > 0. ) {
                    double ZZZZ = D.Z1Z2[i] * D.Z1Z2[i];
                    if( D.debye > D.lTF[i] ) {
                        // Bound-electron (Thomas-Fermi) screening
                        double lnLTF = D.lTF[i] > 7.839*bmin[i] ? log( D.lTF[i] / bmin[i] ) : 2.;
                        // Total screening
                        qqqqlogL = qqqq * lnLD[i] + ( ZZZZ - qqqq ) * lnLTF;
                        logL = qqqqlogL / ZZZZ;
                    } else {
                        logL = lnLD[i];
                        qqqqlogL = ZZZZ * lnLD[i];
                    }
                } else {
                    logL = lnLD[i];
                    qqqqlogL = qqqq * lnLD[i];
                }
                
                // Calculate the collision parameter s12 (similar to number of real collisions)
                D.buffer5[i] = coeff3_ * qqqqlogL * D.buffer2[i] * D.buffer2[i] * D.buffer3[i] / ( D.gamma[0][i] * D.gamma[1][i] );
                logLmean_ += logL;
                smean_    += D.buffer5[i];
            }
        
        }
    
    // if constant coulomb log 
    } else {
        
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( size_t i = 0; i<n; i++ ) {
            // Calculate the collision parameter s12 (similar to number of real collisions)
            D.buffer5[i] = coeff3_ * D.buffer4[i] * D.buffer4[i] * coulomb_log_ * D.buffer2[i] * D.buffer2[i] * D.buffer3[i] / ( D.gamma[0][i] * D.gamma[1][i] );
            logLmean_ += coulomb_log_;
            smean_    += D.buffer5[i];
        }
        
    }
    
    // Low-temperature correction to s
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        double n = D.n123 > D.R[i] * D.n223 ? D.n123 : D.R[i] * D.n223;
        double smax = coeff4_ * ( 1 + D.R[i] ) * D.vrel[i] / n;
        
        double &s = D.buffer5[i];
        if( s > smax ) {
            s = smax;
        }
        D.buffer5[i] *= D.dt_correction[i];
    }
    
    // Pick the deflection angles in the center-of-mass frame.
    // Instead of Nanbu http://dx.doi.org/10.1103/PhysRevE.55.4642
    // and Perez http://dx.doi.org/10.1063/1.4742167
    // we made a new fit (faster and more accurate)
    SMILEI_ACCELERATOR_LOOP_SEQ
    for( size_t i = 0; i<n; i++ ) {
        D.buffer3[i] = random->uniform();
        D.buffer4[i] = random->uniform_2pi();
    }
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        double &s = D.buffer5[i];
        double &U1 = D.buffer3[i];
        if( D.buffer5[i] < 4. ) {
            double s2 = s*s;
            double alpha = 0.37*s - 0.005*s2 - 0.0064*s2*s;
            double sin2X2 = alpha * U1 / sqrt( (1.-U1) + alpha*alpha*U1 ); // sin^2( X^2 )
            D.buffer1[i] = 1. - 2.*sin2X2; // cosX
            D.buffer2[i] = 2.*sqrt( sin2X2 *(1.-sin2X2) ); // sinX
        } else {
            D.buffer1[i] = 2.*U1 - 1.; // cosX
            D.buffer2[i] = sqrt( 1. - D.buffer1[i]*D.buffer1[i] ); // sinX
        }
    }
    
    // Calculate the combination with angle phi
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        double &phi = D.buffer4[i];
        D.buffer3[i] = D.buffer2[i]*cos( phi ); // sinXcosPhi
        D.buffer2[i] = D.buffer2[i]*sin( phi ); // sinXsinPhi
    }
    
    // Apply the deflection
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        double & cosX       = D.buffer1[i];
        double & sinXsinPhi = D.buffer2[i];
        double & sinXcosPhi = D.buffer3[i];
        double p_perp = sqrt( D.px_COM[i]*D.px_COM[i] + D.py_COM[i]*D.py_COM[i] );
        double newpx_COM, newpy_COM, newpz_COM;
        if( p_perp > 1.e-10*D.p_COM[i] ) { // make sure p_perp is not too small
            double inv_p_perp = 1./p_perp;
            newpx_COM = ( D.px_COM[i] * D.pz_COM[i] * sinXcosPhi - D.py_COM[i] * D.p_COM[i] * sinXsinPhi ) * inv_p_perp + D.px_COM[i] * cosX;
            newpy_COM = ( D.py_COM[i] * D.pz_COM[i] * sinXcosPhi + D.px_COM[i] * D.p_COM[i] * sinXsinPhi ) * inv_p_perp + D.py_COM[i] * cosX;
            newpz_COM = -p_perp * sinXcosPhi + D.pz_COM[i] * cosX;
        } else { // if p_perp is too small, we use the limit px->0, py=0
            newpx_COM = D.p_COM[i] * sinXcosPhi;
            newpy_COM = D.p_COM[i] * sinXsinPhi;
            newpz_COM = D.p_COM[i] * cosX;
        }
        D.buffer1[i] = newpx_COM;
        D.buffer2[i] = newpy_COM;
        D.buffer3[i] = newpz_COM;
    }
    
    // Go back to the lab frame and update particles momenta
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        double & newpx_COM = D.buffer1[i];
        double & newpy_COM = D.buffer2[i];
        double & newpz_COM = D.buffer3[i];
        double pp = ( D.px_tot[i] * newpx_COM + D.py_tot[i] * newpy_COM + D.pz_tot[i] * newpz_COM ) / ( D.gamma_tot[i] + D.gamma_tot_COM[i] );
        D.buffer4[i] = ( D.gamma_COM0[i] + pp ) / D.gamma_tot_COM[i];
    }
    SMILEI_ACCELERATOR_LOOP_SEQ
    for( size_t i = 0; i<n; i++ ) {
        D.buffer5[i] = random->uniform();
    }
    SMILEI_ACCELERATOR_LOOP_VECTOR
    for( size_t i = 0; i<n; i++ ) {
        double & newpx_COM = D.buffer1[i];
        double & newpy_COM = D.buffer2[i];
        double & newpz_COM = D.buffer3[i];
        double & U2 = D.buffer5[i];
        if( U2 * D.W[0][i] < D.W[1][i] ) { // deflect particle 1 only with some probability
            D.px[0][i] = newpx_COM + D.buffer4[i] * D.px_tot[i];
            D.py[0][i] = newpy_COM + D.buffer4[i] * D.py_tot[i];
            D.pz[0][i] = newpz_COM + D.buffer4[i] * D.pz_tot[i];
        }
        if( U2 * D.W[1][i] < D.W[0][i] ) { // deflect particle 2 only with some probability
            double m12 = D.m[0][i] / D.m[1][i];
            D.px[1][i] = ( -newpx_COM + ( 1 - D.buffer4[i] ) * D.px_tot[i] ) * m12;
            D.py[1][i] = ( -newpy_COM + ( 1 - D.buffer4[i] ) * D.py_tot[i] ) * m12;
            D.pz[1][i] = ( -newpz_COM + ( 1 - D.buffer4[i] ) * D.pz_tot[i] ) * m12;
        }
    }
    
    npairs_tot_ += n;
}
#ifdef SMILEI_ACCELERATOR_GPU_OMP
#pragma omp end declare target
#endif

void Collisions::finish( Params &, Patch *, std::vector<Diagnostic *> &, bool, std::vector<unsigned int>, std::vector<unsigned int>, int )
{
    if( npairs_tot_>0. ) {
        smean_    /= npairs_tot_;
        logLmean_ /= npairs_tot_;
    }
}



