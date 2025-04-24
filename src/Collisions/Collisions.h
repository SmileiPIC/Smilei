#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <vector>
#include <cmath>

#include "Tools.h"
#include "H5.h"
#include "Random.h"
#include "Params.h"
#include "BinaryProcessData.h"

class Patch;
class Species;
class VectorPatch;

class Collisions
{
public:
    //! Constructor for Collisions between two species
    Collisions(
        Params &params,
        double coulomb_log,
        double coulomb_log_factor
    );
    //! Constructor for Collisions that do nothing (no collisions)
    Collisions();
    
    operator bool() const { 
        return coulomb_log_ >= 0.;
    }
    
    void prepare();
    
    #ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine vector
    #endif
    void apply( Random *random, BinaryProcessData &D, uint32_t n )
    {
        // Buffer intermediate quantities
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            D.buffer1[i] = 1./( D.m[0][i] * D.p_COM[i] );
        }
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            D.buffer2[i] = D.gamma0[i] * D.buffer1[i];
        }
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            D.buffer3[i] = D.R[i] / ( D.p_COM[i] * D.gamma_tot_COM[i] );
        }
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            D.buffer4[i] = ( SMILEI_BINARYPROCESS_FLOAT ) D.q[0][i] * ( SMILEI_BINARYPROCESS_FLOAT ) D.q[1][i];
        }
        
        // Calculate coulomb log
        // and the collision parameter s12 (similar to number of real collisions)
        double logLmean = 0., smean = 0.;
        
        // if auto coulomb log
        if( coulomb_log_ <= 0. ) {
        
            // Calculate the minimum impact parameter and the Debye-screening logarithm
            SMILEI_BINARYPROCESS_FLOAT bmin[SMILEI_BINARYPROCESS_BUFFERSIZE];
            SMILEI_BINARYPROCESS_FLOAT lnLD[SMILEI_BINARYPROCESS_BUFFERSIZE];
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                // Note : 0.00232282 is coeff2 / coeff1
                SMILEI_BINARYPROCESS_FLOAT b1 = 0.00232282*D.buffer4[i]*D.buffer3[i]*D.buffer2[i];
                if( b1 < 0 ) {
                    b1 = -b1;
                }
                bmin[i] = coeff1_ * ( D.buffer1[i] > b1 ? D.buffer1[i] : b1 );
                lnLD[i] = D.debye > 7.839*bmin[i] ? log( D.debye / bmin[i] ) : 2.;
            }
            
            // If no Thomas-Fermi screening
            if( D.screening_group == 0 ) {
                
                
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    SMILEI_BINARYPROCESS_FLOAT qqqqlogL = D.buffer4[i] * D.buffer4[i] * lnLD[i];
                    D.buffer5[i] = coeff3_ * qqqqlogL * D.buffer2[i] * D.buffer2[i] * D.buffer3[i] / ( D.gamma[0][i] * D.gamma[1][i] );
                    logLmean += lnLD[i];
                }
            
            // If Thomas-Fermi screening
            } else {
                
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    SMILEI_BINARYPROCESS_FLOAT logL;
                    SMILEI_BINARYPROCESS_FLOAT qqqqlogL; // q1^2 q2^2 logL
                    SMILEI_BINARYPROCESS_FLOAT qqqq = D.buffer4[i] * D.buffer4[i];
                    // For e-i collisions, consider the bound-electron screening
                    if( D.lTF[i] > 0. ) {
                        SMILEI_BINARYPROCESS_FLOAT ZZZZ = D.Z1Z2[i] * D.Z1Z2[i];
                        if( D.debye > D.lTF[i] ) {
                            // Bound-electron (Thomas-Fermi) screening
                            SMILEI_BINARYPROCESS_FLOAT lnLTF = D.lTF[i] > 7.839*bmin[i] ? log( D.lTF[i] / bmin[i] ) : 2.;
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
                    logLmean += logL;
                }
            
            }
        
        // if constant coulomb log 
        } else {
            
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                // Calculate the collision parameter s12 (similar to number of real collisions)
                D.buffer5[i] = coeff3_ * D.buffer4[i] * D.buffer4[i] * coulomb_log_ * D.buffer2[i] * D.buffer2[i] * D.buffer3[i] / ( D.gamma[0][i] * D.gamma[1][i] );
                logLmean += coulomb_log_;
            }
            
        }
        SMILEI_ACCELERATOR_ATOMIC
        logLmean_ += logLmean;
        
        // Low-temperature correction to s
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            SMILEI_BINARYPROCESS_FLOAT n = D.n123 > D.R[i] * D.n223 ? D.n123 : D.R[i] * D.n223;
            SMILEI_BINARYPROCESS_FLOAT smax = coeff4_ * ( 1 + D.R[i] ) * D.vrel[i] / n;
            
            SMILEI_BINARYPROCESS_FLOAT &s = D.buffer5[i];
            if( s > smax ) {
                s = smax;
            }
            smean += D.buffer5[i];
            D.buffer5[i] *= D.dt_correction[i];
        }
        SMILEI_ACCELERATOR_ATOMIC
        smean_ += smean;
        
        // Pick the deflection angles in the center-of-mass frame.
        // Instead of Nanbu http://dx.doi.org/10.1103/PhysRevE.55.4642
        // and Perez http://dx.doi.org/10.1063/1.4742167
        // we made a new fit (faster and more accurate)
        SMILEI_ACCELERATOR_LOOP_SEQ
        for( uint32_t i = 0; i<n; i++ ) {
            D.buffer3[i] = random->uniform1f();
            D.buffer4[i] = random->uniform_2pi();
        }
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            SMILEI_BINARYPROCESS_FLOAT &s = D.buffer5[i];
            SMILEI_BINARYPROCESS_FLOAT &U1 = D.buffer3[i];
            if( D.buffer5[i] < 4. ) {
                SMILEI_BINARYPROCESS_FLOAT s2 = s*s;
                SMILEI_BINARYPROCESS_FLOAT alpha = 0.37*s - 0.005*s2 - 0.0064*s2*s;
                SMILEI_BINARYPROCESS_FLOAT sin2X2 = alpha * U1 / sqrt( (1.-U1) + alpha*alpha*U1 ); // sin^2( X^2 )
                D.buffer1[i] = 1. - 2.*sin2X2; // cosX
                D.buffer2[i] = 2.*sqrt( sin2X2 *(1.-sin2X2) ); // sinX
            } else {
                D.buffer1[i] = 2.*U1 - 1.; // cosX
                D.buffer2[i] = sqrt( 1. - D.buffer1[i]*D.buffer1[i] ); // sinX
            }
        }
        
        // Calculate the combination with angle phi
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            SMILEI_BINARYPROCESS_FLOAT &phi = D.buffer4[i];
            D.buffer3[i] = D.buffer2[i]*cos( phi ); // sinXcosPhi
            D.buffer2[i] = D.buffer2[i]*sin( phi ); // sinXsinPhi
        }
        
        // Apply the deflection
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            SMILEI_BINARYPROCESS_FLOAT & cosX       = D.buffer1[i];
            SMILEI_BINARYPROCESS_FLOAT & sinXsinPhi = D.buffer2[i];
            SMILEI_BINARYPROCESS_FLOAT & sinXcosPhi = D.buffer3[i];
            SMILEI_BINARYPROCESS_FLOAT p_perp = sqrt( D.px_COM[i]*D.px_COM[i] + D.py_COM[i]*D.py_COM[i] );
            SMILEI_BINARYPROCESS_FLOAT newpx_COM, newpy_COM, newpz_COM;
            if( p_perp > 1.e-10*D.p_COM[i] ) { // make sure p_perp is not too small
                SMILEI_BINARYPROCESS_FLOAT inv_p_perp = 1./p_perp;
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
        for( uint32_t i = 0; i<n; i++ ) {
            SMILEI_BINARYPROCESS_FLOAT & newpx_COM = D.buffer1[i];
            SMILEI_BINARYPROCESS_FLOAT & newpy_COM = D.buffer2[i];
            SMILEI_BINARYPROCESS_FLOAT & newpz_COM = D.buffer3[i];
            SMILEI_BINARYPROCESS_FLOAT pp = ( D.px_tot[i] * newpx_COM + D.py_tot[i] * newpy_COM + D.pz_tot[i] * newpz_COM ) / ( D.gamma_tot[i] + D.gamma_tot_COM[i] );
            D.buffer4[i] = ( D.gamma_COM0[i] + pp ) / D.gamma_tot_COM[i];
        }
        SMILEI_ACCELERATOR_LOOP_SEQ
        for( uint32_t i = 0; i<n; i++ ) {
            D.buffer5[i] = random->uniform();
        }
        SMILEI_ACCELERATOR_LOOP_VECTOR
        for( uint32_t i = 0; i<n; i++ ) {
            SMILEI_BINARYPROCESS_FLOAT & newpx_COM = D.buffer1[i];
            SMILEI_BINARYPROCESS_FLOAT & newpy_COM = D.buffer2[i];
            SMILEI_BINARYPROCESS_FLOAT & newpz_COM = D.buffer3[i];
            SMILEI_BINARYPROCESS_FLOAT & U2 = D.buffer5[i];
            if( U2 * D.W[0][i] < D.W[1][i] ) { // deflect particle 1 only with some probability
                D.px[0][i] = newpx_COM + D.buffer4[i] * D.px_tot[i];
                D.py[0][i] = newpy_COM + D.buffer4[i] * D.py_tot[i];
                D.pz[0][i] = newpz_COM + D.buffer4[i] * D.pz_tot[i];
            }
            if( U2 * D.W[1][i] < D.W[0][i] ) { // deflect particle 2 only with some probability
                SMILEI_BINARYPROCESS_FLOAT m12 = D.m[0][i] / D.m[1][i];
                D.px[1][i] = ( -newpx_COM + ( 1 - D.buffer4[i] ) * D.px_tot[i] ) * m12;
                D.py[1][i] = ( -newpy_COM + ( 1 - D.buffer4[i] ) * D.py_tot[i] ) * m12;
                D.pz[1][i] = ( -newpz_COM + ( 1 - D.buffer4[i] ) * D.pz_tot[i] ) * m12;
            }
        }
        
        SMILEI_ACCELERATOR_ATOMIC
        npairs_tot_ += (unsigned int) n;
    }
    
    void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool intra, std::vector<unsigned int> sg1, std::vector<unsigned int> sg2, int itime );
    
    std::string name() {
        std::ostringstream t;
        t << "Collisions with Coulomb logarithm: ";
        if( coulomb_log_ == 0. ) {
            t << "auto";
        } else {
            t << coulomb_log_;
        }
        if( coulomb_log_factor_ != 1. ) {
            t << " x " << coulomb_log_factor_;
        }
        return t.str();
    };
    
    unsigned int npairs_tot_;
    double smean_, logLmean_;
    
protected:
    
    //! Coulomb logarithm (zero or negative means automatic)
    const double coulomb_log_;
    
    //! Coulomb logarithm factor
    const double coulomb_log_factor_;
    
    const double twoPi = 2. * 3.14159265358979323846;
    const double coeff1_, coeff2_, coeff3_, coeff4_;
    
};


#endif
