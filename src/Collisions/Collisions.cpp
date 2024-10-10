
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
    coulomb_log_factor_( coulomb_log_factor )
{
    coeff1_ = 4.046650232e-21*params.reference_angular_frequency_SI; // h*omega/(2*me*c^2)
    coeff2_ = 2.817940327e-15*params.reference_angular_frequency_SI/299792458.; // re omega / c
    coeff3_ = coeff2_ * coulomb_log_factor_;
    coeff4_ = pow( 3.*coeff2_, -1./3. );
}


// Cloning Constructor
Collisions::Collisions( Collisions *coll )
{
    coulomb_log_        = coll->coulomb_log_       ;
    coulomb_log_factor_ = coll->coulomb_log_factor_;
    coeff1_             = coll->coeff1_            ;
    coeff2_             = coll->coeff2_            ;
    coeff3_             = coll->coeff3_            ;
    coeff4_             = coll->coeff4_            ;
}


Collisions::~Collisions()
{
}

void Collisions::prepare()
{
    npairs_tot_  = 0.;
    smean_       = 0.;
    logLmean_    = 0.;
}

void Collisions::apply( Random *random, BinaryProcessData &D )
{
    
    vector<double> s( D.n );
    
    if( coulomb_log_ <= 0. ) { // if auto coulomb log
    
        for( size_t i = 0; i<D.n; i++ ) {
            double inv_p_COM = 1./( D.m[0][i] * D.p_COM[i] );
            double gamma0_p = D.gamma0[i] * inv_p_COM;
            double m_p_gamma_COM = D.m[1][i] / ( D.m[0][i] * D.p_COM[i] * D.gamma_tot_COM[i] );
            
            // Calculate coulomb log if necessary
            double logL;
            double qqqqlogL; // q1^2 q2^2 logL
            double qq = ( double ) D.q[0][i] * ( double ) D.q[1][i];
            double qqqq = qq * qq;
            // Note : 0.00232282 is coeff2 / coeff1
            double bmin = coeff1_ * std::max( inv_p_COM, std::abs( 0.00232282*qq*m_p_gamma_COM*gamma0_p ) ); // min impact parameter
            // Debye screening
            double lnLD = D.debye > 7.839*bmin ? log( D.debye / bmin ) : 2.;
            // For e-i collisions, consider the bound-electron screening
            if( D.lTF[i] > 0. ) {
                if( D.debye > D.lTF[i] ) {
                    // Bound-electron (Thomas-Fermi) screening
                    double lnLTF = D.lTF[i] > 7.839*bmin ? log( D.lTF[i] / bmin ) : 2.;
                    // Total screening
                    qqqqlogL = qqqq * lnLD + ( D.Z1Z2[i] * D.Z1Z2[i] - qqqq ) * lnLTF;
                    logL = qqqqlogL / ( D.Z1Z2[i] * D.Z1Z2[i] );
                } else {
                    logL = lnLD;
                    qqqqlogL = D.Z1Z2[i] * D.Z1Z2[i] * lnLD;
                }
            } else {
                logL = lnLD;
                qqqqlogL = qqqq * lnLD;
            }
            
            // Calculate the collision parameter s12 (similar to number of real collisions)
            s[i] = coeff3_ * qqqqlogL * gamma0_p * gamma0_p * m_p_gamma_COM / ( D.gamma[0][i] * D.gamma[1][i] );
            logLmean_ += logL;
            smean_    += s[i];
        }
        
    } else {
        
        for( size_t i = 0; i<D.n; i++ ) {
            double inv_p_COM = 1./( D.m[0][i] * D.p_COM[i] );
            double gamma0_p = D.gamma0[i] * inv_p_COM;
            double m_p_gamma_COM = D.m[1][i] / ( D.m[0][i] * D.p_COM[i] * D.gamma_tot_COM[i] );
            
            double qq = ( double ) D.q[0][i] * ( double ) D.q[1][i];
            
            // Calculate the collision parameter s12 (similar to number of real collisions)
            s[i] = coeff3_ * qq * qq * coulomb_log_ * gamma0_p * gamma0_p * m_p_gamma_COM / ( D.gamma[0][i] * D.gamma[1][i] );
            logLmean_ += coulomb_log_;
            smean_    += s[i];
        }
        
    }
    
    for( size_t i = 0; i<D.n; i++ ) {
        // Low-temperature correction
        double smax = coeff4_ * ( D.m[0][i] + D.m[1][i] ) * D.vrel[i] / std::max( D.m[0][i] * D.n123, D.m[1][i] * D.n223 );
        if( s[i] > smax ) {
            s[i] = smax;
        }
        
        s[i] *= D.dt_correction[i];
        
        // Pick the deflection angles in the center-of-mass frame.
        // Instead of Nanbu http://dx.doi.org/10.1103/PhysRevE.55.4642
        // and Perez http://dx.doi.org/10.1063/1.4742167
        // we made a new fit (faster and more accurate)
        double cosX, sinX;
        double U1 = random->uniform();
        if( s[i] < 4. ) {
            double s2 = s[i]*s[i];
            double alpha = 0.37*s[i] - 0.005*s2 - 0.0064*s2*s[i];
            double sin2X2 = alpha * U1 / sqrt( (1.-U1) + alpha*alpha*U1 );
            cosX = 1. - 2.*sin2X2;
            sinX = 2.*sqrt( sin2X2 *(1.-sin2X2) );
        } else {
            cosX = 2.*U1 - 1.;
            sinX = sqrt( 1. - cosX*cosX );
        }
        
        // Calculate combination of angles
        double phi = random->uniform_2pi();
        double sinXcosPhi = sinX*cos( phi );
        double sinXsinPhi = sinX*sin( phi );
        
        // Apply the deflection
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
        
        // Go back to the lab frame and store the results in the particle array
        double pp = ( D.px_tot[i] * newpx_COM + D.py_tot[i] * newpy_COM + D.pz_tot[i] * newpz_COM ) / ( D.gamma_tot[i] + D.gamma_tot_COM[i] );
        double f = ( D.gamma_COM0[i] + pp ) / D.gamma_tot_COM[i];
        double U2 = random->uniform();
        if( U2 * D.W[0][i] < D.W[1][i] ) { // deflect particle 1 only with some probability
            D.px[0][i] = newpx_COM + f * D.px_tot[i];
            D.py[0][i] = newpy_COM + f * D.py_tot[i];
            D.pz[0][i] = newpz_COM + f * D.pz_tot[i];
        }
        if( U2 * D.W[1][i] < D.W[0][i] ) { // deflect particle 2 only with some probability
            double m12 = D.m[0][i] / D.m[1][i];
            D.px[1][i] = ( -newpx_COM + ( 1 - f ) * D.px_tot[i] ) * m12;
            D.py[1][i] = ( -newpy_COM + ( 1 - f ) * D.py_tot[i] ) * m12;
            D.pz[1][i] = ( -newpz_COM + ( 1 - f ) * D.pz_tot[i] ) * m12;
        }
        
    }
    
    npairs_tot_ += D.n;
}

void Collisions::finish( Params &, Patch *, std::vector<Diagnostic *> &, bool, std::vector<unsigned int>, std::vector<unsigned int>, int )
{
    if( npairs_tot_>0. ) {
        smean_    /= npairs_tot_;
        logLmean_ /= npairs_tot_;
    }
}

