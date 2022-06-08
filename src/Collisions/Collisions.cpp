
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
    double qqm  = D.p1->charge( D.i1 ) * D.p2->charge( D.i2 ) / D.m1;
    double qqm2 = qqm * qqm;
    
    // Calculate coulomb log if necessary
    double logL = coulomb_log_;
    if( logL <= 0. ) { // if auto-calculation requested
        // Note : 0.00232282 is coeff2 / coeff1
        double bmin = coeff1_ * std::max( 1./(D.m1*D.p_COM), std::abs( 0.00232282*qqm*D.term3*D.term5 ) ); // min impact parameter
        logL = 0.5*log( 1. + D.debye2/( bmin*bmin ) );
        if( logL < 2. ) {
            logL = 2.;
        }
    }
    
    // Calculate the collision parameter s12 (similar to number of real collisions)
    double s = coeff3_ * logL * qqm2 * D.term3 * D.p_COM * D.term5*D.term5 / ( D.gamma1*D.gamma2 );
    
    // Low-temperature correction
    double smax = coeff4_ * ( D.m12+1. ) * D.vrel / std::max( D.m12*D.n123, D.n223 );
    if( s>smax ) {
        s = smax;
    }
    
    s *= D.dt_correction;
    
    // Pick the deflection angles in the center-of-mass frame.
    // Instead of Nanbu http://dx.doi.org/10.1103/PhysRevE.55.4642
    // and Perez http://dx.doi.org/10.1063/1.4742167
    // we made a new fit (faster and more accurate)
    double cosX, sinX;
    double U1 = random->uniform();
    if( s < 4. ) {
        double s2 = s*s;
        double alpha = 0.37*s - 0.005*s2 - 0.0064*s2*s;
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
    double p_perp = sqrt( D.px_COM*D.px_COM + D.py_COM*D.py_COM );
    double newpx_COM, newpy_COM, newpz_COM;
    if( p_perp > 1.e-10*D.p_COM ) { // make sure p_perp is not too small
        double inv_p_perp = 1./p_perp;
        newpx_COM = ( D.px_COM * D.pz_COM * sinXcosPhi - D.py_COM * D.p_COM * sinXsinPhi ) * inv_p_perp + D.px_COM * cosX;
        newpy_COM = ( D.py_COM * D.pz_COM * sinXcosPhi + D.px_COM * D.p_COM * sinXsinPhi ) * inv_p_perp + D.py_COM * cosX;
        newpz_COM = -p_perp * sinXcosPhi + D.pz_COM * cosX;
    } else { // if p_perp is too small, we use the limit px->0, py=0
        newpx_COM = D.p_COM * sinXcosPhi;
        newpy_COM = D.p_COM * sinXsinPhi;
        newpz_COM = D.p_COM * cosX;
    }
    
    // Go back to the lab frame and store the results in the particle array
    double vcp = D.COM_vx * newpx_COM + D.COM_vy * newpy_COM + D.COM_vz * newpz_COM;
    double U2 = random->uniform();
    if( U2 * D.p1->weight( D.i1 ) < D.p2->weight( D.i2 ) ) { // deflect particle 1 only with some probability
        double term6 = D.term1*vcp + D.gamma1_COM * D.COM_gamma;
        D.p1->momentum( 0, D.i1 ) = newpx_COM + D.COM_vx * term6;
        D.p1->momentum( 1, D.i1 ) = newpy_COM + D.COM_vy * term6;
        D.p1->momentum( 2, D.i1 ) = newpz_COM + D.COM_vz * term6;
    }
    if( U2 * D.p2->weight( D.i2 ) < D.p1->weight( D.i1 ) ) { // deflect particle 2 only with some probability
        double term6 = -D.m12 * D.term1*vcp + D.gamma2_COM * D.COM_gamma;
        D.p2->momentum( 0, D.i2 ) = -D.m12 * newpx_COM + D.COM_vx * term6;
        D.p2->momentum( 1, D.i2 ) = -D.m12 * newpy_COM + D.COM_vy * term6;
        D.p2->momentum( 2, D.i2 ) = -D.m12 * newpz_COM + D.COM_vz * term6;
    }
    
    npairs_tot_ ++;
    smean_    += s;
    logLmean_ += logL;
}

void Collisions::finish( Params &, Patch *, std::vector<Diagnostic *> &, bool intra, std::vector<unsigned int> sg1, std::vector<unsigned int> sg2, int itime )
{
    if( npairs_tot_>0. ) {
        smean_    /= npairs_tot_;
        logLmean_ /= npairs_tot_;
    }
}

