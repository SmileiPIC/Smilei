#include "CollisionalFusionDD.h"

#include "Collisions.h"
#include "Species.h"
#include "Patch.h"

#include <cmath>

using namespace std;

// Coefficients used for energy interpolation
// The list of energies is in logarithmic scale,
//  with Emin=1 keV, Emax=631 MeV and npoints=50.
const int    CollisionalFusionDD::npoints = 50;
const double CollisionalFusionDD::npointsm1 = ( double )( npoints-1 );
const double CollisionalFusionDD::a1 = log(511./(2.*2.013553)); // = ln(me*c^2 / Emin / n_nucleons)
const double CollisionalFusionDD::a2 = 3.669039; // = (npoints-1) / ln( Emax/Emin )
const double CollisionalFusionDD::a3 = log(511./0.7/(2.*2.013553));; // = ln(me*c^2 / Eref / n_nucleons)
// Log of cross-section in units of 4 pi re^2
const double CollisionalFusionDD::DB_log_crossSection[50] = {
    -27.307, -23.595, -20.383, -17.607, -15.216, -13.167, -11.418, -9.930, -8.666, -7.593,
    -6.682, -5.911, -5.260, -4.710, -4.248, -3.858, -3.530, -3.252, -3.015, -2.814, -2.645,
    -2.507, -2.398, -2.321, -2.273, -2.252, -2.255, -2.276, -2.310, -2.352, -2.402, -2.464,
    -2.547, -2.664, -2.831, -3.064, -3.377, -3.773, -4.249, -4.794, -5.390, -6.021, -6.677,
    -7.357, -8.072, -8.835, -9.648, -10.498, -11.387, -12.459
};

// Constructor
CollisionalFusionDD::CollisionalFusionDD(
    Params *params,
    vector<Particles*> product_particles,
    vector<unsigned int> product_species,
    double rate_multiplier
)
: CollisionalNuclearReaction(params, &product_particles, &product_species, rate_multiplier)
{
}

// Cloning constructor
CollisionalFusionDD::CollisionalFusionDD( CollisionalNuclearReaction *NR )
: CollisionalNuclearReaction( NR )
{
}


double CollisionalFusionDD::crossSection( double log_ekin )
{
    // Interpolate the total cross-section at some value of ekin = m1(g1-1) + m2(g2-1)
    double x = a2*( a1 + log_ekin );
    double cs;
    // if energy below Emin, approximate to 0.
    if( x < 0. ) {
        cs = 0.;
    }
    // if energy within table range, interpolate
    else if( x < npointsm1 ) {
        int i = int( x );
        double a = x - ( double )i;
        cs = exp(
            ( DB_log_crossSection[i+1]-DB_log_crossSection[i] )*a + DB_log_crossSection[i]
        );
    }
    // if energy above table range, extrapolate
    else {
        double a = x - npointsm1;
        cs = exp(
            ( DB_log_crossSection[npoints-1]-DB_log_crossSection[npoints-2] )*a + DB_log_crossSection[npoints-1]
        );
    }
    return cs;
}

void CollisionalFusionDD::makeProducts(
    double U, double ekin, double log_ekin, double q, 
    Particles *&p3, Particles *&p4,
    double &p3_COM, double &p4_COM,
    double &q3, double &q4,
    double &sinX, double &cosX
) {
    U = 2*U - 1.;
    double U1 = abs( U );
    bool up = U > 0.;
    
    // Sample the products angle from empirical fits
    double lnE = a3 + log_ekin;
    double alpha = lnE < 0. ? 1. : exp(-0.024*lnE*lnE);
    double one_m_cosX = alpha*U1 / sqrt( (1.-U1) + alpha*alpha*U1 );
    cosX = 1. - one_m_cosX;
    sinX = sqrt( one_m_cosX * (1.+cosX) );
    
    // Calculate the resulting momenta from energy / momentum conservation
    const double Q = 6.397; // Qvalue
    const double m_n = 1838.7;
    const double m_He = 5497.9;
    double p_COM = sqrt( (ekin+Q) * (ekin+Q+2.*m_n) * (ekin+Q+2.*m_He) * (ekin+Q+2.*m_n+2.*m_He) )
        / ( ( ekin+Q+m_n+m_He ) * (2.*m_He) );
    
    // Set particle properties
    if( up ) {
        p3_COM = p_COM;
        q3 = q;
        p3 = product_particles_[0]; // helium3
    } else {
        p4_COM = p_COM;
        q4 = q;
        p4 = product_particles_[0]; // helium3
    }
}
