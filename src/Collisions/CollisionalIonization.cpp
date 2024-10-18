#include "CollisionalIonization.h"

#include "Collisions.h"
#include "Species.h"
#include "Patch.h"
#include "IonizationTables.h"

#include <cmath>


using namespace std;

// Coefficients used for energy interpolation
// The list of energies has to be in logarithmic scale,
//  with Emin=1eV, Emax=10MeV and npoints=100.
const int    CollisionalIonization::npoints = 100;
const double CollisionalIonization::npointsm1 = ( double )( npoints-1 );
const double CollisionalIonization::a1 = 510998.9 ; // = me*c^2/Emin
const double CollisionalIonization::a2 = 6.142165 ; // = (npoints-1) / ln( Emax/Emin )

// Constructor
CollisionalIonization::CollisionalIonization( int Z, Params *params, int ionization_electrons, Particles* particles )
{
    atomic_number = Z;
    rate .resize( Z );
    irate.resize( Z );
    prob .resize( Z );
    ionization_electrons_ = ionization_electrons;
    if( params ) {
        new_electrons.initialize( 0, *particles );
    }
    if( Z>0 ) {
        dataBaseIndex = createDatabase( params->reference_angular_frequency_SI );
    }
}

// Cloning Constructor
CollisionalIonization::CollisionalIonization( CollisionalIonization *CI )
{
    atomic_number = CI->atomic_number;
    rate .resize( atomic_number );
    irate.resize( atomic_number );
    prob .resize( atomic_number );
    ionization_electrons_ = CI->ionization_electrons_;
    new_electrons.initialize( 0, CI->new_electrons );
    
    dataBaseIndex = CI->dataBaseIndex;
}

// Static members
vector<int> CollisionalIonization::DB_Z;
vector<vector<vector<double> > > CollisionalIonization::DB_crossSection;
vector<vector<vector<double> > > CollisionalIonization::DB_transferredEnergy;
vector<vector<vector<double> > > CollisionalIonization::DB_lostEnergy;

// Initializes the databases (by patch master only)
unsigned int CollisionalIonization::createDatabase( double reference_angular_frequency_SI )
{
    // Leave if the database already exists with same atomic number
    for( unsigned int i=0; i<DB_Z.size(); i++ ) {
        if( atomic_number == DB_Z[i] ) {
            return i;
        }
    }
    
    // Otherwise, create the arrays:
    // For each ionization state, calculate the tables of integrated cross-sections
    // Pérez et al., Phys. Plasmas 19, 083104 (2012)
    vector<vector<double> > cs; // cross section
    vector<vector<double> > te; // transferred energy
    vector<vector<double> > le; // lost energy
    cs.resize( atomic_number );
    te.resize( atomic_number );
    le.resize( atomic_number );
    double e, ep, bp, up, ep2, betae2, betab2, betau2, s0, A1, A2, A3, sk, wk, ek;
    int N; // occupation number
    double normalization = 2.81794e-15 * reference_angular_frequency_SI / ( 2.*299792458. ); // r_e omega / 2c
    for( int Zstar=0; Zstar<atomic_number; Zstar++ ) { // For each ionization state
        cs[Zstar].resize( npoints, 0. );
        te[Zstar].resize( npoints, 0. );
        le[Zstar].resize( npoints, 0. );
        for( int i=0; i<npoints; i++ ) { // For each incident electron energy
            ep = exp( double( i )/a2 ) / a1; // = incident electron energy
            N = 1;
            for( int k=0; k<atomic_number-Zstar; k++ ) { // For each orbital
                bp = IonizationTables::binding_energy( atomic_number, Zstar, k );
                // If next orbital is on same level, then continue directly to next
                if( k<atomic_number-Zstar-1 && bp == IonizationTables::binding_energy( atomic_number, Zstar, k+1 ) ) {
                    N++;
                    continue;
                }
                // If electron energy below the ionization energy, then skip to next level
                e = ep/bp;
                if( e>1. ) {
                    up = bp; // we assume up=bp because we don't have exact tables
                    betae2 = 1. - 1./( ( 1.+ep )*( 1.+ep ) );
                    betab2 = 1. - 1./( ( 1.+bp )*( 1.+bp ) );
                    betau2 = 1. - 1./( ( 1.+up )*( 1.+up ) );
                    s0 = normalization * N /( bp * ( betae2 + betab2 + betau2 ) );
                    ep2 = 1./( 1.+ep*0.5 );
                    ep2 *= ep2;
                    A1 = ( 1.+2.*ep )/( 1.+e )*ep2;
                    A2 = ( e-1. )*bp*bp*0.5*ep2;
                    A3 = log( betae2/( 1.-betae2 ) ) - betae2 - log( 2.*bp );
                    sk = s0*( 0.5*A3*( 1.-1./( e*e ) ) + 1. - 1./e + A2 - A1*log( e ) );
                    wk = s0 * ( 0.5*A3*( e-1. )*( e-1. )/e/( e+1. )  + 2.*log( 0.5*( e+1. ) ) - log( e )
                                + 0.25*A2*( e-1. ) - A1*( e*log( e )-( e+1. )*log( 0.5*( e+1. ) ) ) );
                    ek = wk + sk;
                    // Sum these data to the total ones
                    cs[Zstar][i] += sk;
                    te[Zstar][i] += wk * bp;
                    le[Zstar][i] += ek * bp;
                }
                // Reset occupation number for next level
                N = 1;
            }
            // The transferred and lost energies are averages over the orbitals
            if( cs[Zstar][i]>0. ) {
                te[Zstar][i] /= cs[Zstar][i];
                le[Zstar][i] /= cs[Zstar][i];
            }
        }
    }
    
    // Add the new arrays to the static database
    DB_Z                .push_back( atomic_number );
    DB_crossSection     .push_back( cs );
    DB_transferredEnergy.push_back( te );
    DB_lostEnergy       .push_back( le );
    
    return DB_Z.size()-1;
}

// Method to apply the ionization
void CollisionalIonization::apply( Random *random, BinaryProcessData &D, size_t n )
{

// Not ported to GPU yet
#ifndef SMILEI_ACCELERATOR_GPU

    for( size_t i = 0; i<n; i++ ) {
        // Random numbers
        double U1 = random->uniform();
        double U2 = random->uniform();
        // Find which particle in the pair is electron or ion
        size_t which_e = (size_t) !D.electronFirst;
        size_t which_i = (size_t) D.electronFirst;
        Particles *& pe = D.p[which_e][i];
        Particles *& pi = D.p[which_i][i];
        size_t & ie = D.i[which_e][i];
        size_t & ii = D.i[which_i][i];
        double & We = D.W[which_e][i];
        double & Wi = D.W[which_i][i];
        double & gammae = D.gamma[which_e][i];
        double & gammai = D.gamma[which_i][i];
        short & qe = D.q[which_e][i];
        short & qi = D.q[which_i][i];
        
        int Zstar = qi;
        if( Zstar>=atomic_number ) {
            continue;    // if already fully ionized, do nothing
        }
        
        // Calculate coefficient (1-ve.vi)*ve' where ve' is in ion frame
        double K = D.dt_correction[i] * sqrt( D.gamma0[i]*D.gamma0[i]-1. )/gammai;
        
        // Fetch cross sections
        const vector<vector<double> > & crossSection = DB_crossSection[dataBaseIndex];
        const vector<vector<double> > & transferredEnergy = DB_transferredEnergy[dataBaseIndex];
        const vector<vector<double> > & lostEnergy = DB_lostEnergy[dataBaseIndex];
        
        // Loop for multiple ionization
        // k+1 is the number of ionizations
        const int kmax = atomic_number-Zstar-1;
        double cs, w, e, cum_prob = 0;
        for( int k = 0; k <= kmax;  k++ ) {
            // Calculate the location x (~log of energy) in the databases
            const double x = a2*log( a1*( D.gamma0[i]-1. ) );
            
            // Interpolate the databases at location x
            if( x < 0. ) {
                break;    // if energy below Emin, do nothing
            }
            if( x < npointsm1 ) { // if energy within table range, interpolate
                const int j = int( x );
                const double a = x - ( double )j;
                cs = ( crossSection[Zstar][j+1]-crossSection[Zstar][j] )*a + crossSection[Zstar][j];
                w  = ( transferredEnergy[Zstar][j+1]-transferredEnergy[Zstar][j] )*a + transferredEnergy[Zstar][j];
                e  = ( lostEnergy[Zstar][j+1]-lostEnergy[Zstar][j] )*a + lostEnergy[Zstar][j];
            } else { // if energy above table range, extrapolate
                const double a = x - npointsm1;
                cs = ( crossSection[Zstar][npoints-1]-crossSection[Zstar][npoints-2] )*a + crossSection[Zstar][npoints-1];
                w  = transferredEnergy[Zstar][npoints-1];
                e  = lostEnergy[Zstar][npoints-1];
            }
            if( e > D.gamma0[i]-1. ) {
                break;
            }
            
            rate[k] = K*cs/gammae  ; // k-th ionization rate
            irate[k] = 1./rate[k]  ; // k-th ionization inverse rate
            prob[k] = exp( -rate[k] ); // k-th ionization probability
            
            // Calculate the cumulative probability for k-th ionization (Nuter et al, 2011)
            if( k==0 ) {
                cum_prob = prob[k];
            } else {
                for( int p=0; p<k; p++ ) {
                    double cp = 1. - rate[k]*irate[p];
                    for( int j=0  ; j<p; j++ ) {
                        cp *= 1.-rate[p]*irate[j];
                    }
                    for( int j=p+1; j<k; j++ ) {
                        cp *= 1.-rate[p]*irate[j];
                    }
                    cum_prob += ( prob[k]-prob[p] )/cp;
                }
            }
            
            // If no more ionization, leave
            if( U1 < cum_prob ) {
                break;
            }
            
            // Otherwise, we do the ionization
            const double p2 = D.gamma0[i]*D.gamma0[i] - 1.;
            // Ionize the atom and create electron
            if( U2 < We/Wi ) {
                qi++; // increase ion charge
                // Calculate the new electron momentum with correction for moving back to lab frame
                double pr1 = sqrt( w*( w+2. )/p2 );
                double pr2 = w+1. - pr1*D.gamma0[i];
                double newpx = D.px[which_e][i] * pr1 + D.px[which_i][i] * pr2;
                double newpy = D.py[which_e][i] * pr1 + D.py[which_i][i] * pr2;
                double newpz = D.pz[which_e][i] * pr1 + D.pz[which_i][i] * pr2;
                // New electron has ion position
                new_electrons.makeParticleAt( *pi, ii, Wi, qe, newpx, newpy, newpz );
                // If quantum parameter exists for new electron, then calculate it
                if( new_electrons.has_quantum_parameter ) {
                    new_electrons.Chi.back() = pe->chi( ie ) * (w+1.) / gammae;
                }
            }
            // Lose incident electron energy
            if( U2 < Wi/We ) {
                // Calculate the modified electron momentum
                double pr = sqrt( ( pow( D.gamma0[i]-e, 2 )-1. )/p2 );
                D.px[which_e][i] *= pr;
                D.py[which_e][i] *= pr;
                D.pz[which_e][i] *= pr;
                gammae *= pr;
                K *= pr;
                // Correction for moving back to the lab frame
                pr = ( 1.-pr )*D.gamma0[i] - e;
                D.px[which_e][i] += pr * D.px[which_i][i];
                D.py[which_e][i] += pr * D.py[which_i][i];
                D.pz[which_e][i] += pr * D.pz[which_i][i];
                gammae += pr * gammai;
                // Decrease gamma for next ionization
                D.gamma0[i] -= e;
            }
            
            Zstar++;
            
        }
    }
#endif
}


// Finish the ionization (moves new electrons in place)
void CollisionalIonization::finish( Params &params, Patch *patch, std::vector<Diagnostic *> &localDiags, bool, std::vector<unsigned int>, std::vector<unsigned int>, int itime )
{
    patch->vecSpecies[ionization_electrons_]->importParticles( params, patch, new_electrons, localDiags, ( itime + 0.5 ) * params.timestep );
}
