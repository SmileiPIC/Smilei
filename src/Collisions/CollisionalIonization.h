
#ifndef COLLISIONALIONIZATION_H
#define COLLISIONALIONIZATION_H

#include <vector>

#include "Tools.h"
#include "Species.h"
#include "Params.h"
#include "BinaryProcessData.h"

class Patch;

class CollisionalIonization
{

public:
    //! Constructor
    CollisionalIonization( int, Params*, int, Particles* );
    //! Cloning Constructor
    CollisionalIonization( CollisionalIonization * );
    //! destructor
    ~CollisionalIonization() {};
    
    SMILEI_ACCELERATOR_DECLARE_ROUTINE
    void apply( Random *random, BinaryProcessData &D, uint32_t n )
    {
    // Not ported to GPU yet
    #ifndef SMILEI_ACCELERATOR_GPU
    
        for( uint32_t i = 0; i<n; i++ ) {
            // Random numbers
            double U1 = random->uniform();
            double U2 = random->uniform();
            // Find which particle in the pair is electron or ion
            uint32_t which_e = (uint32_t) !D.electronFirst;
            uint32_t which_i = (uint32_t) D.electronFirst;
            Particles *& pe = D.p[which_e][i];
            Particles *& pi = D.p[which_i][i];
            uint32_t & ie = D.i[which_e][i];
            uint32_t & ii = D.i[which_i][i];
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
    SMILEI_ACCELERATOR_DECLARE_ROUTINE_END
    
    void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool intra, std::vector<unsigned int> sg1, std::vector<unsigned int> sg2, int itime );
    std::string name() {
        std:: ostringstream t;
        t << "Collisional ionization with atomic number "<<atomic_number<<" towards species #"<<ionization_electrons_;
        return t.str();
    };
    
    //! Initializes the arrays in the database and returns the index of these arrays in the DB
    unsigned int createDatabase( double );
    
    //! Coefficients used for interpolating the energy over a given initial list
    static const double a1, a2, npointsm1;
    static const int npoints;
    
    //! Local table of integrated cross-section
    std::vector<std::vector<double> > *crossSection;
    //! Local table of average secondary electron energy
    std::vector<std::vector<double> > *transferredEnergy;
    //! Local table of average incident electron energy lost
    std::vector<std::vector<double> > *lostEnergy;
    
    //! New electrons temporary species
    Particles new_electrons;
    
    //! Index of the atomic number in the databases
    unsigned int dataBaseIndex;
    
private:
    
    //! Atomic number
    int atomic_number;
    
    //! Species where new electrons are sent
    int ionization_electrons_;
    
    //! Global table of atomic numbers
    static std::vector<int> DB_Z;
    //! Global table of integrated cross-section
    static std::vector<std::vector<std::vector<double> > > DB_crossSection;
    //! Global table of average secondary electron energy
    static std::vector<std::vector<std::vector<double> > > DB_transferredEnergy;
    //! Global table of average incident electron energy lost
    static std::vector<std::vector<std::vector<double> > > DB_lostEnergy;
    
    //! True if first group of species is the electron
    // bool electronFirst;
    
    //! Current ionization rate array (one cell per number of ionization events)
    std::vector<double> rate;
    std::vector<double> irate;
    //! Current ionization probability array (one cell per number of ionization events)
    std::vector<double> prob;
    
};

#endif
