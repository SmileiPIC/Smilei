
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <fstream>

#include "BinaryProcesses.h"
#include "BinaryProcess.h"
#include "Collisions.h"
#include "CollisionalNuclearReaction.h"
#include "CollisionalFusionDD.h"
#include "CollisionalIonization.h"
#include "H5.h"
#include "Patch.h"
#include "VectorPatch.h"

using namespace std;


// Constructor
BinaryProcesses::BinaryProcesses(
    Params &params,
    vector<unsigned int> species_group1,
    vector<unsigned int> species_group2,
    bool intra,
    vector<BinaryProcess*> processes,
    int every,
    int debug_every,
    double time_frozen,
    string filename
) :
    processes_( processes ),
    species_group1_( species_group1 ),
    species_group2_( species_group2 ),
    intra_( intra ),
    every_( every ),
    debug_every_( debug_every ),
    filename_( filename )
{
    timesteps_frozen_ = time_frozen / params.timestep;
    // Open the HDF5 file
    debug_file_ = NULL;
    if( debug_every > 0 ) {
        MPI_Comm comm = MPI_COMM_WORLD;
        debug_file_ = new H5Write( filename_, &comm );
    }
}


// Cloning Constructor
BinaryProcesses::BinaryProcesses( BinaryProcesses *BPs )
{
    species_group1_     = BPs->species_group1_    ;
    species_group2_     = BPs->species_group2_    ;
    intra_              = BPs->intra_             ;
    every_              = BPs->every_             ;
    debug_every_        = BPs->debug_every_       ;
    timesteps_frozen_   = BPs->timesteps_frozen_  ;
    filename_           = BPs->filename_          ;
    debug_file_         = BPs->debug_file_        ;
    
    processes_.clear();
    for( unsigned int i=0; i<BPs->processes_.size(); i++ ) {
        if( Collisions * coll = dynamic_cast<Collisions*>( BPs->processes_[i] ) ) {
            processes_.push_back( new Collisions( coll ) );
        } else if( CollisionalIonization * CI = dynamic_cast<CollisionalIonization*>( BPs->processes_[i] ) ) {
            processes_.push_back( new CollisionalIonization( CI ) );
        } else if( CollisionalFusionDD * DD = dynamic_cast<CollisionalFusionDD*>( BPs->processes_[i] ) ) {
            processes_.push_back( new CollisionalFusionDD( DD ) );
        } else {
            ERROR( "Undefined binary process" );
        }
    }
    
}


BinaryProcesses::~BinaryProcesses()
{
    for( unsigned int i=0; i<processes_.size(); i++ ) {
        delete processes_[i];
    }
    processes_.clear();
}

// Declare other static variables here
bool BinaryProcesses::debye_length_required_;


// Calculates the debye length squared in each bin
// The formula for the inverse debye length squared is sumOverSpecies(density*charge^2/temperature)
void BinaryProcesses::calculate_debye_length( Params &params, Patch *patch )
{
    Species    *s;
    Particles *p;
    double coeff = 299792458./( 3.*params.reference_angular_frequency_SI*2.8179403267e-15 ); // c / (3 omega re)
    
    unsigned int nspec = patch->vecSpecies.size(); // number of species
    if( nspec==0 ) {
        return;
    }
    unsigned int nbin = patch->vecSpecies[0]->particles->first_index.size();
    
    patch->debye_length_squared.resize( nbin );
    double mean_debye_length = 0.;
    for( unsigned int ibin = 0 ; ibin < nbin ; ibin++ ) {
        double density_max = 0.;
        double inv_D2 = 0.;
        double inv_cell_volume = 0.;
        
        for( unsigned int ispec=0 ; ispec<nspec ; ispec++ ) { // loop all species
            s  = patch->vecSpecies[ispec];
            p  = s->particles;
            
            // Skip when no particles
            if( p->last_index[ibin] <= p->first_index[ibin] ) continue;
            
            if( inv_cell_volume == 0. ) {
                inv_cell_volume = 1. / patch->getPrimalCellVolume( p, p->first_index[ibin], params );
            }
            
            // Calculation of particles density, mean charge, and temperature
            // Density is the sum of weights
            // Temperature definition is the average <v*p> divided by 3
            double density     = 0.;
            double charge      = 0.;
            double temperature = 0.;
            // loop particles to calculate average quantities
            for( unsigned int iPart=p->first_index[ibin]; iPart<( unsigned int )p->last_index[ibin] ; iPart++ ) {
                double p2 = p->momentum( 0, iPart ) * p->momentum( 0, iPart )
                     +p->momentum( 1, iPart ) * p->momentum( 1, iPart )
                     +p->momentum( 2, iPart ) * p->momentum( 2, iPart );
                density     += p->weight( iPart );
                charge      += p->weight( iPart ) * p->charge( iPart );
                temperature += p->weight( iPart ) * p2/sqrt( 1.+p2 );
            }
            if( density > 0. ) {
                charge /= density; // average charge
                temperature *= s->mass_ / ( 3.*density ); // Te in units of me*c^2
                density *= inv_cell_volume; // density in units of critical density
                // compute inverse debye length squared
                if( temperature == 0. ) {
                    inv_D2 += 1e100; // infinite
                } else {
                    inv_D2 += density*charge*charge/temperature;
                }
                // compute maximum density of species
                if( density>density_max ) {
                    density_max = density;
                }
            }
        }
        
        // if there were particles,
        if( inv_D2 > 0. ) {
            // compute debye length squared in code units
            patch->debye_length_squared[ibin] = 1./inv_D2;
            // apply lower limit to the debye length (minimum interatomic distance)
            double rmin2 = pow( coeff*density_max, -2./3. );
            if( patch->debye_length_squared[ibin] < rmin2 ) {
                patch->debye_length_squared[ibin] = rmin2;
            }
        } else {
            patch->debye_length_squared[ibin] = 0.;
        }
        
        mean_debye_length += sqrt( patch->debye_length_squared[ibin] );
    }
    
    mean_debye_length /= nbin;
    DEBUG( "Mean Debye length in code length units = " << scientific << setprecision( 3 ) << mean_debye_length );
    mean_debye_length *= 299792458./params.reference_angular_frequency_SI; // switch to SI
    DEBUG( "Mean Debye length in meters = " << scientific << setprecision( 3 ) << mean_debye_length );
}


// Make the pairing and launch the processes
void BinaryProcesses::apply( Params &params, Patch *patch, int itime, vector<Diagnostic *> &localDiags )
{
    if( itime < timesteps_frozen_ || itime % every_ != 0 ) {
        return;
    }
    
    vector<unsigned int> *sg1, *sg2, index1, index2;
    unsigned int nspec1, nspec2; // numbers of species in each group
    unsigned int npart1, npart2; // numbers of macro-particles in each group
    unsigned int npairs; // number of pairs of macro-particles
    vector<unsigned int> np1, np2; // numbers of macro-particles in each species, in each group
    unsigned int ispec1, ispec2, N2max;
    Species   *s1, *s2;
    Particles *p1=NULL, *p2;
    
    sg1 = &species_group1_;
    sg2 = &species_group2_;
    
    for( unsigned int i=0; i<processes_.size(); i++ ) {
        processes_[i]->prepare();
    }
    
    // Loop bins of particles (typically, cells, but may also be clusters)
    unsigned int nbin = patch->vecSpecies[0]->particles->first_index.size();
    for( unsigned int ibin = 0 ; ibin < nbin ; ibin++ ) {
        
        BinaryProcessData D;
        
        // get number of particles for all necessary species
        for( unsigned int i=0; i<2; i++ ) { // try twice to ensure group 1 has more macro-particles
            nspec1 = sg1->size();
            nspec2 = sg2->size();
            np1.resize( nspec1 ); // number of particles in each species of group 1
            np2.resize( nspec2 ); // number of particles in each species of group 2
            npart1 = 0;
            npart2 = 0;
            for( ispec1=0 ; ispec1<nspec1 ; ispec1++ ) {
                s1 = patch->vecSpecies[( *sg1 )[ispec1]];
                np1[ispec1] = s1->particles->last_index[ibin] - s1->particles->first_index[ibin];
                npart1 += np1[ispec1];
            }
            for( ispec2=0 ; ispec2<nspec2 ; ispec2++ ) {
                s2 = patch->vecSpecies[( *sg2 )[ispec2]];
                np2[ispec2] = s2->particles->last_index[ibin] - s2->particles->first_index[ibin];
                npart2 += np2[ispec2];
            }
            if( npart2 <= npart1 ) {
                break;    // ok if group1 has more macro-particles
            } else { // otherwise, we exchange groups and try again
                swap( sg1, sg2 );
            }
        }
        // now group1 has more macro-particles than group2
        
        // skip to next bin if no particles
        if( npart1==0 || npart2==0 ) {
            continue;
        }
        
        // Shuffle particles to have random pairs
        //    (It does not really exchange them, it is just a temporary re-indexing)
        index1.resize( npart1 );
        for( unsigned int i=0; i<npart1; i++ ) {
            index1[i] = i;    // first, we make an ordered array
        }
        // shuffle the index array
        for( unsigned int i=npart1; i>1; i-- ) {
            unsigned int p = patch->rand_->integer() % i;
            swap( index1[i-1], index1[p] );
        }
        if( intra_ ) { // In the case of pairing within one species
            if( npart1 < 2 ) {
                continue;
            }
            npairs = ( npart1 + 1 ) / 2; // half as many pairs as macro-particles
            index2.resize( npairs );
            for( unsigned int i=0; i<npairs; i++ ) {
                index2[i] = index1[( i+npairs )%npart1];    // index2 is second half
            }
            index1.resize( npairs ); // index1 is first half
            N2max = npart1 - npairs; // number of not-repeated particles (in group 2 only)
        } else { // In the case of pairing between two species
            npairs = npart1; // as many pairs as macro-particles in group 1 (most numerous)
            index2.resize( npairs );
            for( unsigned int i=0; i<npart1; i++ ) {
                index2[i] = i % npart2;
            }
            N2max = npart2; // number of not-repeated particles (in group 2 only)
        }
        
        // Data for ionization
        D.electronFirst = patch->vecSpecies[( *sg1 )[0]]->atomic_number_==0 ? true : false;
        
        // Calculate the densities
        double n1  = 0.; // density of group 1
        for( ispec1=0 ; ispec1<nspec1 ; ispec1++ ) {
            s1 = patch->vecSpecies[( *sg1 )[ispec1]];
            p1 = s1->particles;
            for( int i = p1->first_index[ibin]; i < p1->last_index[ibin]; i++ ) {
                n1 += p1->weight( i );
            }
        }
        double n2  = 0.; // density of group 2
        for( ispec2=0 ; ispec2<nspec2 ; ispec2++ ) {
            s2 = patch->vecSpecies[( *sg2 )[ispec2]];
            p2 = s2->particles;
            for( int i = p2->first_index[ibin]; i < p2->last_index[ibin]; i++ ) {
                n2 += p2->weight( i );
            }
        }
        
        // Get cell volume
        ispec1 = -1;
        do {
            ispec1++;
            p1 = patch->vecSpecies[( *sg1 )[ispec1]]->particles;
        } while( ispec1<nspec1 && p1->first_index[ibin] == p1->last_index[ibin] );
        double inv_cell_volume = 1./patch->getPrimalCellVolume( p1, p1->first_index[ibin], params );
        
        // Set the debye length
        if( BinaryProcesses::debye_length_required_ ) {
            D.debye2 = patch->debye_length_squared[ibin];
        }
        
        // Pre-calculate some numbers before the big loop
        unsigned int ncorr = intra_ ? 2*npairs-1 : npairs;
        double dt_corr = every_ * params.timestep * ((double)ncorr) * inv_cell_volume;
        double weight_correction_1 = 1. / (double)( (npairs-1) / N2max );
        double weight_correction_2 = 1. / (double)( (npairs-1) / N2max + 1 );
        n1  *= inv_cell_volume;
        n2  *= inv_cell_volume;
        D.n123 = pow( n1, 2./3. );
        D.n223 = pow( n2, 2./3. );
        
        // Now start the real loop on pairs of particles
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        for( unsigned int i = 0; i<npairs; i++ ) {
            
            // find species and index i1 of particle "1"
            D.i1 = index1[i];
            for( ispec1=0 ; D.i1>=np1[ispec1]; ispec1++ ) {
                D.i1 -= np1[ispec1];
            }
            // find species and index i2 of particle "2"
            D.i2 = index2[i];
            for( ispec2=0 ; D.i2>=np2[ispec2]; ispec2++ ) {
                D.i2 -= np2[ispec2];
            }
            
            s1 = patch->vecSpecies[( *sg1 )[ispec1]];
            s2 = patch->vecSpecies[( *sg2 )[ispec2]];
            D.i1 += s1->particles->first_index[ibin];
            D.i2 += s2->particles->first_index[ibin];
            D.p1 = s1->particles;
            D.p2 = s2->particles;
            
            D.m1 = s1->mass_;
            D.m2 = s2->mass_;
            D.m12 = s1->mass_ / s2->mass_;
            
            D.minW = D.p1->weight(D.i1);
            D.maxW = D.p2->weight(D.i2);
            if( D.minW > D.maxW ) {
                swap( D.minW, D.maxW );
            }
            
            // If one weight is zero, then skip. Can happen after nuclear reaction
            if( D.minW <= 0. ) continue;
            
            D.dt_correction = D.maxW * dt_corr;
            if( i % N2max <= (npairs-1) % N2max ) {
                D.dt_correction *= weight_correction_2 ;
            } else {
                D.dt_correction *= weight_correction_1;
            }
            
            // Get momenta and calculate gammas
            
            D.gamma1 = D.p1->LorentzFactor( D.i1 );
            D.gamma2 = D.p2->LorentzFactor( D.i2 );
            double gamma12 = D.m12 * D.gamma1 + D.gamma2;
            double gamma12_inv = 1./gamma12;
            
            // Calculate the center-of-mass (COM) frame
            // Quantities starting with "COM" are those of the COM itself, expressed in the lab frame.
            // They are NOT quantities relative to the COM.
            D.COM_vx = ( D.m12 * ( D.p1->momentum( 0, D.i1 ) ) + D.p2->momentum( 0, D.i2 ) ) * gamma12_inv;
            D.COM_vy = ( D.m12 * ( D.p1->momentum( 1, D.i1 ) ) + D.p2->momentum( 1, D.i2 ) ) * gamma12_inv;
            D.COM_vz = ( D.m12 * ( D.p1->momentum( 2, D.i1 ) ) + D.p2->momentum( 2, D.i2 ) ) * gamma12_inv;
            double COM_vsquare = D.COM_vx*D.COM_vx + D.COM_vy*D.COM_vy + D.COM_vz*D.COM_vz;
            
            // Change the momentum to the COM frame (we work only on particle 1)
            // Quantities ending with "COM" are quantities of the particle expressed in the COM frame.
            if( COM_vsquare < 1e-6 ) {
                D.COM_gamma = 1. +0.5 * COM_vsquare;
                D.term1 = 0.5;
            } else {
                D.COM_gamma = 1./sqrt( 1.-COM_vsquare );
                D.term1 = ( D.COM_gamma - 1. ) / COM_vsquare;
            }
            double vcv1g1  = D.COM_vx*( D.p1->momentum( 0, D.i1 ) ) + D.COM_vy*( D.p1->momentum( 1, D.i1 ) ) + D.COM_vz*( D.p1->momentum( 2, D.i1 ) );
            double vcv2g2  = D.COM_vx*( D.p2->momentum( 0, D.i2 ) ) + D.COM_vy*( D.p2->momentum( 1, D.i2 ) ) + D.COM_vz*( D.p2->momentum( 2, D.i2 ) );
            D.gamma1_COM = ( D.gamma1-vcv1g1 )*D.COM_gamma;
            D.gamma2_COM = ( D.gamma2-vcv2g2 )*D.COM_gamma;
            double term2 = D.term1*vcv1g1 - D.COM_gamma * D.gamma1;
            D.px_COM = D.p1->momentum( 0, D.i1 ) + term2*D.COM_vx;
            D.py_COM = D.p1->momentum( 1, D.i1 ) + term2*D.COM_vy;
            D.pz_COM = D.p1->momentum( 2, D.i1 ) + term2*D.COM_vz;
            double p2_COM = D.px_COM*D.px_COM + D.py_COM*D.py_COM + D.pz_COM*D.pz_COM;
            D.p_COM  = sqrt( p2_COM );
            
            // Calculate some intermediate quantities
            D.term3 = D.COM_gamma * gamma12_inv;
            double term4 = D.gamma1_COM * D.gamma2_COM;
            D.term5 = term4/p2_COM + D.m12;
            D.vrel = D.p_COM/D.term3/term4; // relative velocity
            
            for( unsigned int i=0; i<processes_.size(); i++ ) {
                processes_[i]->apply( patch->rand_, D );
            }
            
        } // end loop on pairs of particles
        
    } // end loop on bins
    
    for( unsigned int i=0; i<processes_.size(); i++ ) {
        processes_[i]->finish( params, patch, localDiags, intra_, species_group1_, species_group2_, itime );
    }
}


void BinaryProcesses::debug( Params &params, int itime, unsigned int icoll, VectorPatch &vecPatches )
{

    int debug_every = vecPatches( 0 )->vecBPs[icoll]->debug_every_;
    if( debug_every > 0 && itime % debug_every == 0 ) {
    
        unsigned int npatch = vecPatches.size();
        
        vector<double> smean( npatch, 0. );
        vector<double> logLmean( npatch, 0. );
        vector<double> debye_length( npatch, 0. );
        vector<double> nuclear_reaction_multiplier( npatch, 0. );
        
        // Collect info for all patches
        for( unsigned int ipatch=0; ipatch<npatch; ipatch++ ) {
            
            // debye length
            unsigned int nbin = vecPatches( ipatch )->debye_length_squared.size();
            for( unsigned int ibin=0; ibin<nbin; ibin++ ) {
                debye_length[ipatch] += vecPatches( ipatch )->debye_length_squared[ibin];
            }
            debye_length[ipatch] = sqrt( debye_length[ipatch] / nbin );
            
            // Data from processes
            vector<BinaryProcess*> vBP = vecPatches( ipatch )->vecBPs[icoll]->processes_;
            for( unsigned int iBP=0; iBP<vBP.size(); iBP++ ) {
                if( Collisions *coll = dynamic_cast<Collisions*>(vBP[iBP]) ) {
                    smean   [ipatch] = coll->smean_   ;
                    logLmean[ipatch] = coll->logLmean_;
                } else if( CollisionalNuclearReaction *NR = dynamic_cast<CollisionalNuclearReaction*>(vBP[iBP]) ) {
                    nuclear_reaction_multiplier[ipatch] = NR->rate_multiplier_;
                }
            }
        }
        
        // Create H5 group for the current timestep
        ostringstream name( "" );
        name << "t" << setfill( '0' ) << setw( 8 ) << itime;
        H5Write g = vecPatches( 0 )->vecBPs[icoll]->debug_file_->group( name.str() );
        // Create new datasets for this timestep and write
        g.vect( "s"                          , smean                      [0], params.tot_number_of_patches, H5T_NATIVE_DOUBLE, vecPatches.refHindex_, npatch );
        g.vect( "coulomb_log"                , logLmean                   [0], params.tot_number_of_patches, H5T_NATIVE_DOUBLE, vecPatches.refHindex_, npatch );
        g.vect( "debyelength"                , debye_length               [0], params.tot_number_of_patches, H5T_NATIVE_DOUBLE, vecPatches.refHindex_, npatch );
        g.vect( "nuclear_reaction_multiplier", nuclear_reaction_multiplier[0], params.tot_number_of_patches, H5T_NATIVE_DOUBLE, vecPatches.refHindex_, npatch );
        vecPatches( 0 )->vecBPs[icoll]->debug_file_->flush();
        
    }
    
}
