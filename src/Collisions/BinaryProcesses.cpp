
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <fstream>

#include "BinaryProcesses.h"
#include "Collisions.h"
#include "CollisionalNuclearReaction.h"
#include "CollisionalFusionDD.h"
#include "CollisionalIonization.h"
#include "H5.h"
#include "Patch.h"
#include "VectorPatch.h"
#include "RandomShuffle.h"

using namespace std;


// Constructor
BinaryProcesses::BinaryProcesses(
    Params &params,
    vector<unsigned int> species_group1,
    vector<unsigned int> species_group2,
    bool intra,
    int screening_group,
    CollisionalNuclearReaction * nuclear_reactions,
    double clog,
    double coulomb_log,
    CollisionalIonization * collisional_ionization,
    int every,
    int debug_every,
    double time_frozen,
    string filename
) :
    nuclear_reactions_( nuclear_reactions ),
    collisions_( params, clog, coulomb_log ),
    collisional_ionization_( collisional_ionization ),
    species_group1_( species_group1 ),
    species_group2_( species_group2 ),
    intra_( intra ),
    screening_group_( screening_group ),
    every_( every ),
    debug_every_( debug_every ),
    timesteps_frozen_( time_frozen / params.timestep ),
    filename_( filename ),
    debug_file_( NULL )
{
    // Open the HDF5 file
    if( debug_every > 0 ) {
        MPI_Comm comm = MPI_COMM_WORLD;
        debug_file_ = new H5Write( filename_, &comm );
    }
}


// Cloning Constructor
BinaryProcesses::BinaryProcesses( BinaryProcesses *BPs ) :
    nuclear_reactions_      ( nullptr ),
    collisions_             ( BPs->collisions_ ),
    collisional_ionization_ ( nullptr ),
    species_group1_         ( BPs->species_group1_   ),
    species_group2_         ( BPs->species_group2_   ),
    intra_                  ( BPs->intra_            ),
    screening_group_        ( BPs->screening_group_  ),
    every_                  ( BPs->every_            ),
    debug_every_            ( BPs->debug_every_      ),
    timesteps_frozen_       ( BPs->timesteps_frozen_ ),
    filename_               ( BPs->filename_         ),
    debug_file_             ( BPs->debug_file_       )
{
    if( CollisionalFusionDD * DD = dynamic_cast<CollisionalFusionDD*>( BPs->nuclear_reactions_ ) ) {
        nuclear_reactions_ = new CollisionalFusionDD( DD );
    }
    if( CollisionalIonization * CI = dynamic_cast<CollisionalIonization*>( BPs->collisional_ionization_ ) ) {
        collisional_ionization_ = new CollisionalIonization( CI );
    }
}

BinaryProcesses::~BinaryProcesses()
{
    MESSAGE( nuclear_reactions_ );
    delete nuclear_reactions_;
    delete collisional_ionization_;
}

// Declare other static variables here
bool BinaryProcesses::debye_length_required_;

// Calculates the debye length squared in each bin
// The formula for the inverse debye length squared is sumOverSpecies(density*charge^2/temperature)
void BinaryProcesses::calculate_debye_length( Params &params, Patch *patch )
{
    const size_t nspec = patch->vecSpecies.size(); // number of species
    if( nspec == 0 ) {
        return;
    }
    
    CellVolumeCalculator cellVolume( params, patch );
    const size_t nbin = cellVolume.nbin_;
    
    patch->debye_length_squared.resize( nbin );
    
    const double coeff = 299792458./( 3.*params.reference_angular_frequency_SI*2.8179403267e-15 ); // c / (3 omega re)
    
#ifdef  __DEBUG
    double mean_debye_length = 0.;
#endif
    
    // Make pointers to all the necessary data (needed for GPU)
    double * __restrict__ debye2_ptr = patch->debye_length_squared.data();
    int *__restrict__ last_index_ptr[nspec];
    double *__restrict__ weight_ptr[nspec];
    short *__restrict__ charge_ptr[nspec];
    double *__restrict__ px_ptr[nspec];
    double *__restrict__ py_ptr[nspec];
    double *__restrict__ pz_ptr[nspec];
    double mass[nspec];
    for( size_t ispec = 0; ispec < nspec; ispec ++ ) {
        Particles *p = patch->vecSpecies[ispec]->particles;
        last_index_ptr[ispec] = p->getPtrCellLastIndex();
        weight_ptr[ispec] = p->getPtrWeight();
        charge_ptr[ispec] = p->getPtrCharge();
        px_ptr[ispec] = p->getPtrMomentum( 0 );
        py_ptr[ispec] = p->getPtrMomentum( 1 );
        pz_ptr[ispec] = p->getPtrMomentum( 2 );
        mass[ispec] = patch->vecSpecies[ispec]->mass_;
    }
    
    // Loop bins of particles
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel loop gang worker num_workers(2) vector_length(32) \
        copyout( debye2_ptr[:nbin] ) \
        copyin( cellVolume, nspec, mass[:nspec], last_index_ptr[:nspec], \
            weight_ptr[:nspec], charge_ptr[:nspec], px_ptr[:nspec], py_ptr[:nspec], pz_ptr[:nspec] )
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target teams distribute thread_limit(32)\
        map( from: debye2_ptr[:nbin] ) \
        map( to: cellVolume, nspec, mass[:nspec], last_index_ptr[:nspec], \
            weight_ptr[:nspec], charge_ptr[:nspec], px_ptr[:nspec], py_ptr[:nspec], pz_ptr[:nspec] )
#endif
    for( size_t ibin = 0 ; ibin < nbin ; ibin++ ) {
        double density_tot = 0.;
        double inv_D2 = 0.;
        double inv_cell_volume = 1. / cellVolume( ibin );
        
        for( size_t ispec = 0 ; ispec < nspec ; ispec++ ) { // loop all species
            // Calculation of particles density, mean charge, and temperature
            // Density is the sum of weights
            // Temperature definition is the average <v*p> divided by 3
            double density     = 0.;
            double charge      = 0.;
            double temperature = 0.;
            
            // loop particles to calculate average quantities
            int last = last_index_ptr[ispec][ibin];
            int first = ibin == 0 ? 0 : last_index_ptr[ispec][ibin-1];
            #if defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop vector reduction(+: density, charge, temperature)
            #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp parallel for reduction(+: density, charge, temperature)
            #endif
            for( int iPart = first; iPart < last; iPart++ ) {
                double p2 =
                    px_ptr[ispec][iPart] * px_ptr[ispec][iPart] +
                    py_ptr[ispec][iPart] * py_ptr[ispec][iPart] +
                    pz_ptr[ispec][iPart] * pz_ptr[ispec][iPart];
                density     += weight_ptr[ispec][iPart];
                charge      += weight_ptr[ispec][iPart] * charge_ptr[ispec][iPart];
                temperature += weight_ptr[ispec][iPart] * p2 / sqrt( 1.+p2 );
            }
            
            if( density > 0. ) {
                charge /= density; // average charge
                temperature *= mass[ispec] / ( 3.*density ); // Te in units of me*c^2
                density *= inv_cell_volume; // density in units of critical density
                // compute inverse debye length squared
                if( temperature == 0. ) {
                    inv_D2 += 1e100; // infinite
                } else {
                    inv_D2 += density*charge*charge/temperature;
                }
                // compute total density of species
                density_tot += density;
            }
        }
        
        // if there were particles,
        if( inv_D2 > 0. ) {
            // compute debye length squared in code units
            debye2_ptr[ibin] = 1./inv_D2;
            // apply lower limit to the debye length (minimum interatomic distance)
            double rmin2 = pow( coeff*density_tot, -2./3. );
            if( debye2_ptr[ibin] < rmin2 ) {
                debye2_ptr[ibin] = rmin2;
            }
        } else {
            debye2_ptr[ibin] = 0.;
        }
#ifdef  __DEBUG
        mean_debye_length += sqrt( patch->debye_length_squared[ibin] );
#endif
    }

#ifdef  __DEBUG
    mean_debye_length /= nbin;
    DEBUG( "Mean Debye length in code length units = " << scientific << setprecision( 3 ) << mean_debye_length );
    mean_debye_length *= 299792458./params.reference_angular_frequency_SI; // switch to SI
    DEBUG( "Mean Debye length in meters = " << scientific << setprecision( 3 ) << mean_debye_length );
#endif
}


// Make the pairing and launch the processes
void BinaryProcesses::apply( Params &params, Patch *patch, int itime, vector<Diagnostic *> &localDiags )
{
    if( itime < timesteps_frozen_ || itime % every_ != 0 ) {
        return;
    }
    
    if( nuclear_reactions_ ) {
        nuclear_reactions_->prepare();
    }
    if( collisions_ ) {
        collisions_.prepare();
    }
    
    // Time between binary process events
    const double delta_t = every_ * params.timestep;
    
    // numbers of species in each group
    const size_t nspec1 = species_group1_.size();
    const size_t nspec2 = species_group2_.size();
    
    // Info for ionization
    const bool electronFirst = patch->vecSpecies[species_group1_[0]]->atomic_number_==0 ? true : false;
    
    // Store info for screening (e-i collisions)
    vector<double> screening_Z; // atomic number
    vector<double> lTF; // thomas-fermi length
    if( screening_group_ > 0 ) {
        auto &sg = screening_group_ == 1 ? species_group1_ : species_group2_;
        screening_Z.resize( sg.size() );
        lTF.resize( sg.size() );
        for( size_t i = 0; i < sg.size(); i++ ) {
            screening_Z[i] = (double) patch->vecSpecies[sg[i]]->atomic_number_;
            lTF[i] = 3.1255e-19 // (9*pi^2/16)^(1/3) * a0 /c
                *params.reference_angular_frequency_SI * pow( screening_Z[i], -1/3 );
        }
    }
    screening_Z.push_back( 0. ); // placeholder for no screening
    lTF.push_back( 0. );
    size_t screening_group_size = screening_Z.size();
    
    // Make pointers to all the necessary data (needed for GPU)
    const double *const __restrict__ debye2_ptr = patch->debye_length_squared.data();
    const double *const __restrict__ screening_Z_ptr = screening_Z.data();
    const double *const __restrict__ lTF_ptr = lTF.data();
    // const unsigned int *__restrict__ sg1_ptr = species_group1_.data();
    // const unsigned int *__restrict__ sg2_ptr = species_group2_.data();
    Particles *__restrict__ p1_ptr[nspec1], *__restrict__ p2_ptr[nspec2];
    int *__restrict__ last_index1_ptr[nspec1], *__restrict__ last_index2_ptr[nspec2];
    double *__restrict__ weight1_ptr[nspec1], *__restrict__ weight2_ptr[nspec2];
    short *__restrict__ charge1_ptr[nspec1], *__restrict__ charge2_ptr[nspec2];
    double *__restrict__ px1_ptr[nspec1], *__restrict__ px2_ptr[nspec2];
    double *__restrict__ py1_ptr[nspec1], *__restrict__ py2_ptr[nspec2];
    double *__restrict__ pz1_ptr[nspec1], *__restrict__ pz2_ptr[nspec2];
    double mass1[nspec1], mass2[nspec2];
    for( size_t ispec1 = 0; ispec1 < nspec1; ispec1 ++ ) {
        p1_ptr[ispec1] = patch->vecSpecies[species_group1_[ispec1]]->particles;
        last_index1_ptr[ispec1] = p1_ptr[ispec1]->getPtrCellLastIndex();
        weight1_ptr[ispec1] = p1_ptr[ispec1]->getPtrWeight();
        charge1_ptr[ispec1] = p1_ptr[ispec1]->getPtrCharge();
        px1_ptr[ispec1] = p1_ptr[ispec1]->getPtrMomentum( 0 );
        py1_ptr[ispec1] = p1_ptr[ispec1]->getPtrMomentum( 1 );
        pz1_ptr[ispec1] = p1_ptr[ispec1]->getPtrMomentum( 2 );
        mass1[ispec1] = patch->vecSpecies[species_group1_[ispec1]]->mass_;
    }
    for( size_t ispec2 = 0; ispec2 < nspec2; ispec2 ++ ) {
        p2_ptr[ispec2] = patch->vecSpecies[species_group2_[ispec2]]->particles;
        last_index2_ptr[ispec2] = p2_ptr[ispec2]->getPtrCellLastIndex();
        weight2_ptr[ispec2] = p2_ptr[ispec2]->getPtrWeight();
        charge2_ptr[ispec2] = p2_ptr[ispec2]->getPtrCharge();
        px2_ptr[ispec2] = p2_ptr[ispec2]->getPtrMomentum( 0 );
        py2_ptr[ispec2] = p2_ptr[ispec2]->getPtrMomentum( 1 );
        pz2_ptr[ispec2] = p2_ptr[ispec2]->getPtrMomentum( 2 );
        mass2[ispec2] = patch->vecSpecies[species_group2_[ispec2]]->mass_;
    }
    
    BinaryProcessData D;
    D.screening_group = screening_group_;
    D.electronFirst = electronFirst;
    CellVolumeCalculator cellVolume( params, patch );
    const uint32_t nbin = cellVolume.nbin_;
    SMILEI_ASSERT_VERBOSE( cellVolume.nbin_ < 4294967296, "Too many cells per patch in Collisions" );
    
    // Due to GPU offloading, we must have a different rand object for each bin
    // Each bin gets a seed equal to the `patch->rand_` seed plus `ibin`
    Random rand( patch->rand_ );
    RandomShuffle shuffler( rand, 1 );
    patch->rand_->add( nbin );
    
    
    
    // Loop bins of particles
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    uint32_t np1[nspec1], np2[nspec2];
    #pragma acc parallel loop gang vector_length(32) firstprivate(rand) private(D, np1, np2, shuffler)\
        copyin( cellVolume, delta_t, nspec1, nspec2,/* sg1_ptr[:nspec1], sg2_ptr[:nspec2],*/ \
            screening_group_size, screening_Z_ptr[:screening_group_size], lTF_ptr[:screening_group_size], \
            mass1[:nspec1], mass2[:nspec2], debye2_ptr[:nbin], \
            last_index1_ptr[:nspec1], weight1_ptr[:nspec1], charge1_ptr[:nspec1], px1_ptr[:nspec1], py1_ptr[:nspec1], pz1_ptr[:nspec1], \
            last_index2_ptr[:nspec2], weight2_ptr[:nspec2], charge2_ptr[:nspec2], px2_ptr[:nspec2], py2_ptr[:nspec2], pz2_ptr[:nspec2] )
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    uint32_t np1[10], np2[10]; // The cray compiler crashes if these arrays have non-fixed sizes 
    SMILEI_ASSERT_VERBOSE( nspec1 <= 10 && nspec2 <= 10, "Too many species in Collisions" );
    #pragma omp target teams distribute thread_limit(32) firstprivate(rand, shuffler) private(D, np1, np2) \
        map( to: cellVolume, delta_t, nspec1, nspec2,/* sg1_ptr[:nspec1], sg2_ptr[:nspec2],*/ \
            screening_group_size, screening_Z_ptr[:screening_group_size], lTF_ptr[:screening_group_size], \
            mass1[:nspec1], mass2[:nspec2], debye2_ptr[:nbin], \
            last_index1_ptr[:nspec1], weight1_ptr[:nspec1], charge1_ptr[:nspec1], px1_ptr[:nspec1], py1_ptr[:nspec1], pz1_ptr[:nspec1], \
            last_index2_ptr[:nspec2], weight2_ptr[:nspec2], charge2_ptr[:nspec2], px2_ptr[:nspec2], py2_ptr[:nspec2], pz2_ptr[:nspec2] )
#else
    uint32_t np1[nspec1], np2[nspec2];
#endif
    for( uint32_t ibin = 0 ; ibin < nbin ; ibin++ ) {
        
        // get number of particles for all necessary species
        uint32_t npart1 = 0;
        for( uint32_t ispec1=0 ; ispec1<nspec1 ; ispec1++ ) {
            uint32_t last = last_index1_ptr[ispec1][ibin];
            uint32_t first = ( ibin == 0 ) ? 0 : last_index1_ptr[ispec1][ibin-1];
            np1[ispec1] = last - first;
            npart1 += np1[ispec1];
        }
        uint32_t npart2 = 0;
        for( uint32_t ispec2=0 ; ispec2<nspec2 ; ispec2++ ) {
            uint32_t last = last_index2_ptr[ispec2][ibin];
            uint32_t first = ( ibin == 0 ) ? 0 : last_index2_ptr[ispec2][ibin-1];
            np2[ispec2] = last - first;
            npart2 += np2[ispec2];
        }
        // We need to shuffle the group that has most particles
        bool shuffle1 = npart1 > npart2;
        uint32_t npartmin = npart1;
        uint32_t npartmax = npart2;
        if( shuffle1 ) {
            npartmin = npart2;
            npartmax = npart1;
        }
        
        // skip to next bin if not enough pairs
        if( npartmin == 0 || ( intra_ && npartmin < 2 ) ) {
            continue;
        }
        
        uint32_t npairs, npairs_not_repeated;
        SMILEI_BINARYPROCESS_FLOAT weight_correction_1, weight_correction_2;
        if( intra_ ) { // In the case of pairing within one species
            npairs = ( npartmax + 1 ) / 2; // half as many pairs as macro-particles
            npairs_not_repeated = npartmax - npairs;
            weight_correction_1 = 1.;
            weight_correction_2 = 0.5;
        } else { // In the case of pairing between two species
            npairs = npartmax; // as many pairs as macro-particles in group with more particles
            npairs_not_repeated = npartmin;
            weight_correction_1 = 1. / (SMILEI_BINARYPROCESS_FLOAT)( npairs / npairs_not_repeated );
            weight_correction_2 = 1. / (SMILEI_BINARYPROCESS_FLOAT)( npairs / npairs_not_repeated + 1 );
        }
        
        // Set the shuffler seed and size
        rand.add( ibin );
        shuffler.reinit( rand, npartmax );
        
        // Calculate the densities
        SMILEI_BINARYPROCESS_FLOAT n1  = 0., n2 = 0.;
        #if defined( SMILEI_ACCELERATOR_GPU_OACC )
        #pragma acc loop vector reduction( +:n1 )
        #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp parallel for reduction( +:n1 )
        #endif
        for( uint32_t ispec1=0 ; ispec1<nspec1 ; ispec1++ ) {
            uint32_t start = ibin == 0 ? 0 : last_index1_ptr[ispec1][ibin-1];
            for( uint32_t i = start; i < last_index1_ptr[ispec1][ibin]; i++ ) {
                n1 += weight1_ptr[ispec1][i];
            }
        }
        #if defined( SMILEI_ACCELERATOR_GPU_OACC )
        #pragma acc loop vector reduction( +:n2 )
        #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp parallel for reduction( +:n2 )
        #endif
        for( uint32_t ispec2=0 ; ispec2<nspec2 ; ispec2++ ) {
            uint32_t start = ibin == 0 ? 0 : last_index2_ptr[ispec2][ibin-1];
            for( uint32_t i = start; i < last_index2_ptr[ispec2][ibin]; i++ ) {
                n2 += weight2_ptr[ispec2][i];
            }
        }
        
        // Set the debye length
        if( BinaryProcesses::debye_length_required_ ) {
            D.debye = sqrt( debye2_ptr[ibin] );
        }
        
        // Pre-calculate some numbers before the big loop
        const SMILEI_BINARYPROCESS_FLOAT inv_cell_volume = 1. / cellVolume( ibin );
        const uint32_t ncorr = intra_ ? 2*npairs-1 : npairs;
        const SMILEI_BINARYPROCESS_FLOAT dt_corr = delta_t * ((SMILEI_BINARYPROCESS_FLOAT)ncorr) * inv_cell_volume;
        n1  *= inv_cell_volume;
        n2  *= inv_cell_volume;
        D.n123 = cbrt( n1 * n1 );
        D.n223 = cbrt( n2 * n2 );
        
        // Prepare buffers
        const uint32_t buffer_size = ( npairs_not_repeated < SMILEI_BINARYPROCESS_BUFFERSIZE ) ? npairs_not_repeated : SMILEI_BINARYPROCESS_BUFFERSIZE;
        const uint32_t nbuffers = ( npairs - 1 ) / buffer_size + 1;
        
        // Now start the real loop on pairs of particles
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        
        // Loop on buffers
        for( uint32_t ibuffer = 0; ibuffer < nbuffers; ibuffer++ ) {
            
            const uint32_t start = ibuffer * buffer_size;
            const uint32_t n = ( npairs < start + buffer_size ) ? ( uint32_t ) ( npairs - start ) : buffer_size;
            
            // Determine the shuffled indices in the whole groups of species
            if( intra_ ) {
                shuffler.next( n, &D.i[0][0] );
                shuffler.next( n, &D.i[1][0] );
            } else if( shuffle1 ) {
                shuffler.next( n, &D.i[0][0] );
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    D.i[1][i] = ( i + start ) % npart2;
                }
            } else {
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    D.i[0][i] = ( i + start ) % npart1;
                }
                shuffler.next( n, &D.i[1][0] );
            }
            
            // find species and indices of particles
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                for( D.ispec[0][i] = 0; D.i[0][i]>=np1[D.ispec[0][i]]; D.ispec[0][i]++ ) {
                    D.i[0][i] -= np1[D.ispec[0][i]];
                }
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                for( D.ispec[1][i] = 0; D.i[1][i]>=np2[D.ispec[1][i]]; D.ispec[1][i]++ ) {
                    D.i[1][i] -= np2[D.ispec[1][i]];
                }
            }
            
            // Get screening length & Z
            if( screening_group_ == 0 ) {
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    D.lTF[i] = lTF_ptr[screening_group_size - 1];
                    D.Z1Z2[i] = screening_Z_ptr[screening_group_size - 1];
                }
            } else if( screening_group_ == 1 ) {
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    D.lTF[i] = lTF_ptr[D.ispec[0][i]];
                    D.Z1Z2[i] = screening_Z_ptr[D.ispec[0][i]];
                }
            } else if( screening_group_ == 2 ) {
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    D.lTF[i] = lTF_ptr[D.ispec[1][i]];
                    D.Z1Z2[i] = screening_Z_ptr[D.ispec[1][i]];
                }
            }
            
            // Get particle indices in this bin
            if( ibin > 0 ) {
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    D.i[0][i] += last_index1_ptr[D.ispec[0][i]][ibin-1];
                }
                SMILEI_ACCELERATOR_LOOP_VECTOR
                for( uint32_t i = 0; i<n; i++ ) {
                    D.i[1][i] += last_index2_ptr[D.ispec[1][i]][ibin-1];
                }
            }
#ifndef SMILEI_ACCELERATOR_GPU
            // Get pointers to Particles
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.p[0][i] = p1_ptr[D.ispec[0][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.p[1][i] = p2_ptr[D.ispec[1][i]];
            }
#endif
            // Get masses
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.m[0][i] = mass1[D.ispec[0][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.m[1][i] = mass2[D.ispec[1][i]];
            }
            // Get Weights
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.W[0][i] = weight1_ptr[D.ispec[0][i]][D.i[0][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.W[1][i] = weight2_ptr[D.ispec[1][i]][D.i[1][i]];
            }
            // Get charges
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.q[0][i] = charge1_ptr[D.ispec[0][i]][D.i[0][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.q[1][i] = charge2_ptr[D.ispec[1][i]][D.i[1][i]];
            }
            // Get momenta
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.px[0][i] = px1_ptr[D.ispec[0][i]][D.i[0][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.px[1][i] = px2_ptr[D.ispec[1][i]][D.i[1][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.py[0][i] = py1_ptr[D.ispec[0][i]][D.i[0][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.py[1][i] = py2_ptr[D.ispec[1][i]][D.i[1][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.pz[0][i] = pz1_ptr[D.ispec[0][i]][D.i[0][i]];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.pz[1][i] = pz2_ptr[D.ispec[1][i]][D.i[1][i]];
            }
            
            // Calculate the timestep correction
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.dt_correction[i] = ( D.W[0][i] > D.W[1][i] ? D.W[0][i] : D.W[1][i] ) * dt_corr;
                double corr2 = ( i + start ) % npairs_not_repeated < npairs % npairs_not_repeated;
                double corr1 = 1. - corr2;
                D.dt_correction[i] *= corr1 * weight_correction_1 + corr2 * weight_correction_2;
            }
            
            // Calculate gammas
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.gamma[0][i] = sqrt( 1 + D.px[0][i]*D.px[0][i] + D.py[0][i]*D.py[0][i] + D.pz[0][i]*D.pz[0][i] );
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.gamma[1][i] = sqrt( 1 + D.px[1][i]*D.px[1][i] + D.py[1][i]*D.py[1][i] + D.pz[1][i]*D.pz[1][i] );
            }
            
            // Calculate the mass ratio
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.R[i] = D.m[1][i] / D.m[0][i];
            }
            
            // Calculate the total gamma
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.gamma_tot[i] = D.gamma[0][i] + D.R[i] * D.gamma[1][i];
            }
            
            // Calculate the total momentum
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.px_tot[i] = D.px[0][i] + D.R[i] * D.px[1][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.py_tot[i] = D.py[0][i] + D.R[i] * D.py[1][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.pz_tot[i] = D.pz[0][i] + D.R[i] * D.pz[1][i];
            }
            
            // Calculate the Lorentz invariant gamma1 gamma2 - u1.u2
            // It is equal to the gamma of one particle in the rest frame of the other particle
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.gamma0[i] = D.gamma[0][i] * D.gamma[1][i] - D.px[0][i] * D.px[1][i] - D.py[0][i] * D.py[1][i] - D.pz[0][i] * D.pz[1][i];
            }
            
            // Now we calculate quantities in the center-of-mass frame
            // denoted by the suffix _COM
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                D.gamma_tot_COM[i] = sqrt( 2*D.R[i]*D.gamma0[i] + D.R[i] * D.R[i] + 1 );
                D.gamma_COM0[i] = ( D.R[i] * D.gamma0[i] + 1 ) / D.gamma_tot_COM[i];
            }
            
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                SMILEI_BINARYPROCESS_FLOAT gg = ( D.gamma[0][i] + D.gamma_COM0[i] ) / ( D.gamma_tot[i] + D.gamma_tot_COM[i] );
                D.px_COM[i] = D.px[0][i] - gg * D.px_tot[i];
                D.py_COM[i] = D.py[0][i] - gg * D.py_tot[i];
                D.pz_COM[i] = D.pz[0][i] - gg * D.pz_tot[i];
                D.p_COM[i] = sqrt( D.px_COM[i]*D.px_COM[i] + D.py_COM[i]*D.py_COM[i] + D.pz_COM[i]*D.pz_COM[i] );
            }
            
            // Calculate some intermediate quantities
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                SMILEI_BINARYPROCESS_FLOAT p_gamma_COM = D.p_COM[i] * D.gamma_tot_COM[i];
                D.vrel[i] = p_gamma_COM / ( D.gamma_COM0[i] * ( D.gamma_tot_COM[i] - D.gamma_COM0[i] ) ); // | v2_COM - v1_COM |
            }
            
            // Apply all processes (collisions, ionization, ...)
            #ifndef SMILEI_ACCELERATOR_GPU
            if( nuclear_reactions_ ) {
                nuclear_reactions_->apply( &rand, D, n );
            }
            #endif
            if( collisions_ ) {
                collisions_.apply( &rand, D, n );
            }
            #ifndef SMILEI_ACCELERATOR_GPU
            if( collisional_ionization_ ) {
                collisional_ionization_->apply( &rand, D, n );
            }
            #endif
            
            // Update the particle arrays from the buffers
            // Store Weights
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                weight1_ptr[D.ispec[0][i]][D.i[0][i]] = D.W[0][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                weight2_ptr[D.ispec[1][i]][D.i[1][i]] = D.W[1][i];
            }
            // Store charges
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                charge1_ptr[D.ispec[0][i]][D.i[0][i]] = D.q[0][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                charge2_ptr[D.ispec[1][i]][D.i[1][i]] = D.q[1][i];
            }
            // Store momenta
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                px1_ptr[D.ispec[0][i]][D.i[0][i]] = D.px[0][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                px2_ptr[D.ispec[1][i]][D.i[1][i]] = D.px[1][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                py1_ptr[D.ispec[0][i]][D.i[0][i]] = D.py[0][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                py2_ptr[D.ispec[1][i]][D.i[1][i]] = D.py[1][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                pz1_ptr[D.ispec[0][i]][D.i[0][i]] = D.pz[0][i];
            }
            SMILEI_ACCELERATOR_LOOP_VECTOR
            for( uint32_t i = 0; i<n; i++ ) {
                pz2_ptr[D.ispec[1][i]][D.i[1][i]] = D.pz[1][i];
            }
            
        } // end loop on buffers of particles
        
    } // end loop on bins
    
    // The finishing touch on all processes
    if( nuclear_reactions_ ) {
        nuclear_reactions_->finish( params, patch, localDiags, intra_, species_group1_, species_group2_, itime );
    }
    if( collisions_ ) {
        collisions_.finish( params, patch, localDiags, intra_, species_group1_, species_group2_, itime );
    }
    if( collisional_ionization_ ) {
        collisional_ionization_->finish( params, patch, localDiags, intra_, species_group1_, species_group2_, itime );
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
        for( size_t ipatch=0; ipatch<npatch; ipatch++ ) {

            // debye length
            size_t nbin = vecPatches( ipatch )->debye_length_squared.size();
            for( size_t ibin=0; ibin<nbin; ibin++ ) {
                debye_length[ipatch] += vecPatches( ipatch )->debye_length_squared[ibin];
            }
            debye_length[ipatch] = sqrt( debye_length[ipatch] / nbin );
            
            // Data from processes
            BinaryProcesses * BPs = vecPatches( ipatch )->vecBPs[icoll];
            if( BPs->collisions_ ) {
                smean   [ipatch] = BPs->collisions_.smean_   ;
                logLmean[ipatch] = BPs->collisions_.logLmean_;
            } else if( BPs->nuclear_reactions_ ) {
                nuclear_reaction_multiplier[ipatch] = BPs->nuclear_reactions_->rate_multiplier_;
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
