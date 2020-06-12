
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <fstream>

#include "Collisions.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "H5.h"
#include "Patch.h"
#include "VectorPatch.h"

using namespace std;


// Constructor
Collisions::Collisions(
    Params &params,
    unsigned int n_collisions,
    vector<unsigned int> species_group1,
    vector<unsigned int> species_group2,
    double coulomb_log,
    bool intra_collisions,
    int debug_every,
    CollisionalIonization *ionization,
    CollisionalNuclearReaction *nuclear_reaction,
    string filename
) :
    Ionization( ionization ),
    NuclearReaction( nuclear_reaction ),
    n_collisions_( n_collisions ),
    species_group1_( species_group1 ),
    species_group2_( species_group2 ),
    coulomb_log_( coulomb_log ),
    intra_collisions_( intra_collisions ),
    debug_every_( debug_every ),
    filename_( filename )
{
    coeff1_ = 4.046650232e-21*params.reference_angular_frequency_SI; // h*omega/(2*me*c^2)
    coeff2_ = 2.817940327e-15*params.reference_angular_frequency_SI/299792458.; // re omega / c
}


// Cloning Constructor
Collisions::Collisions( Collisions *coll )
{

    n_collisions_     = coll->n_collisions_    ;
    species_group1_   = coll->species_group1_  ;
    species_group2_   = coll->species_group2_  ;
    coulomb_log_      = coll->coulomb_log_     ;
    intra_collisions_ = coll->intra_collisions_;
    debug_every_      = coll->debug_every_     ;
    filename_         = coll->filename_        ;
    coeff1_           = coll->coeff1_        ;
    coeff2_           = coll->coeff2_        ;
    
    if( dynamic_cast<CollisionalNoIonization *>( coll->Ionization ) ) {
        Ionization = new CollisionalNoIonization();
    } else {
        Ionization = new CollisionalIonization( coll->Ionization );
    }
    
    if( dynamic_cast<CollisionalNoNuclearReaction *>( coll->NuclearReaction ) ) {
        NuclearReaction = new CollisionalNoNuclearReaction();
    } else if ( dynamic_cast<CollisionalFusionDD *>( coll->NuclearReaction ) ) {
        NuclearReaction = new CollisionalFusionDD( coll->NuclearReaction );
    }
}


Collisions::~Collisions()
{
    delete Ionization;
    delete NuclearReaction;
}

// Declare other static variables here
bool   Collisions::debye_length_required;


// Calculates the debye length squared in each patch
// The formula for the inverse debye length squared is sumOverSpecies(density*charge^2/temperature)
void Collisions::calculate_debye_length( Params &params, Patch *patch )
{
    double p2, density, density_max, charge, temperature, rmin2;
    Species    *s;
    Particles *p;
    double coeff = 299792458./( 3.*params.reference_angular_frequency_SI*2.8179403267e-15 ); // c / (3 omega re)
    
    unsigned int nspec = patch->vecSpecies.size(); // number of species
    if( nspec==0 ) {
        return;
    }
    unsigned int nbin = patch->vecSpecies[0]->particles->first_index.size();
    
    patch->debye_length_squared.resize( nbin, 0. );
    double mean_debye_length = 0.;
    for( unsigned int ibin = 0 ; ibin < nbin ; ibin++ ) {
        density_max = 0.;
        
        double inv_cell_volume = 0.;
        
        for( unsigned int ispec=0 ; ispec<nspec ; ispec++ ) { // loop all species
            s  = patch->vecSpecies[ispec];
            p  = s->particles;
            
            // Skip when no particles
            if( s->particles->last_index[ibin] <= s->particles->first_index[ibin] ) continue;
            
            if( inv_cell_volume == 0. ) {
                inv_cell_volume = 1. / patch->getPrimalCellVolume( p, s->particles->first_index[ibin], params );
            }
            
            // Calculation of particles density, mean charge, and temperature
            // Density is the sum of weights
            // Temperature definition is the average <v*p> divided by 3
            density     = 0.;
            charge      = 0.;
            temperature = 0.;
            // loop particles to calculate average quantities
            for( unsigned int iPart=s->particles->first_index[ibin]; iPart<( unsigned int )s->particles->last_index[ibin] ; iPart++ ) {
                p2 = p->momentum( 0, iPart ) * p->momentum( 0, iPart )
                     +p->momentum( 1, iPart ) * p->momentum( 1, iPart )
                     +p->momentum( 2, iPart ) * p->momentum( 2, iPart );
                density     += p->weight( iPart );
                charge      += p->weight( iPart ) * p->charge( iPart );
                temperature += p->weight( iPart ) * p2/sqrt( 1.+p2 );
            }
            if( density <= 0. ) {
                continue;
            }
            charge /= density; // average charge
            temperature *= s->mass_ / ( 3.*density ); // Te in units of me*c^2
            density *= inv_cell_volume; // density in units of critical density
            // compute inverse debye length squared
            if( temperature>0. ) {
                patch->debye_length_squared[ibin] += density*charge*charge/temperature;
            }
            // compute maximum density of species
            if( density>density_max ) {
                density_max = density;
            }
        }
        
        // if there were particles,
        if( patch->debye_length_squared[ibin] > 0. ) {
            // compute debye length squared in code units
            patch->debye_length_squared[ibin] = 1./( patch->debye_length_squared[ibin] );
            // apply lower limit to the debye length (minimum interatomic distance)
            rmin2 = pow( coeff*density_max, -2./3. );
            if( patch->debye_length_squared[ibin] < rmin2 ) {
                patch->debye_length_squared[ibin] = rmin2;
            }
        }
        
        mean_debye_length += sqrt( patch->debye_length_squared[ibin] );
    }
    
    mean_debye_length /= nbin;
    DEBUG( "Mean Debye length in code length units = " << scientific << setprecision( 3 ) << mean_debye_length );
    mean_debye_length *= 299792458./params.reference_angular_frequency_SI; // switch to SI
    DEBUG( "Mean Debye length in meters = " << scientific << setprecision( 3 ) << mean_debye_length );
}

// Calculates the collisions for a given Collisions object
void Collisions::collide( Params &params, Patch *patch, int itime, vector<Diagnostic *> &localDiags )
{

    vector<unsigned int> *sg1, *sg2, index1, index2;
    unsigned int nspec1, nspec2; // numbers of species in each group
    unsigned int npart1, npart2; // numbers of macro-particles in each group
    unsigned int npairs; // number of pairs of macro-particles
    vector<unsigned int> np1, np2; // numbers of macro-particles in each species, in each group
    unsigned int i1=0, i2, ispec1, ispec2, N2max;
    Species   *s1, *s2;
    Particles *p1=NULL, *p2;
    double coeff3, coeff4, logL, s, ncol, debye2=0.;
    
    sg1 = &species_group1_;
    sg2 = &species_group2_;
    
    
    bool debug = ( debug_every_ > 0 && itime % debug_every_ == 0 ); // debug only every N timesteps
    
    ncol = 0.;
    smean_       = 0.;
    logLmean_    = 0.;
    //temperature = 0.;
    
    NuclearReaction->prepare();
    
    // Loop bins of particles (typically, cells, but may also be clusters)
    unsigned int nbin = patch->vecSpecies[0]->particles->first_index.size();
    for( unsigned int ibin = 0 ; ibin < nbin ; ibin++ ) {
    
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
        
        // Set the debye length
        if( Collisions::debye_length_required ) {
            debye2 = patch->debye_length_squared[ibin];
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
        if( intra_collisions_ ) { // In the case of collisions within one species
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
        } else { // In the case of collisions between two species
            npairs = npart1; // as many pairs as macro-particles in group 1 (most numerous)
            index2.resize( npairs );
            for( unsigned int i=0; i<npart1; i++ ) {
                index2[i] = i % npart2;
            }
            N2max = npart2; // number of not-repeated particles (in group 2 only)
        }
        
        // Prepare the ionization
        Ionization->prepare1( patch->vecSpecies[( *sg1 )[0]]->atomic_number_ );
        
        // Calculate the densities
        double n1  = 0.; // density of group 1
        for( ispec1=0 ; ispec1<nspec1 ; ispec1++ ) {
            s1 = patch->vecSpecies[( *sg1 )[ispec1]];
            p1 = s1->particles;
            for( int i = s1->particles->first_index[ibin]; i < s1->particles->last_index[ibin]; i++ ) {
                n1 += p1->weight( i );
            }
        }
        double n2  = 0.; // density of group 2
        for( ispec2=0 ; ispec2<nspec2 ; ispec2++ ) {
            s2 = patch->vecSpecies[( *sg2 )[ispec2]];
            p2 = s2->particles;
            for( int i = s2->particles->first_index[ibin]; i < s2->particles->last_index[ibin]; i++ ) {
                n2 += p2->weight( i );
            }
        }
        
        // Pre-calculate some numbers before the big loop
        s1 = patch->vecSpecies[( *sg1 )[0]];
        p1 = s1->particles;
        double inv_cell_volume = 1./patch->getPrimalCellVolume( p1, s1->particles->first_index[ibin], params );
        unsigned int ncorr = intra_collisions_ ? 2*npairs-1 : npairs;
        double dt_corr = params.timestep * ((double)ncorr) * inv_cell_volume;
        coeff3 = coeff2_ * dt_corr;
        coeff4 = pow( 3.*coeff2_, -1./3. ) * dt_corr;
        double weight_correction_1 = 1. / (double)( (npairs-1) / N2max );
        double weight_correction_2 = 1. / (double)( (npairs-1) / N2max + 1 );
        n1  *= inv_cell_volume;
        n2  *= inv_cell_volume;
        double n123 = pow( n1, 2./3. );
        double n223 = pow( n2, 2./3. );
        
        // Now start the real loop on pairs of particles
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        for( unsigned int i = 0; i<npairs; i++ ) {
            
            // find species and index i1 of particle "1"
            i1 = index1[i];
            for( ispec1=0 ; i1>=np1[ispec1]; ispec1++ ) {
                i1 -= np1[ispec1];
            }
            // find species and index i2 of particle "2"
            i2 = index2[i];
            for( ispec2=0 ; i2>=np2[ispec2]; ispec2++ ) {
                i2 -= np2[ispec2];
            }
            
            s1 = patch->vecSpecies[( *sg1 )[ispec1]];
            s2 = patch->vecSpecies[( *sg2 )[ispec2]];
            i1 += s1->particles->first_index[ibin];
            i2 += s2->particles->first_index[ibin];
            p1 = s1->particles;
            p2 = s2->particles;
            
            double weight_correction = std::max( p1->weight(i1), p2->weight(i2) );
            if( i % N2max <= (npairs-1) % N2max ) {
                weight_correction *= weight_correction_2 ;
            } else {
                weight_correction *= weight_correction_1;
            }
            
            logL = coulomb_log_;
            double U1  = patch->rand_->uniform();
            double U2  = patch->rand_->uniform();
            double phi = patch->rand_->uniform_2pi();
            
            s = one_collision( p1, i1, s1->mass_, p2, i2, s2->mass_, coeff1_, coeff3*weight_correction, coeff4*weight_correction, n123, n223, debye2, logL, U1, U2, phi );
            
            // Handle ionization
            Ionization->apply( patch, p1, i1, p2, i2, dt_corr*weight_correction );
            
            ncol ++;
            if( debug ) {
                smean_    += s;
                logLmean_ += logL;
                //temperature += m1 * (sqrt(1.+pow(p1->momentum(0,i1),2)+pow(p1->momentum(1,i1),2)+pow(p1->momentum(2,i1),2))-1.);
            }
            
        } // end loop on pairs of particles
        
    } // end loop on bins
    
    Ionization->finish( params, patch, localDiags );
    NuclearReaction->finish( params, patch, localDiags, intra_collisions_, species_group1_, species_group2_, ncol, itime );
    
    if( debug && ncol>0. ) {
        smean_    /= ncol;
        logLmean_ /= ncol;
        //temperature /= ncol;
    }
}


void Collisions::debug( Params &params, int itime, unsigned int icoll, VectorPatch &vecPatches )
{

    int debug_every = vecPatches( 0 )->vecCollisions[icoll]->debug_every_;
    if( debug_every > 0 && itime % debug_every == 0 ) {
    
        unsigned int npatch = vecPatches.size();
        
        //vector<double> ncol(npatch, 0.);
        vector<double> smean( npatch, 0. );
        vector<double> logLmean( npatch, 0. );
        //vector<double>  temperature=(npatch, 0.);
        vector<double> debye_length( npatch, 0. );
        vector<double> nuclear_reaction_multiplier( npatch, 0. );
        
        // Collect info for all patches
        for( unsigned int ipatch=0; ipatch<npatch; ipatch++ ) {
            //ncol       [ipatch] = vecPatches(ipatch)->vecCollisions[icoll]->ncol       ;
            smean      [ipatch] = vecPatches( ipatch )->vecCollisions[icoll]->smean_      ;
            logLmean   [ipatch] = vecPatches( ipatch )->vecCollisions[icoll]->logLmean_   ;
            //temperature[ipatch] = vecPatches(ipatch)->vecCollisions[icoll]->temperature;
            unsigned int nbin = vecPatches( ipatch )->debye_length_squared.size();
            for( unsigned int ibin=0; ibin<nbin; ibin++ ) {
                debye_length[ipatch] += vecPatches( ipatch )->debye_length_squared[ibin];
            }
            debye_length[ipatch] = sqrt( debye_length[ipatch] / nbin );
            nuclear_reaction_multiplier[ipatch] = vecPatches( ipatch )->vecCollisions[icoll]->NuclearReaction->rate_multiplier_;
        }
        
        // Open the HDF5 file
        hid_t file_access = H5Pcreate( H5P_FILE_ACCESS );
        H5Pset_fapl_mpio( file_access, MPI_COMM_WORLD, MPI_INFO_NULL );
        hid_t fileId = H5Fopen( vecPatches( 0 )->vecCollisions[icoll]->filename_.c_str(), H5F_ACC_RDWR, file_access );
        H5Pclose( file_access );
        // Create H5 group for the current timestep
        ostringstream name( "" );
        name << "t" << setfill( '0' ) << setw( 8 ) << itime;
        hid_t group = H5::group( fileId, name.str() );
        // Define the size in memory for this MPI
        hsize_t mem_size[1] = {npatch};
        hid_t memspace  = H5Screate_simple( 1, mem_size, NULL );
        // Define size and location in file
        hsize_t dimsf[1] = {( hsize_t )params.tot_number_of_patches};
        hid_t filespace = H5Screate_simple( 1, dimsf, NULL );
        hsize_t offset[1] = {( hsize_t )vecPatches.refHindex_}, stride[1] = {1}, count[1] = {1}, block[1] = {npatch};
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, stride, count, block );
        // Define transfer
        hid_t transfer = H5Pcreate( H5P_DATASET_XFER );
        H5Pset_dxpl_mpio( transfer, H5FD_MPIO_COLLECTIVE );
        // Define dataset property list
        hid_t plist_id = H5Pcreate( H5P_DATASET_CREATE );
        H5Pset_alloc_time( plist_id, H5D_ALLOC_TIME_EARLY );
        // Create new datasets for this timestep and write
        hid_t dset_id;
        dset_id  = H5Dcreate( group, "s", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &smean[0] );
        H5Dclose( dset_id );
        dset_id  = H5Dcreate( group, "coulomb_log", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &logLmean[0] );
        H5Dclose( dset_id );
        //dset_id  = H5Dcreate(group, "temperature", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        //H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &temperature[0] );
        //H5Dclose(dset_id);
        dset_id  = H5Dcreate( group, "debyelength", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &debye_length[0] );
        H5Dclose( dset_id );
        dset_id  = H5Dcreate( group, "nuclear_reaction_multiplier", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &nuclear_reaction_multiplier[0] );
        H5Dclose( dset_id );
        // Close all
        H5Pclose( plist_id );
        H5Pclose( transfer );
        H5Sclose( filespace );
        H5Sclose( memspace );
        H5Gclose( group );
        H5Fclose( fileId );
    }
    
}
