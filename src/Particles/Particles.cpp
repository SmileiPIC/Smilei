// -----------------------------------------------------------------------------
//
//! \file Particles.cpp
//
//! \brief contains the Particles class methods
//
// -----------------------------------------------------------------------------

#include "Particles.h"

#include <cstring>
#include <iostream>

#include "Params.h"
#include "Patch.h"
#include "Species.h"

#include "Particle.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
#include <numeric>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Particle
// ---------------------------------------------------------------------------------------------------------------------
Particles::Particles():
    tracked( false )
{
    Position.resize( 0 );
    Position_old.resize( 0 );
    Momentum.resize( 0 );
    cell_keys.resize( 0 );
    is_test = false;
    isQuantumParameter = false;
    isMonteCarlo = false;

    double_prop_.resize( 0 );
    short_prop_.resize( 0 );
    uint64_prop_.resize( 0 );
}

Particles::~Particles()
{
    clear();
    shrinkToFit();
}

// ---------------------------------------------------------------------------------------------------------------------
// Create nParticles null particles of nDim size
// ---------------------------------------------------------------------------------------------------------------------
void Particles::initialize( unsigned int nParticles, unsigned int nDim, bool keep_position_old )
{
    // if (nParticles > Weight.capacity()) {
    //     WARNING("You should increase c_part_max in specie namelist");
    // }
    // if( Weight.size() == 0 ) {
    //     static constexpr float kGrowthFactor = 1.2F;
    //     // float c_part_max = part.c_part_max;
    //     // float c_part_max = params.species_param[0].c_part_max;
    //     reserve( static_cast<unsigned int>( std::round( kGrowthFactor * static_cast<float>( nParticles ) ) ),
    //              nDim,
    //              keep_position_old );
    // }

    resize( nParticles, nDim, keep_position_old );

    if( double_prop_.empty() ) { // do this just once

        Position.resize( nDim );
        for( unsigned int i = 0; i < nDim; i++ ) {
            double_prop_.push_back( &( Position[i] ) );
        }

        for( unsigned int i = 0; i < 3; i++ ) {
            double_prop_.push_back( &( Momentum[i] ) );
        }

        double_prop_.push_back( &Weight );

        if( keep_position_old ) {
            Position_old.resize( nDim );
            for( unsigned int i = 0; i < nDim; i++ ) {
                double_prop_.push_back( &( Position_old[i] ) );
            }
        }

        short_prop_.push_back( &Charge );
        if( tracked ) {
            uint64_prop_.push_back( &Id );
        }

        // Quantum parameter (for QED effects):
        // - if radiation reaction (continuous or discontinuous)
        // - if multiphoton-Breit-Wheeler if photons
        if( isQuantumParameter ) {
            double_prop_.push_back( &Chi );
        }

        // Optical Depth for Monte-Carlo processes:
        // - if the discontinuous (Monte-Carlo) radiation reaction
        // are activated, tau is the incremental optical depth to emission
        if( isMonteCarlo ) {
            double_prop_.push_back( &Tau );
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// copy properties from another Particles
// ---------------------------------------------------------------------------------------------------------------------
void Particles::initialize( unsigned int nParticles, Particles &part )
{
    is_test=part.is_test;

    tracked=part.tracked;

    isQuantumParameter=part.isQuantumParameter;

    isMonteCarlo=part.isMonteCarlo;

    initialize( nParticles, part.Position.size(), part.Position_old.size() > 0 );
}

// ---------------------------------------------------------------------------------------------------------------------
//! Set capacity of Particles vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::reserve( unsigned int reserved_particles,
                         unsigned int nDim,
                         bool         keep_position_old )
{
    Position.resize( nDim );

    for( unsigned int i = 0; i < nDim; i++ ) {
        Position[i].reserve( reserved_particles );
    }

    if( keep_position_old ) {
        Position_old.resize( nDim );

        for( unsigned int i = 0; i < nDim; i++ ) {
            Position_old[i].reserve( reserved_particles );
        }
    }

    Momentum.resize( 3 );
    for( unsigned int i = 0; i < 3; i++ ) {
        Momentum[i].reserve( reserved_particles );
    }

    Weight.reserve( reserved_particles );
    Charge.reserve( reserved_particles );

    if( tracked ) {
        Id.reserve( reserved_particles );
    }

    if( isQuantumParameter ) {
        Chi.reserve( reserved_particles );
    }

    if( isMonteCarlo ) {
        Tau.reserve( reserved_particles );
    }

    cell_keys.reserve( reserved_particles );
}

// ---------------------------------------------------------------------------------------------------------------------
//! Set capacity of Particles vectors and keep dimensionality
// ---------------------------------------------------------------------------------------------------------------------
void Particles::reserve( unsigned int reserved_particles )
{
    reserve( reserved_particles, Position.size(), Position_old.size() > 0 );
}

// ---------------------------------------------------------------------------------------------------------------------
//! Initialize like Particles object part with 0 particles and reserve space for reserved_particles particles
// ---------------------------------------------------------------------------------------------------------------------
void Particles::initializeReserve( unsigned int npart_max, Particles &part )
{
    initialize( 0, part );
    reserve( npart_max, part.dimension() );
}

// ---------------------------------------------------------------------------------------------------------------------
//! Resize Particle vectors and change dimensionality according to nDim
//
void Particles::resize( unsigned int nParticles,
                        unsigned int nDim,
                        bool         keep_position_old )
{
    Position.resize( nDim );

    for( unsigned int i = 0; i < nDim; i++ ) {
        Position[i].resize( nParticles, 0. );
    }

    if( keep_position_old ) {
        Position_old.resize( nDim );

        for( unsigned int i = 0; i < nDim; i++ ) {
            Position_old[i].resize( nParticles, 0. );
        }
    }

    Momentum.resize( 3 );

    for( unsigned int i = 0; i < 3; i++ ) {
        Momentum[i].resize( nParticles, 0. );
    }

    Weight.resize( nParticles, 0. );
    Charge.resize( nParticles, 0 );

    if( tracked ) {
        Id.resize( nParticles, 0 );
    }

    if( isQuantumParameter ) {
        Chi.resize( nParticles, 0. );
    }

    if( isMonteCarlo ) {
        Tau.resize( nParticles, 0. );
    }

    cell_keys.resize( nParticles, 0. );
}

// ---------------------------------------------------------------------------------------------------------------------
// Resize Particle vectors with nParticles
// ---------------------------------------------------------------------------------------------------------------------
void Particles::resize( unsigned int nParticles)
{

    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] ).resize( nParticles, 0. );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] ).resize( nParticles, 0 );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] ).resize( nParticles, 0 );
    }

    cell_keys.resize( nParticles, 0. );

}

// ---------------------------------------------------------------------------------------------------------------------
//! Resize the cell_keys vector only
// ---------------------------------------------------------------------------------------------------------------------
void Particles::resizeCellKeys(unsigned int nParticles)
{
    cell_keys.resize( nParticles, 0. );
}

// ---------------------------------------------------------------------------------------------------------------------
//! Remove extra capacity of Particles vectors
//! params [in] compute_cell_keys: if true, cell_keys is affected (default is false)
// ---------------------------------------------------------------------------------------------------------------------
void Particles::shrinkToFit(const bool compute_cell_keys)
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        std::vector<double>( *double_prop_[iprop] ).swap( *double_prop_[iprop] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        std::vector<short>( *short_prop_[iprop] ).swap( *short_prop_[iprop] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        std::vector<uint64_t>( *uint64_prop_[iprop] ).swap( *uint64_prop_[iprop] );
    }

    if (compute_cell_keys) {
        cell_keys.swap(cell_keys);
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! Reset of Particles vectors
//! params [in] compute_cell_keys: if true, cell_keys is affected (default is false)
// ---------------------------------------------------------------------------------------------------------------------
void Particles::clear(const bool compute_cell_keys)
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        double_prop_[iprop]->clear();
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        short_prop_[iprop]->clear();
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        uint64_prop_[iprop]->clear();
    }

    if (compute_cell_keys) {
        cell_keys.clear();
    }

}

//! copy particles ipart at the end of the particle vector
//! cell keys not affected
void Particles::copyParticle( unsigned int ipart )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        double_prop_[iprop]->push_back( ( *double_prop_[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        short_prop_[iprop]->push_back( ( *short_prop_[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        uint64_prop_[iprop]->push_back( ( *uint64_prop_[iprop] )[ipart] );
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! Copy particle iPart at the end of dest_parts
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::copyParticle( unsigned int ipart, Particles &dest_parts )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        dest_parts.double_prop_[iprop]->push_back( ( *double_prop_[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        dest_parts.short_prop_[iprop]->push_back( ( *short_prop_[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        dest_parts.uint64_prop_[iprop]->push_back( ( *uint64_prop_[iprop] )[ipart] );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Insert particle iPart at dest_id in dest_parts
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::copyParticle( unsigned int ipart, Particles &dest_parts, int dest_id )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        dest_parts.double_prop_[iprop]->insert( dest_parts.double_prop_[iprop]->begin() + dest_id, ( *double_prop_[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        dest_parts.short_prop_[iprop]->insert( dest_parts.short_prop_[iprop]->begin() + dest_id, ( *short_prop_[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        dest_parts.uint64_prop_[iprop]->insert( dest_parts.uint64_prop_[iprop]->begin() + dest_id, ( *uint64_prop_[iprop] )[ipart] );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
//! Insert nPart particles starting at ipart to dest_id in dest_parts
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::copyParticles( unsigned int iPart, unsigned int nPart, Particles &dest_parts, int dest_id )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        dest_parts.double_prop_[iprop]->insert( dest_parts.double_prop_[iprop]->begin() + dest_id, double_prop_[iprop]->begin()+iPart, double_prop_[iprop]->begin()+iPart+nPart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        dest_parts.short_prop_[iprop]->insert( dest_parts.short_prop_[iprop]->begin() + dest_id, short_prop_[iprop]->begin()+iPart, short_prop_[iprop]->begin()+iPart+nPart );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        dest_parts.uint64_prop_[iprop]->insert( dest_parts.uint64_prop_[iprop]->begin() + dest_id, uint64_prop_[iprop]->begin()+iPart, uint64_prop_[iprop]->begin()+iPart+nPart );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Make a new particle at the position of another
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::makeParticleAt( Particles &source_particles, unsigned int ipart, double w, short q, double px, double py, double pz )
{
    for( unsigned int i=0 ; i<Position.size() ; i++ ) {
        Position[i].push_back( source_particles.Position[i][ipart] );
    }

    if( Position_old.size() > 0. ) {
        for( unsigned int i=0 ; i<Position_old.size() ; i++ ) {
            Position_old[i].push_back( source_particles.Position_old[i][ipart] );
        }
    }

    Momentum[0].push_back( px );
    Momentum[1].push_back( py );
    Momentum[2].push_back( pz );

    Weight.push_back( w );
    Charge.push_back( q );

    if( tracked ) {
        Id.push_back( 0 );
    }

    if( isQuantumParameter ) {
        Chi.push_back( 0. );
    }

    if( isMonteCarlo ) {
        Tau.push_back( 0. );
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! Suppress particle iPart
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::eraseParticle( unsigned int ipart, bool compute_cell_keys )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] ).erase( ( *double_prop_[iprop] ).begin()+ipart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] ).erase( ( *short_prop_[iprop] ).begin()+ipart );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] ).erase( ( *uint64_prop_[iprop] ).begin()+ipart );
    }

    if (compute_cell_keys) {
        cell_keys.erase(cell_keys.begin() + ipart);
    }

}

// ---------------------------------------------------------------------------------------------------------------------
//! Suppress all particles from iPart to the end of particle array
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::eraseParticleTrail( unsigned int ipart, bool compute_cell_keys )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] ).erase( ( *double_prop_[iprop] ).begin()+ipart, ( *double_prop_[iprop] ).end() );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] ).erase( ( *short_prop_[iprop] ).begin()+ipart, ( *short_prop_[iprop] ).end() );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] ).erase( ( *uint64_prop_[iprop] ).begin()+ipart, ( *uint64_prop_[iprop] ).end() );
    }

    if (compute_cell_keys) {
        cell_keys.erase( cell_keys.begin()+ipart, cell_keys.end() );
    }

}
// ---------------------------------------------------------------------------------------------------------------------
//! Suppress npart particles from ipart
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::eraseParticle( unsigned int ipart, unsigned int npart, bool compute_cell_keys )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] ).erase( ( *double_prop_[iprop] ).begin()+ipart, ( *double_prop_[iprop] ).begin()+ipart+npart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] ).erase( ( *short_prop_[iprop] ).begin()+ipart, ( *short_prop_[iprop] ).begin()+ipart+npart );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] ).erase( ( *uint64_prop_[iprop] ).begin()+ipart, ( *uint64_prop_[iprop] ).begin()+ipart+npart );
    }

    if (compute_cell_keys) {
        cell_keys.erase( cell_keys.begin()+ipart, cell_keys.begin()+ipart+npart );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
// Print parameters of particle iPart
// ---------------------------------------------------------------------------------------------------------------------
void Particles::print( unsigned int iPart )
{
    for( unsigned int i=0; i<Position.size(); i++ ) {
        cout << Position[i][iPart] << " ";
        //cout << Position_old[i][iPart] << " ";
    }
    for( unsigned int i=0; i<3; i++ ) {
        cout << Momentum[i][iPart] << " ";
    }
    cout << Weight[iPart] << " ";
    cout << Charge[iPart] << endl;;

    if( tracked ) {
        cout << Id[iPart] << endl;
    }

    if( isQuantumParameter ) {
        cout << Chi[iPart] << endl;
    }

    if( isMonteCarlo ) {
        cout << Tau[iPart] << endl;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Print parameters of particle iPart
// ---------------------------------------------------------------------------------------------------------------------
ostream &operator << ( ostream &out, const Particles &particles )
{
    for( unsigned int iPart=0; iPart<particles.Weight.size(); iPart++ ) {

        for( unsigned int i=0; i<particles.Position.size(); i++ ) {
            out << particles.Position[i][iPart] << " ";
            //out << particles.Position_old[i][iPart] << " ";
        }
        for( unsigned int i=0; i<3; i++ ) {
            out << particles.Momentum[i][iPart] << " ";
        }
        out << particles.Weight[iPart] << " ";
        out << particles.Charge[iPart] << endl;;

        if( particles.tracked ) {
            out << particles.Id[iPart] << endl;
        }

        if( particles.isQuantumParameter ) {
            out << particles.Chi[iPart] << endl;
        }

        if( particles.isMonteCarlo ) {
            out << particles.Tau[iPart] << endl;
        }
    }

    return ( out );
}


// ---------------------------------------------------------------------------------------------------------------------
// Exchange particles part1 & part2 memory location
// ---------------------------------------------------------------------------------------------------------------------
void Particles::swapParticle( unsigned int part1, unsigned int part2 )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        std::swap( ( *double_prop_[iprop] )[part1], ( *double_prop_[iprop] )[part2] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        std::swap( ( *short_prop_[iprop] )[part1], ( *short_prop_[iprop] )[part2] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        std::swap( ( *uint64_prop_[iprop] )[part1], ( *uint64_prop_[iprop] )[part2] );
    }
}


void Particles::swapParticle3( unsigned int part1, unsigned int part2, unsigned int part3 )
{
    // 1 ==> 2 ==> 3 ==> 1
    double temp;
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        temp = ( *double_prop_[iprop] )[part1];
        ( *double_prop_[iprop] )[part1] = ( *double_prop_[iprop] )[part3];
        ( *double_prop_[iprop] )[part3] = ( *double_prop_[iprop] )[part2];
        ( *double_prop_[iprop] )[part2] = temp;
    }

    short stemp;
    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        stemp = ( *short_prop_[iprop] )[part1];
        ( *short_prop_[iprop] )[part1] = ( *short_prop_[iprop] )[part3];
        ( *short_prop_[iprop] )[part3] = ( *short_prop_[iprop] )[part2];
        ( *short_prop_[iprop] )[part2] = stemp;
    }

    unsigned int uitemp;
    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        uitemp = ( *short_prop_[iprop] )[part1];
        ( *uint64_prop_[iprop] )[part1] = ( *uint64_prop_[iprop] )[part3];
        ( *uint64_prop_[iprop] )[part3] = ( *uint64_prop_[iprop] )[part2];
        ( *uint64_prop_[iprop] )[part2] = uitemp;
    }

}


void Particles::swapParticle4( unsigned int part1, unsigned int part2, unsigned int part3, unsigned int part4 )
{
    double temp;
    // 1 ==> 2 ==> 3 ==> 4 ==> 1
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        temp = ( *double_prop_[iprop] )[part1];
        ( *double_prop_[iprop] )[part1] = ( *double_prop_[iprop] )[part4];
        ( *double_prop_[iprop] )[part4] = ( *double_prop_[iprop] )[part3];
        ( *double_prop_[iprop] )[part3] = ( *double_prop_[iprop] )[part2];
        ( *double_prop_[iprop] )[part2] = temp;
    }

    short stemp;
    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        stemp = ( *short_prop_[iprop] )[part1];
        ( *short_prop_[iprop] )[part1] = ( *short_prop_[iprop] )[part4];
        ( *short_prop_[iprop] )[part4] = ( *short_prop_[iprop] )[part3];
        ( *short_prop_[iprop] )[part3] = ( *short_prop_[iprop] )[part2];
        ( *short_prop_[iprop] )[part2] = stemp;
    }

    unsigned int uitemp;
    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        uitemp = ( *short_prop_[iprop] )[part1];
        ( *uint64_prop_[iprop] )[part1] = ( *uint64_prop_[iprop] )[part4];
        ( *uint64_prop_[iprop] )[part4] = ( *uint64_prop_[iprop] )[part3];
        ( *uint64_prop_[iprop] )[part3] = ( *uint64_prop_[iprop] )[part2];
        ( *uint64_prop_[iprop] )[part2] = uitemp;
    }

}


void Particles::swapParticles( std::vector<unsigned int> parts )
{
    // parts[0] ==> parts[1] ==> parts[2] ==> parts[parts.size()-1] ==> parts[0]

    copyParticle( parts.back() );
    translateParticles( parts );
    overwriteParticle( size()-1, parts[0] );
    eraseParticle( size()-1 );

}


void Particles::translateParticles( std::vector<unsigned int> parts )
{
    // parts[0] ==> parts[1] ==> parts[2] ==> parts[parts.size()-1]

    for( int icycle = parts.size()-2; icycle >=0; icycle-- ) {
        overwriteParticle( parts[icycle], parts[icycle+1] );
    }

}


// ---------------------------------------------------------------------------------------------------------------------
//! Move particle src_particle into dest_particle memory location, erasing dest_particle.
//! Warning: do not update first_index and last_index
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwriteParticle( unsigned int src_particle, unsigned int dest_particle, bool compute_cell_keys )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] )[dest_particle] = ( *double_prop_[iprop] )[src_particle];
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] )[dest_particle] = ( *short_prop_[iprop] )[src_particle];
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] )[dest_particle] = ( *uint64_prop_[iprop] )[src_particle];
    }

    if (compute_cell_keys) {
        cell_keys[dest_particle] = cell_keys[src_particle];
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! Move particle part1->part1+N into part2->part2+N memory location erasing part2->part2+N.
//! Warning: do not update first_index and last_index
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwriteParticle( unsigned int part1,
                                   unsigned int part2,
                                   unsigned int N,
                                   bool compute_cell_keys)
{
    unsigned int sizepart = N*sizeof( Position[0][0] );
    unsigned int sizecharge = N*sizeof( Charge[0] );
    unsigned int sizeid = N*sizeof( Id[0] );

    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        memcpy( & ( *double_prop_[iprop] )[part2],  &( *double_prop_[iprop] )[part1], sizepart );
        // std::copy( *double_prop_[iprop]->begin()+part1,  *double_prop_[iprop]->begin() + part1 + N, *double_prop_[iprop]->begin() + part2 );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        memcpy( & ( *short_prop_[iprop] )[part2],  &( *short_prop_[iprop] )[part1], sizecharge );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        memcpy( & ( *uint64_prop_[iprop] )[part2],  &( *uint64_prop_[iprop] )[part1], sizeid );
    }

    if (compute_cell_keys) {
        //std::copy( cell_keys.begin()+part1,  cell_keys.begin() + part1 + N, cell_keys.begin() + part2 );
        memcpy( &cell_keys[part2],  &cell_keys[part1], N*sizeof( cell_keys[0] ) );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Move particle part1 into part2 memory location of dest vector, erasing part2.
//! Warning: do not update first_index and last_index
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwriteParticle( unsigned int part1, Particles &dest_parts, unsigned int part2 )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *dest_parts.double_prop_[iprop] )[part2] = ( *double_prop_[iprop] )[part1];
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *dest_parts.short_prop_[iprop] )[part2] = ( *short_prop_[iprop] )[part1];
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *dest_parts.uint64_prop_[iprop] )[part2] = ( *uint64_prop_[iprop] )[part1];
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location of dest vector, erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwriteParticle( unsigned int part1, Particles &dest_parts, unsigned int part2, unsigned int N )
{
    unsigned int sizepart = N*sizeof( Position[0][0] );
    unsigned int sizecharge = N*sizeof( Charge[0] );
    unsigned int sizeid = N*sizeof( Id[0] );

    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        memcpy( & ( *dest_parts.double_prop_[iprop] )[part2],  &( *double_prop_[iprop] )[part1], sizepart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        memcpy( & ( *dest_parts.short_prop_[iprop] )[part2],  &( *short_prop_[iprop] )[part1], sizecharge );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        memcpy( & ( *dest_parts.uint64_prop_[iprop] )[part2],  &( *uint64_prop_[iprop] )[part1], sizeid );
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Exchange N particles part1->part1+N & part2->part2+N memory location
// ---------------------------------------------------------------------------------------------------------------------
void Particles::swapParticle( unsigned int part1, unsigned int part2, unsigned int N )
{
    double *buffer[N];

    unsigned int sizepart = N*sizeof( Position[0][0] );
    unsigned int sizecharge = N*sizeof( Charge[0] );
    unsigned int sizeid = N*sizeof( Id[0] );

    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        memcpy( buffer, &( ( *double_prop_[iprop] )[part1] ), sizepart );
        memcpy( &( ( *double_prop_[iprop] )[part1] ), &( ( *double_prop_[iprop] )[part2] ), sizepart );
        memcpy( &( ( *double_prop_[iprop] )[part2] ), buffer, sizepart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        memcpy( buffer, &( ( *short_prop_[iprop] )[part1] ), sizecharge );
        memcpy( &( ( *short_prop_[iprop] )[part1] ), &( ( *short_prop_[iprop] )[part2] ), sizecharge );
        memcpy( &( ( *short_prop_[iprop] )[part2] ), buffer, sizecharge );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        memcpy( buffer, &( ( *uint64_prop_[iprop] )[part1] ), sizeid );
        memcpy( &( ( *uint64_prop_[iprop] )[part1] ), &( ( *uint64_prop_[iprop] )[part2] ), sizeid );
        memcpy( &( ( *uint64_prop_[iprop] )[part2] ), buffer, sizeid );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
// Move iPart at the end of vectors (to do for MPI)
// ---------------------------------------------------------------------------------------------------------------------
void Particles::pushToEnd( unsigned int iPart )
{

}

// ---------------------------------------------------------------------------------------------------------------------
// Create a new particle at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::createParticle()
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] ).push_back( 0. );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] ).push_back( 0 );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] ).push_back( 0 );
    }
//MESSAGE("create1");
}

// ---------------------------------------------------------------------------------------------------------------------
//! Create n_additional_particles new particles at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::createParticles( int n_additional_particles )
{
    int nParticles = size();
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] ).resize( nParticles+n_additional_particles, 0. );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] ).resize( nParticles+n_additional_particles, 0 );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] ).resize( nParticles+n_additional_particles, 0 );
    }

    cell_keys.resize( nParticles+n_additional_particles, 0);

//MESSAGE("create2");
}

// ---------------------------------------------------------------------------------------------------------------------
//! Create n_additional_particles new particles at the position pstart
// ---------------------------------------------------------------------------------------------------------------------
void Particles::createParticles( int n_additional_particles, int pstart )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] ).insert( ( *double_prop_[iprop] ).begin()+pstart, n_additional_particles, 0. );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] ).insert( ( *short_prop_[iprop] ).begin()+pstart, n_additional_particles, 0 );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] ).insert( ( *uint64_prop_[iprop] ).begin()+pstart, n_additional_particles, 0 );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! This method erases some particles of the particles vector using a mask vector between istart and iend.
//! This function is optimized.
//! The mask determines which particles to keep (>= 0) and which to delete (< 0)
//! Particle order is kept in case the vector is sorted
//! Warning: This method do not update count, first_index and last_index in Species
// ---------------------------------------------------------------------------------------------------------------------
void Particles::eraseParticlesWithMask( int istart, int iend, vector <int> & mask ) {

    unsigned int idest = (unsigned int) istart;
    unsigned int isrc = (unsigned int) istart;
    while (isrc < (unsigned int) iend) {
        if (mask[idest] < 0) {
            if (mask[isrc] >= 0) {
                overwriteParticle( isrc, idest);
                cell_keys[idest] = cell_keys[isrc];
                mask[idest] = 1;
                mask[isrc] = -1;
                idest++;
            } else {
                isrc++;
            }
        } else {
            idest++;
            isrc = idest;
        }
    }

    // At the end we resize particles
    resize(idest);
    //cell_keys.resize(idest);
}

// ---------------------------------------------------------------------------------------------------------------------
//! This method erases some particles of the particles vector using cell_keys as a mask vector between istart and iend.
//! This function is optimized.
//! Warning: This method do not update count, first_index and last_index in Species
// ---------------------------------------------------------------------------------------------------------------------
void Particles::eraseParticlesWithMask( int istart, int iend) {

    unsigned int idest = (unsigned int) istart;
    unsigned int isrc = (unsigned int) istart;
    while (isrc < (unsigned int) iend) {
        if (cell_keys[idest] < 0) {
            if (cell_keys[isrc] >= 0) {
                overwriteParticle( isrc, idest);
                cell_keys[idest] = cell_keys[isrc];
                idest++;
            } else {
                isrc++;
            }
        } else {
            idest++;
            isrc = idest;
        }
    }

    // At the end we resize particles
    resize(idest);
    //cell_keys.resize(idest);
}

// ---------------------------------------------------------------------------------------------------------------------
//! This method erases some particles of the particles vector using a mask vector.
//! This function is not optimized.
//! If the value of the mask is negative, then the particle is deleted.
//! Warning: This method do not update count, first_index and last_index in Species
// ---------------------------------------------------------------------------------------------------------------------
// void Particles::eraseParticlesWithMask( int istart, int iend, vector <int> & mask) {
//
//     //unsigned int deleted_particles = 0;
//
//     for (int ip = istart ; ip < iend ; ip++) {
//         if (mask[ip] < 0) {
//             eraseParticle(ip);
//             cell_keys.erase(cell_keys.begin()+ip);
//             //deleted_particles += 1;
//         }
//     }
// }

// ---------------------------------------------------------------------------------------------------------------------
//! This method erases particles according to the provided mask
//! between istart and iend
// ---------------------------------------------------------------------------------------------------------------------
// void Particles::eraseParticlesWithMask( int istart, int iend, vector <bool> & to_be_erased) {
//     int last_index = iend;
//     for (int ip = iend ; ip >= istart ; ip-- )
//     {
//         if (to_be_erased[ip]) {
//             if (last_index != ip) {
//                 overwriteParticle( last_index, ip);
//             }
//             last_index--;
//         }
//     }
// }

// ---------------------------------------------------------------------------------------------------------------------
//! This method eliminates the space between the bins
//! (presence of empty particles between the bins)
// ---------------------------------------------------------------------------------------------------------------------
void Particles::compress(bool compute_cell_keys) {

    unsigned int nbin = numberOfBins();

    for (int ibin = 0 ; ibin < nbin-1 ; ibin++) {

        // Removal of the photons
        const unsigned int nb_deleted_photon = first_index[ibin+1] - last_index[ibin];

        if( nb_deleted_photon > 0 ) {
            eraseParticle( last_index[ibin], nb_deleted_photon, compute_cell_keys );

            for( int ii=ibin+1; ii<nbin; ii++ ) {
                first_index[ii] -= nb_deleted_photon;
                last_index[ii] -= nb_deleted_photon;
            }
        }
    }
    eraseParticleTrail( last_index[nbin-1], true );
}

// ---------------------------------------------------------------------------------------------------------------------
//! This method eliminates the space between the bins
//! (presence of empty particles between the bins)
// ---------------------------------------------------------------------------------------------------------------------
// void Particles::compress() {
//     for (auto ibin = 1 ; ibin < first_index.size() ; ibin++) {
//
//         // Compute the number of particles
//         unsigned int particle_number = last_index[ibin] - first_index[ibin];
//
//         // Compute the space between the bins
//         unsigned int bin_space = first_index[ibin] - last_index[ibin-1];
//
//         // Determine first index and number of particles to copy.
//         // We copy from first index to the end to limit the number of copy (more efficient than copying the full bin to keep the same order)
//
//         // Compute the number of particles
//         unsigned int copy_particle_number = 0;
//
//         if (bin_space > 0) {
//
//             // if last_index[ibin] - bin_space < first_index[ibin], it means that the empty space is larger than the number of particles in ibin
//             // then we move the full bin
//             // Else we only move the particles from copy_first_bin to last_index[ibin]
//             unsigned int copy_first_index = last_index[ibin] - bin_space;
//
//             if (copy_first_index < first_index[ibin]) {
//                 copy_first_index = first_index[ibin];
//                 copy_particle_number = particle_number;
//             } else {
//                 copy_particle_number = bin_space;
//             }
//
//             if (copy_particle_number>0) {
//                 overwriteParticle(copy_first_index, last_index[ibin-1], copy_particle_number, true );
//             }
//
//             //Update bin indexes
//             first_index[ibin] = last_index[ibin-1];
//             last_index[ibin] = first_index[ibin] + copy_particle_number;
//
//         }
//     } // for bin
//
//     unsigned int particles_to_erase = Weight.size() - last_index[first_index.size()-1];
//
//     if (particles_to_erase > 0) {
//
//         // Particles is rezised to fit the real number of particles
//         // resize(last_index[first_index.size()]-1);
//         eraseParticleTrail( last_index[first_index.size()-1], true );
//     }
// }

void Particles::sum(int ibin_min, int ibin_max) {

    double sum_px = 0;
    double sum_py = 0;
    double sum_mx = 0;
    double sum_my = 0;
    int iterations = 0;
    int nb_particles_total = 0;

    for (int ibin = ibin_min ; ibin < ibin_max ; ibin++) {
        for (int ipart = first_index[ibin] ; ipart < last_index[ibin] ; ipart++) {
            if (Weight[ipart] > 0) {
                sum_px += Position[0][ipart];
                sum_py += Position[1][ipart];
                sum_mx += Momentum[0][ipart];
                sum_my += Momentum[1][ipart];
                iterations += 1;
            }
            nb_particles_total ++;
        }
    }
    std::cerr << " Position_x: " << sum_px
              << " Position_y: " << sum_py
              << " Momentum_x: " << sum_mx
              << " Momentum_y: " << sum_my
              << " nb real particles: " << iterations
              << " nb particles: " << nb_particles_total
              << std::endl;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Move ipart at new_pos in the particles data structure
// ---------------------------------------------------------------------------------------------------------------------
void Particles::moveParticles( int iPart, int new_pos )
{
    for( unsigned int iprop=0 ; iprop<double_prop_.size() ; iprop++ ) {
        ( *double_prop_[iprop] ).insert( ( *double_prop_[iprop] ).begin()+new_pos,( *double_prop_[iprop] )[iPart]  );
    }

    for( unsigned int iprop=0 ; iprop<short_prop_.size() ; iprop++ ) {
        ( *short_prop_[iprop] ).insert( ( *short_prop_[iprop] ).begin()+new_pos, ( *short_prop_[iprop] )[iPart] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop_.size() ; iprop++ ) {
        ( *uint64_prop_[iprop] ).insert( ( *uint64_prop_[iprop] ).begin()+new_pos,( *uint64_prop_[iprop] )[iPart]  );
    }

    eraseParticle( iPart+1 );
}

// ---------------------------------------------------------------------------------------------------------------------
// Create nParticles new particles at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
//void Particles::createParticles(int nAdditionalParticles )
//{
//    int nParticles = size();
//    for (unsigned int i=0; i<Position.size(); i++) {
//        Position[i].resize(nParticles+nAdditionalParticles,0.);
//        Position_old[i].resize(nParticles+nAdditionalParticles,0.);
//    }
//
//    for (unsigned int i=0; i<3; i++) {
//        Momentum[i].resize(nParticles+nAdditionalParticles,0.);
//    }
//    Weight.resize(nParticles+nAdditionalParticles,0.);
//    Charge.resize(nParticles+nAdditionalParticles,0);
//
//    if (tracked)
//        Id.resize(nParticles+nAdditionalParticles,0);
//
//    if (isQuantumParameter)
//        Chi.resize(nParticles+nAdditionalParticles,0.);
//
//}

// ---------------------------------------------------------------------------------------------------------------------
// Test if ipart is in the local patch
//---------------------------------------------------------------------------------------------------------------------
bool Particles::isParticleInDomain( unsigned int ipart, Patch *patch )
{
    for( unsigned int i=0; i<Position.size(); i++ ) {
        if( Position[i][ipart] <  patch->getDomainLocalMin( i )
         || Position[i][ipart] >= patch->getDomainLocalMax( i ) ) {
            return false;
        }
    }
    return true;
}


void Particles::sortById()
{
    if( !tracked ) {
        ERROR( "Impossible" );
        return;
    }
    int nParticles( Weight.size() );

    bool stop;
    int jPart( 0 );
    do {
        stop = true;
        for( int iPart = nParticles-1 ; iPart > jPart ; --iPart ) {
            if( Id[iPart] < Id[iPart-1] ) {
                swapParticle( iPart, jPart );
                stop = false;
            }
        }
        jPart++;
    } while( !stop );

}

void Particles::initializeDataOnDevice()
{
    ERROR( "Device only feature, should not have come here!" );
}
void Particles::syncGPU()
{
    ERROR( "Device only feature, should not have come here!" );
}
void Particles::syncCPU()
{
    ERROR( "Device only feature, should not have come here!" );
}

void Particles::extractParticles( Particles* particles_to_move )
{
    particles_to_move->clear();
    for ( int ipart=0 ; ipart<size() ; ipart++ ) {
        if ( cell_keys[ipart] == -1 ) {
            copyParticle( ipart, *particles_to_move );
        }
    }
}

void Particles::savePositions() {
    unsigned int ndim = Position.size(), npart = size();
    double *p[3], *pold[3];
    for( unsigned int i = 0 ; i<ndim ; i++ ) {
        p[i] =  &( Position[i][0] );
        pold[i] =  &( Position_old[i][0] );
    }
    if (ndim == 1) {
        #pragma omp simd
        for( unsigned int ipart=0 ; ipart<npart; ipart++ ) {
            pold[0][ipart] = p[0][ipart];
        }
    } else if (ndim == 2) {
        #pragma omp simd
        for( unsigned int ipart=0 ; ipart<npart; ipart++ ) {
            pold[0][ipart] = p[0][ipart];
            pold[1][ipart] = p[1][ipart];
        }
    } else if (ndim == 3) {
        #pragma omp simd
        for( unsigned int ipart=0 ; ipart<npart; ipart++ ) {
            pold[0][ipart] = p[0][ipart];
            pold[1][ipart] = p[1][ipart];
            pold[2][ipart] = p[2][ipart];
        }
    }

    // --- old version that only vectorizes with Intel ---
    // #pragma omp simd
    // for( unsigned int ipart=0 ; ipart<npart; ipart++ ) {
    //     for( unsigned int i = 0 ; i<ndim ; i++ ) {
    //         pold[i][ipart] = p[i][ipart];
    //     }
    // }
    // -----------------------------------------------------

}

int Particles::eraseLeavingParticles()
{
    ERROR( "Device only feature, should not have come here!" );
    return 0;
}

int Particles::injectParticles( Particles *particles_to_inject )
{
    ERROR( "Device only feature, should not have come here! On CPU it's done in sortParticles." );
    return 0;
}

void Particles::importAndSortParticles( Particles *particles_to_inject )
{
    ERROR( "Device only feature, should not have come here! On CPU it's done in sortParticles." );
}

unsigned int Particles::deviceCapacity() const
{
    ERROR( "deviceCapacity is a feature only available for accelerator device" );
    return 0;
}

#ifdef __DEBUG
bool Particles::testMove( int iPartStart, int iPartEnd, Params &params )
{
    for( int iDim = 0 ; iDim < Position.size() ; iDim++ ) {
        double dx2 = params.cell_length[iDim];//*params.cell_length[iDim];
        for( int iPart = iPartStart ; iPart < iPartEnd ; iPart++ ) {
            if( dist( iPart, iDim ) > dx2 ) {
                ERROR( "Too large displacment for particle : " << iPart << "\t: " << ( *this )( iPart ) );
                return false;
            }
        }
    }
    return true;

}
#endif

Particle Particles::operator()( unsigned int iPart )
{
    return  Particle( *this, iPart );
}
