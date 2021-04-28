#include "Particles.h"

#include <cstring>
#include <iostream>

#include "Params.h"
#include "Patch.h"
#include "Species.h"

#include "Particle.h"

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

    double_prop.resize( 0 );
    short_prop.resize( 0 );
    uint64_prop.resize( 0 );
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
    //if (nParticles > Weight.capacity()) {
    //    WARNING("You should increase c_part_max in specie namelist");
    //}
    if( Weight.size()==0 ) {
        float c_part_max =1.2;
        //float c_part_max = part.c_part_max;
        //float c_part_max = params.species_param[0].c_part_max;
        reserve( round( c_part_max * nParticles ), nDim );
    }

    resize( nParticles, nDim, keep_position_old );
    //cell_keys.resize( nParticles );

    if( double_prop.empty() ) {  // do this just once

        Position.resize( nDim );
        for( unsigned int i=0 ; i< nDim ; i++ ) {
            double_prop.push_back( &( Position[i] ) );
        }

        for( unsigned int i=0 ; i< 3 ; i++ ) {
            double_prop.push_back( &( Momentum[i] ) );
        }

        double_prop.push_back( &Weight );

        if( keep_position_old ) {
            Position_old.resize( nDim );
            for( unsigned int i=0 ; i< nDim ; i++ ) {
                double_prop.push_back( &( Position_old[i] ) );
            }
        }

        short_prop.push_back( &Charge );
        if( tracked ) {
            uint64_prop.push_back( &Id );
        }

        // Quantum parameter (for QED effects):
        // - if radiation reaction (continuous or discontinuous)
        // - if multiphoton-Breit-Wheeler if photons
        if( isQuantumParameter ) {
            double_prop.push_back( &Chi );
        }

        // Optical Depth for Monte-Carlo processes:
        // - if the discontinuous (Monte-Carlo) radiation reaction
        // are activated, tau is the incremental optical depth to emission
        if( isMonteCarlo ) {
            double_prop.push_back( &Tau );
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
// Set capacity of Particles vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::reserve( unsigned int n_part_max, unsigned int nDim )
{
    return;

    Position.resize( nDim );
    Position_old.resize( nDim );
    for( unsigned int i=0 ; i< nDim ; i++ ) {
        Position[i].reserve( n_part_max );
        Position_old[i].reserve( n_part_max );
    }
    Momentum.resize( 3 );
    for( unsigned int i=0 ; i< 3 ; i++ ) {
        Momentum[i].reserve( n_part_max );
    }
    Weight.reserve( n_part_max );
    Charge.reserve( n_part_max );

    if( tracked ) {
        Id.reserve( n_part_max );
    }

    if( isQuantumParameter ) {
        Chi.reserve( n_part_max );
    }

    if( isMonteCarlo ) {
        Tau.reserve( n_part_max );
    }

    cell_keys.reserve( n_part_max );

}

void Particles::initializeReserve( unsigned int npart_max, Particles &part )
{
    initialize( 0, part );
    reserve( npart_max, part.dimension() );
}



// ---------------------------------------------------------------------------------------------------------------------
//Resize Particle vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::resize( unsigned int nParticles, unsigned int nDim, bool keep_position_old )
{
    Position.resize( nDim );
    for( unsigned int i=0 ; i<nDim ; i++ ) {
        Position[i].resize( nParticles, 0. );
    }
    
    if( keep_position_old ) {
        Position_old.resize( nDim );
        for( unsigned int i=0 ; i<nDim ; i++ ) {
            Position_old[i].resize( nParticles, 0. );
        }
    }

    Momentum.resize( 3 );
    for( unsigned int i=0 ; i< 3 ; i++ ) {
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

    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] ).resize( nParticles, 0. );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] ).resize( nParticles, 0 );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] ).resize( nParticles, 0 );
    }

    cell_keys.resize( nParticles, 0. );

}

// ---------------------------------------------------------------------------------------------------------------------
//! Resize the cell_keys vector
// ---------------------------------------------------------------------------------------------------------------------
void Particles::resizeCellKeys(unsigned int nParticles)
{
    cell_keys.resize( nParticles, 0. );
}

// ---------------------------------------------------------------------------------------------------------------------
// Remove extra capacity of Particles vectors
// Cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::shrinkToFit()
{

    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        std::vector<double>( *double_prop[iprop] ).swap( *double_prop[iprop] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        std::vector<short>( *short_prop[iprop] ).swap( *short_prop[iprop] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        std::vector<uint64_t>( *uint64_prop[iprop] ).swap( *uint64_prop[iprop] );
    }
    
    //cell_keys.swap(cell_keys);
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Reset of Particles vectors
// Cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::clear()
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        double_prop[iprop]->clear();
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        short_prop[iprop]->clear();
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        uint64_prop[iprop]->clear();
    }
    
    //cell_keys.clear();
    
}

//! copy particles ipart at the end of the particle vector
//! cell keys not affected
void Particles::copyParticle( unsigned int ipart )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        double_prop[iprop]->push_back( ( *double_prop[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        short_prop[iprop]->push_back( ( *short_prop[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        uint64_prop[iprop]->push_back( ( *uint64_prop[iprop] )[ipart] );
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! Copy particle iPart at the end of dest_parts
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::copyParticle( unsigned int ipart, Particles &dest_parts )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        dest_parts.double_prop[iprop]->push_back( ( *double_prop[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        dest_parts.short_prop[iprop]->push_back( ( *short_prop[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        dest_parts.uint64_prop[iprop]->push_back( ( *uint64_prop[iprop] )[ipart] );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Insert particle iPart at dest_id in dest_parts
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::copyParticle( unsigned int ipart, Particles &dest_parts, int dest_id )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        dest_parts.double_prop[iprop]->insert( dest_parts.double_prop[iprop]->begin() + dest_id, ( *double_prop[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        dest_parts.short_prop[iprop]->insert( dest_parts.short_prop[iprop]->begin() + dest_id, ( *short_prop[iprop] )[ipart] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        dest_parts.uint64_prop[iprop]->insert( dest_parts.uint64_prop[iprop]->begin() + dest_id, ( *uint64_prop[iprop] )[ipart] );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
//! Insert nPart particles starting at ipart to dest_id in dest_parts
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::copyParticles( unsigned int iPart, unsigned int nPart, Particles &dest_parts, int dest_id )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        dest_parts.double_prop[iprop]->insert( dest_parts.double_prop[iprop]->begin() + dest_id, double_prop[iprop]->begin()+iPart, double_prop[iprop]->begin()+iPart+nPart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        dest_parts.short_prop[iprop]->insert( dest_parts.short_prop[iprop]->begin() + dest_id, short_prop[iprop]->begin()+iPart, short_prop[iprop]->begin()+iPart+nPart );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        dest_parts.uint64_prop[iprop]->insert( dest_parts.uint64_prop[iprop]->begin() + dest_id, uint64_prop[iprop]->begin()+iPart, uint64_prop[iprop]->begin()+iPart+nPart );
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! Copy particle iPart at the end of dest_parts -- safe
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::copyParticleSafe( unsigned int ipart, Particles &dest_parts )
{
    unsigned int ndouble = double_prop.size();
    if( dest_parts.double_prop.size() < ndouble ) {
        ndouble = dest_parts.double_prop.size();
    }
    unsigned int nuint = uint64_prop.size();
    if( dest_parts.uint64_prop.size() < nuint ) {
        nuint = dest_parts.uint64_prop.size();
    }
    
    for( unsigned int iprop=0 ; iprop<ndouble ; iprop++ ) {
        dest_parts.double_prop[iprop]->push_back( ( *double_prop[iprop] )[ipart] );
    }
    
    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        dest_parts.short_prop[iprop]->push_back( ( *short_prop[iprop] )[ipart] );
    }
    
    for( unsigned int iprop=0 ; iprop<nuint ; iprop++ ) {
        dest_parts.uint64_prop[iprop]->push_back( ( *uint64_prop[iprop] )[ipart] );
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! Suppress particle iPart
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::eraseParticle( unsigned int ipart )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] ).erase( ( *double_prop[iprop] ).begin()+ipart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] ).erase( ( *short_prop[iprop] ).begin()+ipart );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] ).erase( ( *uint64_prop[iprop] ).begin()+ipart );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
//! Suppress all particles from iPart to the end of particle array
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::eraseParticleTrail( unsigned int ipart )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] ).erase( ( *double_prop[iprop] ).begin()+ipart, ( *double_prop[iprop] ).end() );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] ).erase( ( *short_prop[iprop] ).begin()+ipart, ( *short_prop[iprop] ).end() );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] ).erase( ( *uint64_prop[iprop] ).begin()+ipart, ( *uint64_prop[iprop] ).end() );
    }

}
// ---------------------------------------------------------------------------------------------------------------------
//! Suppress npart particles from ipart
//! cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void Particles::eraseParticle( unsigned int ipart, unsigned int npart )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] ).erase( ( *double_prop[iprop] ).begin()+ipart, ( *double_prop[iprop] ).begin()+ipart+npart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] ).erase( ( *short_prop[iprop] ).begin()+ipart, ( *short_prop[iprop] ).begin()+ipart+npart );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] ).erase( ( *uint64_prop[iprop] ).begin()+ipart, ( *uint64_prop[iprop] ).begin()+ipart+npart );
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
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        std::swap( ( *double_prop[iprop] )[part1], ( *double_prop[iprop] )[part2] );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        std::swap( ( *short_prop[iprop] )[part1], ( *short_prop[iprop] )[part2] );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        std::swap( ( *uint64_prop[iprop] )[part1], ( *uint64_prop[iprop] )[part2] );
    }
}


void Particles::swapParticle3( unsigned int part1, unsigned int part2, unsigned int part3 )
{
    // 1 ==> 2 ==> 3 ==> 1
    double temp;
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        temp = ( *double_prop[iprop] )[part1];
        ( *double_prop[iprop] )[part1] = ( *double_prop[iprop] )[part3];
        ( *double_prop[iprop] )[part3] = ( *double_prop[iprop] )[part2];
        ( *double_prop[iprop] )[part2] = temp;
    }

    short stemp;
    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        stemp = ( *short_prop[iprop] )[part1];
        ( *short_prop[iprop] )[part1] = ( *short_prop[iprop] )[part3];
        ( *short_prop[iprop] )[part3] = ( *short_prop[iprop] )[part2];
        ( *short_prop[iprop] )[part2] = stemp;
    }

    unsigned int uitemp;
    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        uitemp = ( *short_prop[iprop] )[part1];
        ( *uint64_prop[iprop] )[part1] = ( *uint64_prop[iprop] )[part3];
        ( *uint64_prop[iprop] )[part3] = ( *uint64_prop[iprop] )[part2];
        ( *uint64_prop[iprop] )[part2] = uitemp;
    }

}


void Particles::swapParticle4( unsigned int part1, unsigned int part2, unsigned int part3, unsigned int part4 )
{
    double temp;
    // 1 ==> 2 ==> 3 ==> 4 ==> 1
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        temp = ( *double_prop[iprop] )[part1];
        ( *double_prop[iprop] )[part1] = ( *double_prop[iprop] )[part4];
        ( *double_prop[iprop] )[part4] = ( *double_prop[iprop] )[part3];
        ( *double_prop[iprop] )[part3] = ( *double_prop[iprop] )[part2];
        ( *double_prop[iprop] )[part2] = temp;
    }

    short stemp;
    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        stemp = ( *short_prop[iprop] )[part1];
        ( *short_prop[iprop] )[part1] = ( *short_prop[iprop] )[part4];
        ( *short_prop[iprop] )[part4] = ( *short_prop[iprop] )[part3];
        ( *short_prop[iprop] )[part3] = ( *short_prop[iprop] )[part2];
        ( *short_prop[iprop] )[part2] = stemp;
    }

    unsigned int uitemp;
    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        uitemp = ( *short_prop[iprop] )[part1];
        ( *uint64_prop[iprop] )[part1] = ( *uint64_prop[iprop] )[part4];
        ( *uint64_prop[iprop] )[part4] = ( *uint64_prop[iprop] )[part3];
        ( *uint64_prop[iprop] )[part3] = ( *uint64_prop[iprop] )[part2];
        ( *uint64_prop[iprop] )[part2] = uitemp;
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
// Move particle src_particle into dest_particle memory location, erasing dest_particle.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwriteParticle( unsigned int src_particle, unsigned int dest_particle )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] )[dest_particle] = ( *double_prop[iprop] )[src_particle];
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] )[dest_particle] = ( *short_prop[iprop] )[src_particle];
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] )[dest_particle] = ( *uint64_prop[iprop] )[src_particle];
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwriteParticle( unsigned int part1, unsigned int part2, unsigned int N )
{
    unsigned int sizepart = N*sizeof( Position[0][0] );
    unsigned int sizecharge = N*sizeof( Charge[0] );
    unsigned int sizeid = N*sizeof( Id[0] );

    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        memcpy( & ( *double_prop[iprop] )[part2],  &( *double_prop[iprop] )[part1], sizepart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        memcpy( & ( *short_prop[iprop] )[part2],  &( *short_prop[iprop] )[part1], sizecharge );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        memcpy( & ( *uint64_prop[iprop] )[part2],  &( *uint64_prop[iprop] )[part1], sizeid );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1 into part2 memory location of dest vector, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwriteParticle( unsigned int part1, Particles &dest_parts, unsigned int part2 )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *dest_parts.double_prop[iprop] )[part2] = ( *double_prop[iprop] )[part1];
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *dest_parts.short_prop[iprop] )[part2] = ( *short_prop[iprop] )[part1];
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *dest_parts.uint64_prop[iprop] )[part2] = ( *uint64_prop[iprop] )[part1];
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

    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        memcpy( & ( *dest_parts.double_prop[iprop] )[part2],  &( *double_prop[iprop] )[part1], sizepart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        memcpy( & ( *dest_parts.short_prop[iprop] )[part2],  &( *short_prop[iprop] )[part1], sizecharge );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        memcpy( & ( *dest_parts.uint64_prop[iprop] )[part2],  &( *uint64_prop[iprop] )[part1], sizeid );
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

    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        memcpy( buffer, &( ( *double_prop[iprop] )[part1] ), sizepart );
        memcpy( &( ( *double_prop[iprop] )[part1] ), &( ( *double_prop[iprop] )[part2] ), sizepart );
        memcpy( &( ( *double_prop[iprop] )[part2] ), buffer, sizepart );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        memcpy( buffer, &( ( *short_prop[iprop] )[part1] ), sizecharge );
        memcpy( &( ( *short_prop[iprop] )[part1] ), &( ( *short_prop[iprop] )[part2] ), sizecharge );
        memcpy( &( ( *short_prop[iprop] )[part2] ), buffer, sizecharge );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        memcpy( buffer, &( ( *uint64_prop[iprop] )[part1] ), sizeid );
        memcpy( &( ( *uint64_prop[iprop] )[part1] ), &( ( *uint64_prop[iprop] )[part2] ), sizeid );
        memcpy( &( ( *uint64_prop[iprop] )[part2] ), buffer, sizeid );
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
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] ).push_back( 0. );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] ).push_back( 0 );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] ).push_back( 0 );
    }
//MESSAGE("create1");
}

// ---------------------------------------------------------------------------------------------------------------------
// Create nParticles new particles at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::createParticles( int nAdditionalParticles )
{
    int nParticles = size();
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] ).resize( nParticles+nAdditionalParticles, 0. );
    }

    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] ).resize( nParticles+nAdditionalParticles, 0 );
    }

    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] ).resize( nParticles+nAdditionalParticles, 0 );
    }
    
    cell_keys.resize( nParticles+nAdditionalParticles, 0);

//MESSAGE("create2");
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
// Create nParticles new particles at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::createParticles( int nAdditionalParticles, int pstart )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] ).insert( ( *double_prop[iprop] ).begin()+pstart, nAdditionalParticles, 0. );
    }
    
    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] ).insert( ( *short_prop[iprop] ).begin()+pstart, nAdditionalParticles, 0 );
    }
    
    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] ).insert( ( *uint64_prop[iprop] ).begin()+pstart, nAdditionalParticles, 0 );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Move ipart at new_pos in the particles data structure
// ---------------------------------------------------------------------------------------------------------------------
void Particles::moveParticles( int iPart, int new_pos )
{
    for( unsigned int iprop=0 ; iprop<double_prop.size() ; iprop++ ) {
        ( *double_prop[iprop] ).insert( ( *double_prop[iprop] ).begin()+new_pos,( *double_prop[iprop] )[iPart]  );
    }
    
    for( unsigned int iprop=0 ; iprop<short_prop.size() ; iprop++ ) {
        ( *short_prop[iprop] ).insert( ( *short_prop[iprop] ).begin()+new_pos, ( *short_prop[iprop] )[iPart] );
    }
    
    for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
        ( *uint64_prop[iprop] ).insert( ( *uint64_prop[iprop] ).begin()+new_pos,( *uint64_prop[iprop] )[iPart]  );
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

void Particles::savePositions() {
    unsigned int ndim = Position.size(), npart = size();
    double *p[3], *pold[3];
    for( unsigned int i = 0 ; i<ndim ; i++ ) {
        p[i] =  &( Position[i][0] );
        pold[i] =  &( Position_old[i][0] );
    }
    #pragma omp simd
    for( unsigned int ipart=0 ; ipart<npart; ipart++ ) {
        for( unsigned int i = 0 ; i<ndim ; i++ ) {
            pold[i][ipart] = p[i][ipart];
        }
    }
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
