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
    is_test = false;
    isQuantumParameter = false;
    isMonteCarlo = false;

    double_prop.resize( 0 );
    short_prop.resize( 0 );
    uint64_prop.resize( 0 );
}

// ---------------------------------------------------------------------------------------------------------------------
// Create nParticles null particles of nDim size
// ---------------------------------------------------------------------------------------------------------------------
void Particles::initialize( unsigned int nParticles, unsigned int nDim )
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

    resize( nParticles, nDim );
    cell_keys.resize( nParticles );

    if( double_prop.empty() ) {  // do this just once

        Position.resize( nDim );
        for( unsigned int i=0 ; i< nDim ; i++ ) {
            double_prop.push_back( &( Position[i] ) );
        }

        for( unsigned int i=0 ; i< 3 ; i++ ) {
            double_prop.push_back( &( Momentum[i] ) );
        }

        double_prop.push_back( &Weight );

#ifdef  __DEBUG
        Position_old.resize( nDim );
        for( unsigned int i=0 ; i< nDim ; i++ ) {
            double_prop.push_back( &( Position_old[i] ) );
        }
#endif

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

    initialize( nParticles, part.Position.size() );
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

}

// ---------------------------------------------------------------------------------------------------------------------
//Resize Particle vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::resize( unsigned int nParticles, unsigned int nDim )
{
    Position.resize( nDim );
    for( unsigned int i=0 ; i<nDim ; i++ ) {
        Position[i].resize( nParticles, 0. );
    }
#ifdef  __DEBUG
    Position_old.resize( nDim );
    for( unsigned int i=0 ; i<nDim ; i++ ) {
        Position_old[i].resize( nParticles, 0. );
    }
#endif

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

}

// ---------------------------------------------------------------------------------------------------------------------
// Remove extra capacity of Particles vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::shrink_to_fit()
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
}


// ---------------------------------------------------------------------------------------------------------------------
// Reset of Particles vectors
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
}


void Particles::cp_particle( unsigned int ipart )
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
// Copy particle iPart at the end of dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particle( unsigned int ipart, Particles &dest_parts )
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
// Insert particle iPart at dest_id in dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particle( unsigned int ipart, Particles &dest_parts, int dest_id )
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
// Insert nPart particles starting at ipart to dest_id in dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particles( unsigned int iPart, unsigned int nPart, Particles &dest_parts, int dest_id )
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
// Suppress particle iPart
// ---------------------------------------------------------------------------------------------------------------------
void Particles::erase_particle( unsigned int ipart )
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
// Suppress all particles from iPart to the end of particle array
// ---------------------------------------------------------------------------------------------------------------------
void Particles::erase_particle_trail( unsigned int ipart )
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
// Suppress npart particles from ipart
// ---------------------------------------------------------------------------------------------------------------------
void Particles::erase_particle( unsigned int ipart, unsigned int npart )
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
void Particles::swap_part( unsigned int part1, unsigned int part2 )
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


void Particles::swap_part3( unsigned int part1, unsigned int part2, unsigned int part3 )
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


void Particles::swap_part4( unsigned int part1, unsigned int part2, unsigned int part3, unsigned int part4 )
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


void Particles::swap_parts( std::vector<unsigned int> parts )
{
    // parts[0] ==> parts[1] ==> parts[2] ==> parts[parts.size()-1] ==> parts[0]

    cp_particle( parts.back() );
    translate_parts( parts );
    overwrite_part( size()-1, parts[0] );
    erase_particle( size()-1 );

}


void Particles::translate_parts( std::vector<unsigned int> parts )
{
    // parts[0] ==> parts[1] ==> parts[2] ==> parts[parts.size()-1]

    for( int icycle = parts.size()-2; icycle >=0; icycle-- ) {
        overwrite_part( parts[icycle], parts[icycle+1] );
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Move particle src_particle into dest_particle memory location, erasing dest_particle.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part( unsigned int src_particle, unsigned int dest_particle )
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
void Particles::overwrite_part( unsigned int part1, unsigned int part2, unsigned int N )
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
void Particles::overwrite_part( unsigned int part1, Particles &dest_parts, unsigned int part2 )
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
void Particles::overwrite_part( unsigned int part1, Particles &dest_parts, unsigned int part2, unsigned int N )
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
void Particles::swap_part( unsigned int part1, unsigned int part2, unsigned int N )
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
void Particles::push_to_end( unsigned int iPart )
{

}

// ---------------------------------------------------------------------------------------------------------------------
// Create a new particle at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::create_particle()
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
void Particles::create_particles( int nAdditionalParticles )
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

//MESSAGE("create2");
}

// ---------------------------------------------------------------------------------------------------------------------
//! This method compresses the particles vectors according to the provided mask
//! between istart and iend
//! The mask determines which particles to keep (>= 0) and which to delete (< 0)
//! Particle order is kept in case the vector is sorted
//! Warning: This method do not update count, first_index and last_index in Species
// ---------------------------------------------------------------------------------------------------------------------
void Particles::compressParticles( int istart, int iend, vector <int> & mask ) {

    // unsigned int ipp = istart;
    // for (unsigned int ip = (unsigned int)(istart+1) ; ip < (unsigned int) iend ; ip ++) {
    //     if (mask[ipp] < 0) {
    //         if (mask[ip] >= 0) {
    //             overwrite_part( ip, ipp);
    //             mask[ip] = -1;
    //             mask[ipp] = 1;
    //             cell_keys[ipp] = cell_keys[ip];
    //             ipp++;
    //         }
    //     } else {
    //         ipp++;
    //     }
    // }

    unsigned int idest = istart;
    unsigned int ip = istart;
    while (ip < (unsigned int) iend) {
        if (mask[idest] < 0) {
            if (mask[ip] >= 0) {
                overwrite_part( ip, idest);
                mask[ip] = -1;
                mask[idest] = 1;
                cell_keys[idest] = cell_keys[ip];
                idest++;
            } else {
                ip++;
            }
        } else {
            idest++;
            ip = idest;
        }
    }


    // At the end we resize particles
    resize(idest);
    cell_keys.resize(idest);
}

// ---------------------------------------------------------------------------------------------------------------------
// Create nParticles new particles at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
//void Particles::create_particles(int nAdditionalParticles )
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
bool Particles::is_part_in_domain( unsigned int ipart, Patch *patch )
{
    for( unsigned int i=0; i<Position.size(); i++ ) {
        if( Position[i][ipart] <  patch->getDomainLocalMin( i ) ) {
            return false;
        }
        if( Position[i][ipart] >= patch->getDomainLocalMax( i ) ) {
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
                swap_part( iPart, jPart );
                stop = false;
            }
        }
        jPart++;
    } while( !stop );

}

#ifdef __DEBUG
bool Particles::test_move( int iPartStart, int iPartEnd, Params &params )
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
