#include "Particles.h"
#include "PicParams.h"

#include <iostream>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Particle
// ---------------------------------------------------------------------------------------------------------------------
Particles::Particles()
{
    Position.resize(0);
    Position_old.resize(0);
    Momentum.resize(0);
}

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Particle
// ---------------------------------------------------------------------------------------------------------------------
Particles::~Particles()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// Create nParticles null particles of nDim size
// ---------------------------------------------------------------------------------------------------------------------
void Particles::initialize( int nParticles, int nDim )
{
	Position.resize(nDim);
    Position_old.resize(nDim);
    for (int i=0 ; i< nDim ; i++) {
        Position[i].resize(nParticles, 0.);
        Position_old[i].resize(nParticles, 0.);
    }
    Momentum.resize(3);
    for (int i=0 ; i< 3 ; i++) {
        Momentum[i].resize(nParticles, 0.);
    }
    Weight.resize(nParticles, 0.);
    Charge.resize(nParticles, 0);

}


// ---------------------------------------------------------------------------------------------------------------------
// Set capacity of Particles vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::reserve( unsigned int n_part_max, int nDim )
{
    Position.resize(nDim);
    Position_old.resize(nDim);
    for (int i=0 ; i< nDim ; i++) {
        Position[i].reserve(n_part_max);
        Position_old[i].reserve(n_part_max);
    }
    Momentum.resize(3);
    for (int i=0 ; i< 3 ; i++) {
        Momentum[i].reserve(n_part_max);
    }
    Weight.reserve(n_part_max);
    Charge.reserve(n_part_max);

}

// ---------------------------------------------------------------------------------------------------------------------
// Reset of Particles vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::clear()
{
    for (unsigned int i=0 ; i< Position.size() ; i++) {
        Position[i].clear();
        Position_old[i].clear();
    }
    for (int i=0 ; i< 3 ; i++) {
        Momentum[i].clear();
    }
    Weight.clear();
    Charge.clear();

}

// ---------------------------------------------------------------------------------------------------------------------
// Copy particle iPart at the end of dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particle(int ipart, Particles &dest_parts )
{
    for (unsigned int i=0; i<Position.size(); i++) {
        dest_parts.Position[i].push_back(Position[i][ipart]);
        dest_parts.Position_old[i].push_back(Position_old[i][ipart]);
    }

    for (unsigned int i=0; i<3; i++) {
        dest_parts.Momentum[i].push_back( Momentum[i][ipart] );
    }
    dest_parts.Weight.push_back( Weight[ipart] );
    dest_parts.Charge.push_back( Charge[ipart] );

}

// ---------------------------------------------------------------------------------------------------------------------
// Copy particle iPart at dest_id in dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particle(int ipart, Particles &dest_parts, int dest_id )
{
    for (unsigned int i=0; i<Position.size(); i++) {
        dest_parts.Position[i].insert( dest_parts.Position[i].begin() + dest_id, Position[i][ipart] );
        dest_parts.Position_old[i].insert( dest_parts.Position_old[i].begin() + dest_id, Position_old[i][ipart] );
    }

    for (unsigned int i=0; i<3; i++) {
        dest_parts.Momentum[i].insert( dest_parts.Momentum[i].begin() + dest_id, Momentum[i][ipart] );
    }
    dest_parts.Weight.insert( dest_parts.Weight.begin() + dest_id, Weight[ipart] );
    dest_parts.Charge.insert( dest_parts.Charge.begin() + dest_id, Charge[ipart] );

}

// ---------------------------------------------------------------------------------------------------------------------
// Suppress particle iPart
// ---------------------------------------------------------------------------------------------------------------------
void Particles::erase_particle(int ipart )
{
    for (unsigned int i=0; i<Position.size(); i++) {

        Position[i].erase(Position[i].begin()+ipart);
        Position_old[i].erase(Position_old[i].begin()+ipart);
    }

    for (unsigned int i=0; i<3; i++) {
        Momentum[i].erase( Momentum[i].begin()+ipart );
    }
    Weight.erase( Weight.begin()+ipart );
    Charge.erase( Charge.begin()+ipart );
}

// ---------------------------------------------------------------------------------------------------------------------
// Print parameters of particle iPart
// ---------------------------------------------------------------------------------------------------------------------
void Particles::print(int iPart) {
    for (unsigned int i=0; i<Position.size(); i++) {
        cout << Position[i][iPart] << " ";
        cout << Position_old[i][iPart] << " ";
    }
    for (unsigned int i=0; i<3; i++)
        cout << Momentum[i][iPart] << " ";
    cout << Weight[iPart] << " ";
    cout << Charge[iPart] << endl;;
}


// ---------------------------------------------------------------------------------------------------------------------
// Exchange particles part1 & part2 memory location
// ---------------------------------------------------------------------------------------------------------------------
void Particles::swap_part(int part1, int part2)
{
    for (unsigned int i=0; i<Position.size(); i++) {
        std::swap( Position[i][part1], Position[i][part2] );
        std::swap( Position_old[i][part1], Position_old[i][part2] );
    }
    for (unsigned int i=0; i<3; i++)
        std::swap( Momentum[i][part1], Momentum[i][part2] );
    std::swap( Charge[part1], Charge[part2] );
    std::swap( Weight[part1], Weight[part2] );

}

// ---------------------------------------------------------------------------------------------------------------------
// Move iPart at the end of vectors (to do for MPI)
// ---------------------------------------------------------------------------------------------------------------------
void Particles::push_to_end(int iPart )
{

}

// ---------------------------------------------------------------------------------------------------------------------
// Create a new particle at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::create_particle()
{
    for (unsigned int i=0; i<Position.size(); i++) {
        Position[i].push_back(0.);
        Position_old[i].push_back(0.);
    }

    for (unsigned int i=0; i<3; i++) {
        Momentum[i].push_back(0.);
    }
    Weight.push_back(0.);
    Charge.push_back(0);

}

// ---------------------------------------------------------------------------------------------------------------------
// Create nParticles new particles at the end of vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::create_particles(int nAdditionalParticles )
{
    int nParticles = size();
    for (unsigned int i=0; i<Position.size(); i++) {
        Position[i].resize(nParticles+nAdditionalParticles,0.);
        Position_old[i].resize(nParticles+nAdditionalParticles,0.);
    }

    for (unsigned int i=0; i<3; i++) {
        Momentum[i].resize(nParticles+nAdditionalParticles,0.);
    }
    Weight.resize(nParticles+nAdditionalParticles,0.);
    Charge.resize(nParticles+nAdditionalParticles,0);

}

