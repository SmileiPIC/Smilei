#include "Particles.h"

#include <cstring>
#include <iostream>

#include "PicParams.h"

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
	//if (nParticles > Weight.capacity()) {
	//	WARNING("You should increase c_part_max in specie namelist");
	//}

    if (Weight.size()==0) {
	float c_part_max = 1.0;
	//reserve( round( params->species_param[speciesNumber].c_part_max * nParticles ), nDim );
	reserve( round( c_part_max * nParticles ), nDim );
    }

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
// Insert particle iPart at dest_id in dest_parts
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
// Insert first nPart particles at dest_id in dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particles(int nPart, Particles &dest_parts, int dest_id )
{
    for (unsigned int i=0; i<Position.size(); i++) {
        dest_parts.Position[i].insert( dest_parts.Position[i].begin() + dest_id, Position[i].begin(), Position[i].begin()+nPart );
        dest_parts.Position_old[i].insert( dest_parts.Position_old[i].begin() + dest_id, Position_old[i].begin(), Position_old[i].begin()+nPart );
    }

    for (unsigned int i=0; i<3; i++) {
        dest_parts.Momentum[i].insert( dest_parts.Momentum[i].begin() + dest_id, Momentum[i].begin(), Momentum[i].begin()+nPart );
    }
    dest_parts.Weight.insert( dest_parts.Weight.begin() + dest_id, Weight.begin(), Weight.begin()+nPart );
    dest_parts.Charge.insert( dest_parts.Charge.begin() + dest_id, Charge.begin(), Charge.begin()+nPart );

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
// Suppress all particles from iPart to the end of particle array  
// ---------------------------------------------------------------------------------------------------------------------
void Particles::erase_particle_trail(int ipart)
{
    for (unsigned int i=0; i<Position.size(); i++) {
        Position[i].erase(Position[i].begin()+ipart,Position[i].end() );
        Position_old[i].erase(Position_old[i].begin()+ipart,Position_old[i].end() );
    }

    for (unsigned int i=0; i<3; i++) {
        Momentum[i].erase( Momentum[i].begin()+ipart,Momentum[i].end() );
    }
    Weight.erase( Weight.begin()+ipart,Weight.end() );
    Charge.erase( Charge.begin()+ipart,Charge.end() );
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
// Move particle part1 into part2 memory location, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part1D(int part1, int part2)
{
            Position[0][part2]     = Position[0][part1];
        Position_old[0][part2] = Position_old[0][part1];
            Momentum[0][part2] =     Momentum[0][part1];
            Momentum[1][part2] =     Momentum[1][part1];
            Momentum[2][part2] =     Momentum[2][part1];
              Charge[part2]      =        Charge[part1];
              Weight[part2]      =        Weight[part1];      
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1 into part2 memory location, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part2D(int part1, int part2)
{
            Position[0][part2]     = Position[0][part1];
            Position[1][part2]     = Position[1][part1];
        Position_old[0][part2] = Position_old[0][part1];
        Position_old[1][part2] = Position_old[1][part1];
            Momentum[0][part2] =     Momentum[0][part1];
            Momentum[1][part2] =     Momentum[1][part1];
            Momentum[2][part2] =     Momentum[2][part1];
              Charge[part2]      =        Charge[part1];
              Weight[part2]      =        Weight[part1];      
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part2D(int part1, int part2, int N)
{
    unsigned int sizepart,sizecharge;
    sizepart = N*sizeof(Position[0][0]);
    sizecharge = N*sizeof(Charge[0]);

    memcpy(&Position[0][part2]     ,  &Position[0][part1]     , sizepart)    ;
    memcpy(&Position[1][part2]     ,  &Position[1][part1]     , sizepart)    ;
    memcpy(&Position_old[0][part2] ,  &Position_old[0][part1] , sizepart)    ;
    memcpy(&Position_old[1][part2] ,  &Position_old[1][part1] , sizepart)    ;
    memcpy(&Momentum[0][part2]     ,  &Momentum[0][part1]     , sizepart)    ;
    memcpy(&Momentum[1][part2]     ,  &Momentum[1][part1]     , sizepart)    ;
    memcpy(&Momentum[2][part2]     ,  &Momentum[2][part1]     , sizepart)    ;
    memcpy(&Charge[part2]          ,  &Charge[part1]          , sizecharge)    ;
    memcpy(&Weight[part2]          ,  &Weight[part1]          , sizepart)    ;      
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1 into part2 memory location of dest vector, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part2D(int part1, Particles &dest_parts, int part2)
{

            dest_parts.Position[0][part2]     = Position[0][part1];
            dest_parts.Position[1][part2]     = Position[1][part1];
            dest_parts.Position_old[0][part2] = Position_old[0][part1];
            dest_parts.Position_old[1][part2] = Position_old[1][part1];
            dest_parts.Momentum[0][part2] =     Momentum[0][part1];
            dest_parts.Momentum[1][part2] =     Momentum[1][part1];
            dest_parts.Momentum[2][part2] =     Momentum[2][part1];
            dest_parts.Charge[part2]      =     Charge[part1];
            dest_parts.Weight[part2]      =     Weight[part1];      
    }

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location of dest vector, erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part2D(int part1, Particles &dest_parts, int part2, int N)
{
    unsigned int sizepart,sizecharge;
    sizepart = N*sizeof(Position[0][0]);
    sizecharge = N*sizeof(Charge[0]);

    memcpy(&dest_parts.Position[0][part2]     ,  &Position[0][part1]     , sizepart)    ;
    memcpy(&dest_parts.Position[1][part2]     ,  &Position[1][part1]     , sizepart)    ;
    memcpy(&dest_parts.Position_old[0][part2] ,  &Position_old[0][part1] , sizepart)    ;
    memcpy(&dest_parts.Position_old[1][part2] ,  &Position_old[1][part1] , sizepart)    ;
    memcpy(&dest_parts.Momentum[0][part2]     ,  &Momentum[0][part1]     , sizepart)    ;
    memcpy(&dest_parts.Momentum[1][part2]     ,  &Momentum[1][part1]     , sizepart)    ;
    memcpy(&dest_parts.Momentum[2][part2]     ,  &Momentum[2][part1]     , sizepart)    ;
    memcpy(&dest_parts.Charge[part2]          ,  &Charge[part1]         , sizecharge)    ;
    memcpy(&dest_parts.Weight[part2]          ,  &Weight[part1]         , sizepart)    ;      
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part1D(int part1, int part2, int N)
{
        for (unsigned int j=0; j< N; j++) {
                Position[0][part2+j]   =    Position[0][part1+j];
            Position_old[0][part2+j] =  Position_old[0][part1+j];
                Momentum[0][part2+j]     =  Momentum[0][part1+j];
                Momentum[1][part2+j]     =  Momentum[1][part1+j];
                Momentum[2][part2+j]     =  Momentum[2][part1+j];
                     Charge[part2+j]          =  Charge[part1+j];
                     Weight[part2+j]          =  Weight[part1+j];      
        }
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location of dest vector, erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part1D(int part1, Particles &dest_parts, int part2, int N)
{
        for (unsigned int j=0; j< N; j++) {
                dest_parts.Position[0][part2+j]     = Position[0][part1+j];
            dest_parts.Position_old[0][part2+j] = Position_old[0][part1+j];
                dest_parts.Momentum[0][part2+j] =     Momentum[0][part1+j];
                dest_parts.Momentum[1][part2+j] =     Momentum[1][part1+j];
                dest_parts.Momentum[2][part2+j] =     Momentum[2][part1+j];
                     dest_parts.Charge[part2+j]      =     Charge[part1+j];
                     dest_parts.Weight[part2+j]      =     Weight[part1+j];      
        }
}

// ---------------------------------------------------------------------------------------------------------------------
// Exchange N particles part1->part1+N & part2->part2+N memory location
// ---------------------------------------------------------------------------------------------------------------------
void Particles::swap_part(int part1, int part2, int N)
{
    double* buffer[N];
    unsigned int sizepart,sizecharge;

    sizepart = N*sizeof(Position[0][0]);
    sizecharge = N*sizeof(Charge[0]);

    for (unsigned int i=0; i<Position.size(); i++) {
        memcpy(buffer,&Position[i][part1], sizepart);
        memcpy(&Position[i][part1],&Position[i][part2], sizepart);
        memcpy(&Position[i][part2],buffer, sizepart);

        memcpy(buffer,&Position_old[i][part1], sizepart);
        memcpy(&Position_old[i][part1],&Position_old[i][part2], sizepart);
        memcpy(&Position_old[i][part2],buffer, sizepart);
    }
    for (unsigned int i=0; i<3; i++){
        memcpy(buffer,&Momentum[i][part1], sizepart);
        memcpy(&Momentum[i][part1],&Momentum[i][part2], sizepart);
        memcpy(&Momentum[i][part2],buffer, sizepart);
    }
        memcpy(buffer,&Charge[part1], sizecharge);
        memcpy(&Charge[part1],&Charge[part2], sizecharge);
        memcpy(&Charge[part2],buffer, sizecharge);

        memcpy(buffer,&Weight[part1], sizepart);
        memcpy(&Weight[part1],&Weight[part2], sizepart);
        memcpy(&Weight[part2],buffer, sizepart);
    
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

