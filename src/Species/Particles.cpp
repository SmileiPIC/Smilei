#include "Particles.h"

#include <cstring>
#include <iostream>

#include "PicParams.h"
#include "SmileiMPI.h"

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Particle
// ---------------------------------------------------------------------------------------------------------------------
Particles::Particles()
{
    Position.resize(0);
    Position_old.resize(0);
    Momentum.resize(0);
    isTestParticles = false;
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
void Particles::initialize( int nParticles, PicParams &params)
{
    //if (nParticles > Weight.capacity()) {
    //    WARNING("You should increase c_part_max in specie namelist");
    //}
    
    if (Weight.size()==0) {
        float c_part_max = 1.0;
        //reserve( round( params->species_param[speciesNumber].c_part_max * nParticles ), params.nDim );
        reserve( round( c_part_max * nParticles ), params.nDim_particle );
    }
    
    Position.resize(params.nDim_particle);
    Position_old.resize(params.nDim_particle);
    for (unsigned int i=0 ; i< params.nDim_particle ; i++) {
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
// Create nParticles null particles of nDim size
// ---------------------------------------------------------------------------------------------------------------------
void Particles::initialize( int nParticles, PicParams &params, int speciesNumber)
{
    initialize( nParticles, params );
    
    if (params.species_param[speciesNumber].isTest) {
        isTestParticles = true;
        Id.resize(nParticles, 0);
    }

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
    
    if (isTestParticles)
        Id.reserve(n_part_max);

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
    
    if (isTestParticles)
        Id.clear();
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
    
    if (isTestParticles)
        dest_parts.Id.push_back( Id[ipart] );
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
    
    if (isTestParticles)
        dest_parts.Id.insert( dest_parts.Id.begin() + dest_id, Id[ipart] );

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
    
    if (isTestParticles)
        dest_parts.Id.insert( dest_parts.Id.begin() + dest_id, Id.begin(), Id.begin()+nPart );

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
    
    if (isTestParticles)
    Id.erase( Id.begin()+ipart );

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
    
    if (isTestParticles)
        Id.erase( Id.begin()+ipart,Id.end() );

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
    
    if (isTestParticles)
        cout << Id[iPart] << endl;
}


// ---------------------------------------------------------------------------------------------------------------------
// Print parameters of particle iPart
// ---------------------------------------------------------------------------------------------------------------------
ostream& operator << (ostream& out, const Particles& particles) {
    for (unsigned int iPart=0;iPart<particles.Weight.size();iPart++) {
        
        for (unsigned int i=0; i<particles.Position.size(); i++) {
            out << particles.Position[i][iPart] << " ";
            out << particles.Position_old[i][iPart] << " ";
        }
        for (unsigned int i=0; i<3; i++)
            out << particles.Momentum[i][iPart] << " ";
        out << particles.Weight[iPart] << " ";
        out << particles.Charge[iPart] << endl;;
        
        if (particles.isTestParticles)
            out << particles.Id[iPart] << endl;
    }
    
    return (out);
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
    
    if (isTestParticles)
        std::swap( Id[part1], Id[part2] );
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1 into part2 memory location, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part1D(int part1, int part2)
{
    Position    [0][part2] = Position    [0][part1];
    Position_old[0][part2] = Position_old[0][part1];
    Momentum    [0][part2] = Momentum    [0][part1];
    Momentum    [1][part2] = Momentum    [1][part1];
    Momentum    [2][part2] = Momentum    [2][part1];
    Charge         [part2] = Charge         [part1];
    Weight         [part2] = Weight         [part1];
    
    if (isTestParticles)
        Id[part2] = Id[part1];
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1 into part2 memory location, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part2D(int part1, int part2)
{
    Position    [0][part2] = Position    [0][part1];
    Position    [1][part2] = Position    [1][part1];
    Position_old[0][part2] = Position_old[0][part1];
    Position_old[1][part2] = Position_old[1][part1];
    Momentum    [0][part2] = Momentum    [0][part1];
    Momentum    [1][part2] = Momentum    [1][part1];
    Momentum    [2][part2] = Momentum    [2][part1];
    Charge         [part2] = Charge         [part1];
    Weight         [part2] = Weight         [part1];
    
    if (isTestParticles)
        Id[part2] = Id[part1];

}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part2D(int part1, int part2, int N)
{
    unsigned int sizepart,sizecharge;
    sizepart = N*sizeof(Position[0][0]);
    sizecharge = N*sizeof(Charge[0]);
    
    memcpy( &Position    [0][part2] , &Position    [0][part1] , sizepart  );
    memcpy( &Position    [1][part2] , &Position    [1][part1] , sizepart  );
    memcpy( &Position_old[0][part2] , &Position_old[0][part1] , sizepart  );
    memcpy( &Position_old[1][part2] , &Position_old[1][part1] , sizepart  );
    memcpy( &Momentum    [0][part2] , &Momentum    [0][part1] , sizepart  );
    memcpy( &Momentum    [1][part2] , &Momentum    [1][part1] , sizepart  );
    memcpy( &Momentum    [2][part2] , &Momentum    [2][part1] , sizepart  );
    memcpy( &Charge         [part2] , &Charge         [part1] , sizecharge);
    memcpy( &Weight         [part2] , &Weight         [part1] , sizepart  );
      
    if (isTestParticles)
        memcpy(&Id[part2]          ,  &Id[part1]              , sizecharge)    ;
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1 into part2 memory location of dest vector, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part2D(int part1, Particles &dest_parts, int part2)
{
    dest_parts.Position    [0][part2] = Position    [0][part1];
    dest_parts.Position    [1][part2] = Position    [1][part1];
    dest_parts.Position_old[0][part2] = Position_old[0][part1];
    dest_parts.Position_old[1][part2] = Position_old[1][part1];
    dest_parts.Momentum    [0][part2] = Momentum    [0][part1];
    dest_parts.Momentum    [1][part2] = Momentum    [1][part1];
    dest_parts.Momentum    [2][part2] = Momentum    [2][part1];
    dest_parts.Charge         [part2] = Charge         [part1];
    dest_parts.Weight         [part2] = Weight         [part1];
    
    if (isTestParticles)
        dest_parts.Id[part2] = Id[part1];
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location of dest vector, erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part2D(int part1, Particles &dest_parts, int part2, int N)
{
    unsigned int sizepart,sizecharge;
    sizepart = N*sizeof(Position[0][0]);
    sizecharge = N*sizeof(Charge[0]);
    
    memcpy( &dest_parts.Position    [0][part2] , &Position    [0][part1] , sizepart  );
    memcpy( &dest_parts.Position    [1][part2] , &Position    [1][part1] , sizepart  );
    memcpy( &dest_parts.Position_old[0][part2] , &Position_old[0][part1] , sizepart  );
    memcpy( &dest_parts.Position_old[1][part2] , &Position_old[1][part1] , sizepart  );
    memcpy( &dest_parts.Momentum    [0][part2] , &Momentum    [0][part1] , sizepart  );
    memcpy( &dest_parts.Momentum    [1][part2] , &Momentum    [1][part1] , sizepart  );
    memcpy( &dest_parts.Momentum    [2][part2] , &Momentum    [2][part1] , sizepart  );
    memcpy( &dest_parts.Charge         [part2] , &Charge         [part1] , sizecharge);
    memcpy( &dest_parts.Weight         [part2] , &Weight         [part1] , sizepart  );
    
    if (isTestParticles)
        memcpy(&dest_parts.Id[part2],  &Id[part1], sizecharge);
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part1D(int part1, int part2, int N)
{
    for (unsigned int j=0; j< (unsigned int) N; j++) {
        Position    [0][part2+j] =  Position    [0][part1+j];
        Position_old[0][part2+j] =  Position_old[0][part1+j];
        Momentum    [0][part2+j] =  Momentum    [0][part1+j];
        Momentum    [1][part2+j] =  Momentum    [1][part1+j];
        Momentum    [2][part2+j] =  Momentum    [2][part1+j];
        Charge         [part2+j] =  Charge         [part1+j];
        Weight         [part2+j] =  Weight         [part1+j];
    }
    
    if (isTestParticles)
        for (int j=0; j< N; j++)
            Id[part2+j] = Id[part1+j];
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location of dest vector, erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part1D(int part1, Particles &dest_parts, int part2, int N)
{
    for (unsigned int j=0; j< (unsigned int) N; j++) {
        dest_parts.Position    [0][part2+j] = Position    [0][part1+j];
        dest_parts.Position_old[0][part2+j] = Position_old[0][part1+j];
        dest_parts.Momentum    [0][part2+j] = Momentum    [0][part1+j];
        dest_parts.Momentum    [1][part2+j] = Momentum    [1][part1+j];
        dest_parts.Momentum    [2][part2+j] = Momentum    [2][part1+j];
        dest_parts.Charge         [part2+j] = Charge         [part1+j];
        dest_parts.Weight         [part2+j] = Weight         [part1+j];
    }
    
    if (isTestParticles)
        for (int j=0; j< N; j++)
            dest_parts.Id[part2+j] = Id[part1+j];
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
    
    if (isTestParticles) {
        memcpy(buffer,&Id[part1], sizecharge);
        memcpy(&Id[part1],&Id[part2], sizecharge);
        memcpy(&Id[part2],buffer, sizecharge);
    }
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
    
    if (isTestParticles)
    Id.push_back(0);
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
    
    if (isTestParticles)
    Id.resize(nParticles+nAdditionalParticles,0);

}

// ---------------------------------------------------------------------------------------------------------------------
// Test if ipart is in the local MPI subdomain
//---------------------------------------------------------------------------------------------------------------------
bool Particles::is_part_in_domain(int ipart, SmileiMPI* smpi)
{
    for (unsigned int i=0; i<Position.size(); i++) {
        if (Position[i][ipart] < smpi->getDomainLocalMin(i) ) return false;
        if (Position[i][ipart] >= smpi->getDomainLocalMax(i) ) return false;
    }
    return true;
}


void Particles::sortById() {
    if (!isTestParticles) {
        ERROR("Impossible");
        return;
    }
    int nParticles(Weight.size());
    
    bool stop;
    int jPart(0);
    do {
        stop = true;
        for ( int iPart = nParticles-1 ; iPart > jPart ; --iPart ) {
            if ( Id[iPart] < Id[iPart-1] ) {
                swap_part(iPart,jPart);
                stop = false;
            }
        }
        jPart++;
    } while(!stop);
    
}
