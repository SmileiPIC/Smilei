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
tracked(false)
{
    Position.resize(0);
    Position_old.resize(0);
    Momentum.resize(0);
    isTest = false;
    isRadReaction = false;

    double_prop.resize(0);
    short_prop.resize(0);
    uint_prop.resize(0);
}

// ---------------------------------------------------------------------------------------------------------------------
// Create nParticles null particles of nDim size
// ---------------------------------------------------------------------------------------------------------------------
void Particles::initialize(unsigned int nParticles, unsigned int nDim)
{
    //if (nParticles > Weight.capacity()) {
    //    WARNING("You should increase c_part_max in specie namelist");
    //}
    if (Weight.size()==0) {
        float c_part_max =1.2;
        //float c_part_max = part.c_part_max;
        //float c_part_max = params.species_param[0].c_part_max;
        reserve( round( c_part_max * nParticles ), nDim );
    }
    
    resize(nParticles, nDim);
    
    if ( double_prop.empty() ) { // do this just once 
        for (unsigned int i=0 ; i< Position.size() ; i++)
            double_prop.push_back( &(Position[i]) );
#ifdef  __DEBUG
        for (unsigned int i=0 ; i< Position_old.size() ; i++)
            double_prop.push_back( &(Position_old[i]) );
#endif
        for (unsigned int i=0 ; i< 3 ; i++)
            double_prop.push_back( &(Momentum[i]) );
        double_prop.push_back( &Weight );
        short_prop.push_back( &Charge );
        if (tracked) {
            uint_prop.push_back( &Id );
        }
        
        if (isRadReaction) {
            double_prop.push_back( &Chi );
        }
        
    }
    
}

// copy properties from another Particles
void Particles::initialize(unsigned int nParticles, Particles &part)
{
    isTest=part.isTest;
    
    tracked=part.tracked;
    
    isRadReaction=part.isRadReaction;
    
    initialize(nParticles, part.Position.size());
}



// ---------------------------------------------------------------------------------------------------------------------
// Set capacity of Particles vectors
// ---------------------------------------------------------------------------------------------------------------------
void Particles::reserve( unsigned int n_part_max, unsigned int nDim )
{
    Position.resize(nDim);
    Position_old.resize(nDim);
    for (unsigned int i=0 ; i< nDim ; i++) {
        Position[i].reserve(n_part_max);
        Position_old[i].reserve(n_part_max);
    }
    Momentum.resize(3);
    for (unsigned int i=0 ; i< 3 ; i++) {
        Momentum[i].reserve(n_part_max);
    }
    Weight.reserve(n_part_max);
    Charge.reserve(n_part_max);
    
    if (tracked)
        Id.reserve(n_part_max);
    
    if (isRadReaction)
        Chi.reserve(n_part_max);

}

void Particles::resize( unsigned int nParticles, unsigned int nDim )
{
    Position.resize(nDim);
    Position_old.resize(nDim);
    for (unsigned int i=0 ; i<nDim ; i++) {
        Position[i].resize(nParticles, 0.);
        Position_old[i].resize(nParticles, 0.);
    }
    
    Momentum.resize(3);
    for (unsigned int i=0 ; i< 3 ; i++) {
        Momentum[i].resize(nParticles, 0.);
    }
    
    Weight.resize(nParticles, 0.);
    Charge.resize(nParticles, 0);
    
    if (tracked) {
        Id.resize(nParticles, 0);
    }
    
    if (isRadReaction) {
        Chi.resize(nParticles, 0.);
    }
}

void Particles::shrink_to_fit( unsigned int nDim )
{
    Position.resize(nDim);
    Position_old.resize(nDim);
    for (unsigned int i=0 ; i<nDim; i++) {
        std::vector<double>(Position[i]).swap(Position[i]);
        std::vector<double>(Position_old[i]).swap(Position_old[i]);
    }
    
    Momentum.resize(3);
    for (unsigned int i=0 ; i< 3 ; i++) {
        std::vector<double>(Momentum[i]).swap(Momentum[i]);
    }
    
    std::vector<double>(Weight).swap(Weight);
    std::vector<short>(Charge).swap(Charge);
    
    if (tracked) {
        std::vector<unsigned int>(Id).swap(Id);
    }
    
    if (isRadReaction) {
        std::vector<double>(Chi).swap(Chi);
    }
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
    for (unsigned int i=0 ; i< 3 ; i++) {
        Momentum[i].clear();
    }
    Weight.clear();
    Charge.clear();
    
    if (tracked)
        Id.clear();
    
    if (isRadReaction)
        Chi.clear();
}

// ---------------------------------------------------------------------------------------------------------------------
// Copy particle iPart at the end of dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particle(unsigned int ipart, Particles &dest_parts )
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
    
    if (tracked)
        dest_parts.Id.push_back( Id[ipart] );
    
    if (isRadReaction)
        dest_parts.Chi.push_back( Chi[ipart] );
}

// ---------------------------------------------------------------------------------------------------------------------
// Insert particle iPart at dest_id in dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particle(unsigned int ipart, Particles &dest_parts, int dest_id )
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
    
    if (tracked)
        dest_parts.Id.insert( dest_parts.Id.begin() + dest_id, Id[ipart] );
    
    if (isRadReaction)
        dest_parts.Chi.insert( dest_parts.Chi.begin() + dest_id, Chi[ipart] );

}

// ---------------------------------------------------------------------------------------------------------------------
// Insert nPart particles starting at ipart to dest_id in dest_parts
// ---------------------------------------------------------------------------------------------------------------------
void Particles::cp_particles(unsigned int iPart, unsigned int nPart, Particles &dest_parts, int dest_id )
{
    for (unsigned int i=0; i<Position.size(); i++) {
        dest_parts.Position[i].insert( dest_parts.Position[i].begin() + dest_id, Position[i].begin()+iPart, Position[i].begin()+iPart+nPart );
        dest_parts.Position_old[i].insert( dest_parts.Position_old[i].begin() + dest_id, Position_old[i].begin()+iPart, Position_old[i].begin()+iPart+nPart );
    }
    
    for (unsigned int i=0; i<3; i++) {
        dest_parts.Momentum[i].insert( dest_parts.Momentum[i].begin() + dest_id, Momentum[i].begin()+iPart, Momentum[i].begin()+iPart+nPart );
    }
    dest_parts.Weight.insert( dest_parts.Weight.begin() + dest_id, Weight.begin()+iPart, Weight.begin()+iPart+nPart );
    dest_parts.Charge.insert( dest_parts.Charge.begin() + dest_id, Charge.begin()+iPart, Charge.begin()+iPart+nPart );
    
    if (tracked)
        dest_parts.Id.insert( dest_parts.Id.begin() + dest_id, Id.begin()+iPart, Id.begin()+iPart+nPart );
    
    if (isRadReaction)
        dest_parts.Chi.insert( dest_parts.Chi.begin() + dest_id, Chi.begin()+iPart, Chi.begin()+iPart+nPart );

}

// ---------------------------------------------------------------------------------------------------------------------
// Suppress particle iPart
// ---------------------------------------------------------------------------------------------------------------------
void Particles::erase_particle(unsigned int ipart )
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
    
    if (tracked)
        Id.erase( Id.begin()+ipart );
    
    if (isRadReaction)
        Chi.erase( Chi.begin()+ipart );

}

// ---------------------------------------------------------------------------------------------------------------------
// Suppress all particles from iPart to the end of particle array  
// ---------------------------------------------------------------------------------------------------------------------
void Particles::erase_particle_trail(unsigned int ipart)
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
    
    if (tracked)
        Id.erase( Id.begin()+ipart,Id.end() );
    
    if (isRadReaction)
        Chi.erase( Chi.begin()+ipart,Chi.end() );

}
// ---------------------------------------------------------------------------------------------------------------------
// Suppress npart particles from ipart  
// ---------------------------------------------------------------------------------------------------------------------
void Particles::erase_particle(unsigned int ipart, unsigned int npart)
{
    for (unsigned int i=0; i<Position.size(); i++) {
        Position[i].erase(Position[i].begin()+ipart,Position[i].begin()+ipart+npart );
        Position_old[i].erase(Position_old[i].begin()+ipart,Position_old[i].begin()+ipart+npart );
    }

    for (unsigned int i=0; i<3; i++) {
        Momentum[i].erase( Momentum[i].begin()+ipart,Momentum[i].begin()+ipart+npart );
    }
    Weight.erase( Weight.begin()+ipart,Weight.begin()+ipart+npart );
    Charge.erase( Charge.begin()+ipart,Charge.begin()+ipart+npart );
    
    if (tracked)
        Id.erase( Id.begin()+ipart,Id.begin()+ipart+npart );
}

// ---------------------------------------------------------------------------------------------------------------------
// Print parameters of particle iPart
// ---------------------------------------------------------------------------------------------------------------------
void Particles::print(unsigned int iPart) {
    for (unsigned int i=0; i<Position.size(); i++) {
        cout << Position[i][iPart] << " ";
        cout << Position_old[i][iPart] << " ";
    }
    for (unsigned int i=0; i<3; i++)
        cout << Momentum[i][iPart] << " ";
    cout << Weight[iPart] << " ";
    cout << Charge[iPart] << endl;;
    
    if (tracked)
        cout << Id[iPart] << endl;
    
    if (isRadReaction)
        cout << Chi[iPart] << endl;
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
        
        if (particles.tracked)
            out << particles.Id[iPart] << endl;
        
        if (particles.isRadReaction)
            out << particles.Chi[iPart] << endl;
    }
    
    return (out);
}


// ---------------------------------------------------------------------------------------------------------------------
// Exchange particles part1 & part2 memory location
// ---------------------------------------------------------------------------------------------------------------------
void Particles::swap_part(unsigned int part1, unsigned int part2)
{
    for (unsigned int i=0; i<Position.size(); i++) {
        std::swap( Position[i][part1], Position[i][part2] );
        std::swap( Position_old[i][part1], Position_old[i][part2] );
    }
    for (unsigned int i=0; i<3; i++)
        std::swap( Momentum[i][part1], Momentum[i][part2] );
    std::swap( Charge[part1], Charge[part2] );
    std::swap( Weight[part1], Weight[part2] );
    
    if (tracked)
        std::swap( Id[part1], Id[part2] );
    
    if (isRadReaction)
        std::swap( Chi[part1], Chi[part2] );
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1 into part2 memory location, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part(unsigned int part1, unsigned int part2)
{
    for (unsigned int i=0; i<Position.size(); i++) {
        Position[i][part2]     = Position[i][part1];
        Position_old[i][part2] = Position_old[i][part1];
    }
    Momentum[0][part2] = Momentum[0][part1];
    Momentum[1][part2] = Momentum[1][part1];
    Momentum[2][part2] = Momentum[2][part1];
    Charge[part2]      = Charge[part1];
    Weight[part2]      = Weight[part1];      
    
    if (tracked)
        Id[part2] = Id[part1];
    
    if (isRadReaction)
        Chi[part2] = Chi[part1];
}


// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part(unsigned int part1, unsigned int part2, unsigned int N)
{
    unsigned int sizepart = N*sizeof(Position[0][0]);
    unsigned int sizecharge = N*sizeof(Charge[0]);

    for (unsigned int i=0; i<Position.size(); i++) {
        memcpy(&Position[i][part2]     ,  &Position[i][part1]     , sizepart)    ;
        memcpy(&Position_old[i][part2] ,  &Position_old[i][part1] , sizepart)    ;
    }
    memcpy(&Momentum[0][part2]     ,  &Momentum[0][part1]     , sizepart)    ;
    memcpy(&Momentum[1][part2]     ,  &Momentum[1][part1]     , sizepart)    ;
    memcpy(&Momentum[2][part2]     ,  &Momentum[2][part1]     , sizepart)    ;
    memcpy(&Charge[part2]          ,  &Charge[part1]          , sizecharge)    ;
    memcpy(&Weight[part2]          ,  &Weight[part1]          , sizepart)    ;      

    if (tracked) {
        unsigned int sizeid = N*sizeof(Id[0]);
        memcpy(&Id[part2]          ,  &Id[part1]              , sizeid);
    }
    
    if (isRadReaction)
        memcpy(&Chi[part2]          ,  &Chi[part1]              , sizepart);
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1 into part2 memory location of dest vector, erasing part2.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part(unsigned int part1, Particles &dest_parts, unsigned int part2)
{
    for (unsigned int i=0; i<Position.size(); i++) {
        dest_parts.Position[i][part2]     = Position[i][part1];
        dest_parts.Position_old[i][part2] = Position_old[i][part1];
    }
    dest_parts.Momentum[0][part2] = Momentum[0][part1];
    dest_parts.Momentum[1][part2] = Momentum[1][part1];
    dest_parts.Momentum[2][part2] = Momentum[2][part1];
    dest_parts.Charge[part2]      = Charge[part1];
    dest_parts.Weight[part2]      = Weight[part1];      

    if (tracked)
        dest_parts.Id[part2] = Id[part1];
    
    if (isRadReaction)
        dest_parts.Chi[part2] = Chi[part1];
}

// ---------------------------------------------------------------------------------------------------------------------
// Move particle part1->part1+N into part2->part2+N memory location of dest vector, erasing part2->part2+N.
// ---------------------------------------------------------------------------------------------------------------------
void Particles::overwrite_part(unsigned int part1, Particles &dest_parts, unsigned int part2, unsigned int N)
{
    unsigned int sizepart = N*sizeof(Position[0][0]);
    unsigned int sizecharge = N*sizeof(Charge[0]);
    
    for (unsigned int i=0; i<Position.size(); i++) {
        memcpy(&dest_parts.Position[i][part2]     ,  &Position[i][part1]     , sizepart)    ;
        memcpy(&dest_parts.Position_old[i][part2] ,  &Position_old[i][part1] , sizepart)    ;
    }

    memcpy(&dest_parts.Momentum[0][part2]     ,  &Momentum[0][part1]     , sizepart)    ;
    memcpy(&dest_parts.Momentum[1][part2]     ,  &Momentum[1][part1]     , sizepart)    ;
    memcpy(&dest_parts.Momentum[2][part2]     ,  &Momentum[2][part1]     , sizepart)    ;
    memcpy(&dest_parts.Charge[part2]          ,  &Charge[part1]         , sizecharge)    ;
    memcpy(&dest_parts.Weight[part2]          ,  &Weight[part1]         , sizepart)    ;      

   if (tracked) {
        unsigned int sizeid = N*sizeof(Id[0]);
        memcpy(&dest_parts.Id[part2],  &Id[part1], sizeid);
    }
    
    if (isRadReaction)
        memcpy(&dest_parts.Chi[part2],  &Chi[part1], sizepart);

}


// ---------------------------------------------------------------------------------------------------------------------
// Exchange N particles part1->part1+N & part2->part2+N memory location
// ---------------------------------------------------------------------------------------------------------------------
void Particles::swap_part(unsigned int part1, unsigned int part2, unsigned int N)
{
    double* buffer[N];
    
    unsigned int sizepart = N*sizeof(Position[0][0]);
    unsigned int sizecharge = N*sizeof(Charge[0]);
    
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
    
    if (tracked) {
        unsigned int sizeid = N*sizeof(Id[0]);
        memcpy(buffer,&Id[part1], sizeid);
        memcpy(&Id[part1],&Id[part2], sizeid);
        memcpy(&Id[part2],buffer, sizeid);
    }
    
    if (isRadReaction) {
        memcpy(buffer,&Chi[part1], sizepart);
        memcpy(&Chi[part1],&Chi[part2], sizepart);
        memcpy(&Chi[part2],buffer, sizepart);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Move iPart at the end of vectors (to do for MPI)
// ---------------------------------------------------------------------------------------------------------------------
void Particles::push_to_end(unsigned int iPart )
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
    
    if (tracked) {
        Id.push_back(0);
    }
    if (isRadReaction) {
        Chi.push_back(0.);
    }
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
//    if (isRadReaction)
//        Chi.resize(nParticles+nAdditionalParticles,0.);
//
//}

// ---------------------------------------------------------------------------------------------------------------------
// Test if ipart is in the local patch
//---------------------------------------------------------------------------------------------------------------------
bool Particles::is_part_in_domain(unsigned int ipart, Patch* patch)
{
    for (unsigned int i=0; i<Position.size(); i++) {
        if (Position[i][ipart] <  patch->getDomainLocalMin(i) ) return false;
        if (Position[i][ipart] >= patch->getDomainLocalMax(i) ) return false;
    }
    return true;
}


void Particles::sortById() {
    if (!tracked) {
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

//bool Particles::test_move( int iPartStart, int iPartEnd, Params& params )
//{
//    for ( int iDim = 0 ; iDim < Position.size() ; iDim++ ) {
//        double dx2 = params.cell_length[iDim]*params.cell_length[iDim];
//        for (int iPart = iPartStart ; iPart < iPartEnd ; iPart++ ) {
//            if ( dist(iPart,iDim) > dx2 ) {
//                ERROR( "Too large displacment for particle : " << iPart << "\t: " << (*this)(iPart) );
//                return false;
//            }
//        }
//    }
//
//}

Particle Particles::operator()(unsigned int iPart)
{
    return  Particle( *this, iPart);
}
