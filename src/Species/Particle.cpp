
#include "Particle.h"

#include "Particles.h"

using namespace std;

Particle::Particle( Particles &parts, int iPart )
{
    Position.resize( parts.Position.size() );
    Position_old.resize( parts.Position.size() );
    Momentum.resize( 3 );
    for( unsigned int iDim = 0 ; iDim < parts.Position.size() ; iDim++ ) {
        Position[iDim]     = parts.position( iDim, iPart );
        Position_old[iDim] = parts.position_old( iDim, iPart );
    }
    for( int iDim = 0 ; iDim < 3 ; iDim++ ) {
        Momentum[iDim]     = parts.momentum( iDim, iPart );
    }
    Weight = parts.weight( iPart );
    Charge = parts.charge( iPart );
    
    if( parts.Chi.size() ) {
        Chi = parts.chi( iPart );
    }
    if( parts.Tau.size() ) {
        Tau = parts.tau( iPart );
    }
    if( parts.Id.size() ) {
        Id  = parts.id( iPart );
    }
};


ostream &operator << ( ostream &out, const Particle &particle )
{
    for( unsigned int i=0; i<particle.Position.size(); i++ ) {
        out << particle.Position[i] << " ";
        out << particle.Position_old[i] << " ";
    }
    for( unsigned int i=0; i<3; i++ ) {
        out << particle.Momentum[i] << " ";
    }
    
    out << particle.Weight << " ";
    out << particle.Charge << " ";
    
    if( 0 ) {
        out << particle.Id << " " ;
    }
    if( 0 ) {
        out << particle.Chi << " " ;
    }
    if( 0 ) {
        out << particle.Tau << " " ;
    }
    return ( out );
}
