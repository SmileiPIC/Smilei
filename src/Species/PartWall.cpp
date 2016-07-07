#include "PartWall.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include <cmath>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "Patch.h"
#include "Tools.h"

using namespace std;

PartWall::PartWall(double pos, unsigned short dir, string kind) :
    position(pos),
    direction(dir)
{
    // Define the "wall" function pointer
    if (kind == "refl" ) {
        wall = &refl_particle;
    } else if (kind == "supp" ) {
        wall = &supp_particle;
    } else if (kind == "stop" ) {
        wall = &stop_particle;
    } else if (kind == "thermalize" ) {
        wall = &thermalize_particle;
    }
}

int PartWall::apply( Particles &particles, int ipart, Species * species, double &nrj_iPart) {
    if( (position-particles.position_old(direction, ipart))
       *(position-particles.position    (direction, ipart))<0.) {
        return (*wall)( particles, ipart, direction, 2.*position, species, nrj_iPart );
    } else {
        return 1;
    }
}




// Reads the input file and creates the ParWall objects accordingly
PartWalls::PartWalls(Params& params, Patch* patch)
{
    if (patch->isMaster()) MESSAGE(1,"Adding particle walls:");
    
    resize(0);
    unsigned int numpartwall=PyTools::nComponents("PartWall");
    direction.resize(numpartwall);
    position .resize(numpartwall);
    kind     .resize(numpartwall);
    
    // Loop over each wall component and parse info
    for (unsigned int iwall = 0; iwall < numpartwall; iwall++) {
        
        // Extract the direction of the wall
        direction[iwall] = -1;
        string dirstring;
        if (PyTools::extract("x",position[iwall],"PartWall",iwall)) {
            direction[iwall]=0;
            dirstring="x";
        }
        if (PyTools::extract("y",position[iwall],"PartWall",iwall)) {
            if (direction[iwall]>=0)
                ERROR("For PartWall #" << iwall << ", cannot have several locations (x, y or z)");
            if (params.nDim_particle < 2)
                ERROR("PartWall #" << iwall << " cannot have y-location in 1D");
            direction[iwall]=1;
            dirstring="y";
        }
        if (PyTools::extract("z",position[iwall],"PartWall",iwall)) {
            if (direction[iwall]>=0)
                ERROR("For PartWall #" << iwall << ", cannot have several locations (x, y or z)");
            if (params.nDim_particle < 3)
                ERROR("PartWall #" << iwall << " cannot have z-location y in 1D or 2D");
            direction[iwall]=2;
            dirstring="z";
        }
        if( direction[iwall] < 0 ) {
            ERROR("PartWall #" << iwall << " must have one location (x, y or z)");
        }
        
        // Ewtract the kind of wall
        PyTools::extract("kind",kind[iwall],"PartWall",iwall);
        if (kind[iwall].empty() || (kind[iwall]!="refl" && kind[iwall]!="supp" && kind[iwall]!="stop" && kind[iwall]!="thermalize")) {
            ERROR("For PartWall #" << iwall << ", `kind` must be one of refl, supp, stop, thermalize");
        }
        
        // Find out wether this proc has the wall or not
        if ( position[iwall] >= patch->getDomainLocalMin(direction[iwall])
          && position[iwall] <= patch->getDomainLocalMax(direction[iwall])) {
            push_back( new PartWall(position[iwall], direction[iwall], kind[iwall]) );
        }
        
        // Create new wall
        MESSAGE(2,"Adding a wall at "<<dirstring<<" = "<< position[iwall] << ", kind:" << kind[iwall] << (kind[iwall]=="thermalize" ? " thermCond" : ""));
    }
    
    if (!direction.size()) {
        if (patch->isMaster()) MESSAGE(2,"Nothing to do");
    }
}

// Clones an existing vector of partWalls
PartWalls::PartWalls(PartWalls* partWalls, Patch* patch)
{
    resize(0);
    direction = partWalls->direction;
    position  = partWalls->position ;
    kind      = partWalls->kind     ;
    
    for (unsigned int iwall = 0; iwall < direction.size(); iwall++) {
        if ( position[iwall] >= patch->getDomainLocalMin(direction[iwall])
          && position[iwall] <= patch->getDomainLocalMax(direction[iwall])) {
            push_back( new PartWall(position[iwall], direction[iwall], kind[iwall]) );
        }
    }
}


// Destructor
PartWalls::~PartWalls()
{
    int n=size();
    for( int i=0; i<n; i++ ) delete vecPartWall[i];
    clear();
}


