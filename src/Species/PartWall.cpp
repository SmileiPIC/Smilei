#include "PartWall.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include <cmath>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

PartWall::PartWall( short direction, string kind, double position )
{
    pos = position;
    
    // Define the "apply" function pointer
    if( direction==0 ) {
        apply = &PartWall::apply_x;
    } else if( direction == 1 ) {
        apply = &PartWall::apply_y;
    } else if( direction == 2 ) {
        apply = &PartWall::apply_z;
    }
    
    // Define the "wall" function pointer
    bool thermCond = false;
    if (kind == "refl" ) {
        wall = &refl_particle;
    } else if (kind == "supp" ) {
        wall = &supp_particle;
    } else if (kind == "stop" ) {
        wall = &stop_particle;
    } else if (kind == "thermalize" ) {
        thermCond = true;
        wall = &thermalize_particle;
    }
}

// Reads the input file and creates the ParWall objects accordingly
vector<PartWall*> PartWall::create(Params& params, SmileiMPI* smpi)
{
    PartWall* partwall;
    vector<PartWall*> vecPartWall;
    
    // Loop over each wall component and parse info
    unsigned int numpartwall=PyTools::nComponents("PartWall");
    for (unsigned int iwall = 0; iwall < numpartwall; iwall++) {
        
        // Extract the direction of the wall
        short direction = -1;
        double position;
        if (PyTools::extract("x",position,"PartWall",iwall)) {
            direction=0;
        }
        if (PyTools::extract("y",position,"PartWall",iwall)) {
            if (direction>=0)
                ERROR("For PartWall #" << iwall << ", cannot have several locations (x, y or z)");
            if (params.nDim_particle < 2)
                ERROR("PartWall #" << iwall << " cannot have y-location in 1D");
            direction=1;
        }
        if (PyTools::extract("z",position,"PartWall",iwall)) {
            if (direction>=0)
                ERROR("For PartWall #" << iwall << ", cannot have several locations (x, y or z)");
            if (params.nDim_particle < 3)
                ERROR("PartWall #" << iwall << " cannot have z-location y in 1D or 2D");
            direction=2;
        }
        if( direction < 0 ) {
            ERROR("PartWall #" << iwall << " must have one location (x, y or z)");
        }
        
        // Find out wether this proc has the wall or not
        bool is_here = (  position > smpi->getDomainLocalMin(direction)
                       && position < smpi->getDomainLocalMax(direction));
        // Skip this wall if the proc does not contain it
        if (!is_here) continue;
        
        // Ewtract the kind of wall
        string kind("");
        PyTools::extract("kind",kind,"PartWall",iwall);
        if (kind.empty() || (kind!="refl" && kind!="supp" && kind!="stop" && kind!="thermalize")) {
            ERROR("For PartWall #" << iwall << ", `kind` must be one of refl, supp, stop, thermalize");
        }
        
        // Create new wall
        vecPartWall.push_back(new PartWall(direction, kind, position));
        
    }
    
    return vecPartWall;
}




