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

PartWall::PartWall(double pos, unsigned short dir) :
position(pos),
direction(dir) {}

// Reads the input file and creates the ParWall objects accordingly
vector<PartWall*> PartWall::create(Params& params, Patch* patch)
{
    vector<PartWall*> vecPartWall;
    
    if (patch->isMaster()) MESSAGE(1,"Adding particle walls:");
    
    // Loop over each wall component and parse info
    unsigned int numpartwall=PyTools::nComponents("PartWall");
    for (unsigned int iwall = 0; iwall < numpartwall; iwall++) {
    
        // Extract the direction of the wall
        short direction = -1;
        double position;
        string dirstring;
        if (PyTools::extract("x",position,"PartWall",iwall)) {
            direction=0;
            dirstring="x";
        }
        if (PyTools::extract("y",position,"PartWall",iwall)) {
            if (direction>=0)
                ERROR("For PartWall #" << iwall << ", cannot have several locations (x, y or z)");
            if (params.nDim_particle < 2)
                ERROR("PartWall #" << iwall << " cannot have y-location in 1D");
            direction=1;
            dirstring="y";
        }
        if (PyTools::extract("z",position,"PartWall",iwall)) {
            if (direction>=0)
                ERROR("For PartWall #" << iwall << ", cannot have several locations (x, y or z)");
            if (params.nDim_particle < 3)
                ERROR("PartWall #" << iwall << " cannot have z-location y in 1D or 2D");
            direction=2;
            dirstring="z";
        }
        if( direction < 0 ) {
            ERROR("PartWall #" << iwall << " must have one location (x, y or z)");
        }
        
        // Find out wether this proc has the wall or not
        if ( position > patch->getDomainLocalMin(direction) && position < patch->getDomainLocalMax(direction)) {
            
            // Ewtract the kind of wall
            string kind("");
            PyTools::extract("kind",kind,"PartWall",iwall);
            if (kind.empty() || (kind!="refl" && kind!="supp" && kind!="stop" && kind!="thermalize")) {
                ERROR("For PartWall #" << iwall << ", `kind` must be one of refl, supp, stop, thermalize");
            }
            
            // Create new wall
            PartWall * tmpWall = new PartWall(position, direction);
            
            // Define the "wall" function pointer
            bool thermCond = false;
            if (kind == "refl" ) {
                tmpWall->wall = &refl_particle;
            } else if (kind == "supp" ) {
                tmpWall->wall = &supp_particle;
            } else if (kind == "stop" ) {
                tmpWall->wall = &stop_particle;
            } else if (kind == "thermalize" ) {
                thermCond = true;
                tmpWall->wall = &thermalize_particle;
            }
            
            if (patch->isMaster()) MESSAGE(2,"Adding a wall at " << position << " in " <<  dirstring << " direction kind:" << kind << (thermCond? "thermCond" : ""));
            vecPartWall.push_back(tmpWall);
        }
        
    }
    if (!vecPartWall.size()) {
        if (patch->isMaster()) MESSAGE(2,"Nothing to do");
    }
    
    return vecPartWall;
}

int PartWall::apply( Particles &particles, int ipart, Species * species, double &nrj_iPart) {
    if( (position-particles.position_old(direction, ipart))
       *(position-particles.position    (direction, ipart))<0.) {
        return (*wall)( particles, ipart, direction, 2.*position, species, nrj_iPart );
    } else {
        return 1;
    }
}



