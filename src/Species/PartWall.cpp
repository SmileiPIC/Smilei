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

PartWall::PartWall( Params& params, SmileiMPI* smpi ) :
is_here(false)
{
    // number of dimensions for the particle
    //!\todo (MG to JD) isn't it always 3?
    nDim_particle = params.nDim_particle;
    
    
    if (PyTools::extract("x",pos,"PartWall")) {
        direction=0;
    } else if (PyTools::extract("y",pos,"PartWall")) {
        direction=1;
        if (params.nDim_particle < 2) {
            ERROR("position of wall can't be in y direction");
        }
    } else if (PyTools::extract("z",pos,"PartWall")) {
        direction=2;
        if (params.nDim_particle < 3) {
            ERROR("position of wall can't be in z direction");
        }
    } else {
        ERROR("position of wall must be one of x, y, z");
    }
    
    is_here = (pos > smpi->getDomainLocalMin(direction) && pos < smpi->getDomainLocalMax(direction));

    string kind("");
    PyTools::extract("kind",kind,"PartWall");
    if (kind.empty()) {
        ERROR("kind of wall must be one of refl, supp, stop, thermalize");
    }
    
    bool thermCond = false;
    
    if (direction==0) {
        if (kind == "refl" ) {
            wall_x = &refl_particle;
        }
        else if (kind == "supp" ) {
            wall_x = &supp_particle;
        }
        else if (kind == "stop" ) {
            wall_x = &stop_particle;
        }
        else if (kind == "thermalize" ) {
            thermCond = true;
            wall_x = &thermalize_particle;
        }
    } else if (direction==1) {
        if (kind == "refl" ) {
            wall_y = &refl_particle;
        }
        else if (kind == "supp" ) {
            wall_y = &supp_particle;
        }
        else if (kind == "stop" ) {
            wall_y = &stop_particle;
        }
        else if (kind == "thermalize" ) {
            thermCond = true;
            wall_y = &thermalize_particle;
        }
    } else if (direction==2) {
        if (kind == "refl" ) {
            wall_z = &refl_particle;
        }
        else if (kind == "supp" ) {
            wall_z = &supp_particle;
        }
        else if (kind == "stop" ) {
            wall_z = &stop_particle;
        }
        else if (kind == "thermalize" ) {
            thermCond = true;
            wall_z = &thermalize_particle;
        }
    }
}


