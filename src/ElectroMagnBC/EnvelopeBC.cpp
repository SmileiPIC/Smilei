
#include <cstdlib>
#include <iostream>
#include <string>

#include "EnvelopeBC.h"

#include "Params.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;

// Constructor for EnvelopeBC
EnvelopeBC::EnvelopeBC( Params &params, Patch *patch, unsigned int i_boundary ) :
    i_boundary_( i_boundary )
{

    // time step
    dt = params.timestep;
}

// Destructor for EnvelopeBC
EnvelopeBC::~EnvelopeBC()
{
}
