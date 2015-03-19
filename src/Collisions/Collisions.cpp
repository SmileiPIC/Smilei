#include "Collisions.h"

#include <cmath>
#include <ostream>

using namespace std;


// Constructor
Collisions::Collisions(PicParams& param, unsigned int ncol, vector<unsigned int> sgroup1, vector<unsigned int>sgroup2, double clog, bool intra)
{

    n_collisions     = ncol;
    species_group1   = sgroup1;
    species_group2   = sgroup2;
    coulomb_log      = clog;
    intra_collisions = intra;
    
}


// Declare static variables here
bool               Collisions::debye_length_required;
vector<double>     Collisions::debye_length_squared;

