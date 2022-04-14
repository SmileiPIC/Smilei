#ifndef NUCLEARREACTIONPRODUCTS_H
#define NUCLEARREACTIONPRODUCTS_H

#include <vector>

class Particles;

//! Contains the data for creating products after a nuclear reaction
struct NuclearReactionProducts
{
    //! Vector of Particles objects that will contain each of the new product particles
    std::vector<Particles*> particles;
    
    //! Vector of the new momentum in the COM frame for each product
    std::vector<double> new_p_COM;
    
    //! Vectors of the sin and cos of the deflection angle for each product
    std::vector<double> sinX, cosX;
    
    //! Vectors of the sin and cos of the azimuthal angle for each product
    std::vector<double> sinPhi, cosPhi;
    
    //! Vector of the charge for each product
    std::vector<short> q;
    
};

#endif
