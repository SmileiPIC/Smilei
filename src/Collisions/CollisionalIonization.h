
#ifndef COLLISIONALIONIZATION_H
#define COLLISIONALIONIZATION_H

#include <vector>

#include "Tools.h"

class CollisionalIonization
{

public:
    //! Constructor
    CollisionalIonization(int, double);
    ~CollisionalIonization();
    
    //! Gets the k-th binding energy of any neutral or ionized atom with atomic number Z and charge Zstar
    double binding_energy(int Zstar, int k);
    
    //! Coefficients used for interpolating the energy over a given initial list
    static const double a1, a2;
    static const int npoints;
    
private:
    
    //! Atomic number
    int atomic_number;
    
    //! Table of first ionization energies of ions
    static const std::vector<std::vector<double> > ionizationEnergy;
    //! Table of binding energies of all electrons in neutral atoms
    static const std::vector<std::vector<double> > bindingEnergy;
    
    //! Table of integrated cross-section
    std::vector<std::vector<double> > crossSection;
    //! Table of average secondary electron energy
    std::vector<std::vector<double> > transferredEnergy;
    //! Table of average incident electron energy lost
    std::vector<std::vector<double> > lostEnergy;
    
    
    
};
    
#endif
