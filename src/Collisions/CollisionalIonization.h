
#ifndef COLLISIONALIONIZATION_H
#define COLLISIONALIONIZATION_H

#include <vector>

#include "Tools.h"
#include "Species.h"

class CollisionalIonization
{

public:
    //! Constructor
    CollisionalIonization(int, double,SmileiMPI*);
    virtual ~CollisionalIonization() {};
    
    //! Gets the k-th binding energy of any neutral or ionized atom with atomic number Z and charge Zstar
    double binding_energy(int Zstar, int k);
    
    //! Coefficients used for interpolating the energy over a given initial list
    static const double a1, a2, npointsm1;
    static const int npoints;
    
    //! Methods to prepare the ionization
    inline void prepare1(int Z_firstgroup) {
        electronFirst = Z_firstgroup==0 ? true : false;
        ne = 0.; ni = 0.; nei = 0.;
        };
    virtual void prepare2(Particles *p1, int i1, Particles *p2, int i2, bool);
    virtual void prepare3(double, int);
    //! Method to apply the ionization
    virtual void apply(double, double, Particles *p1, int i1, Particles *p2, int i2);
    
    //! Table of integrated cross-section
    std::vector<std::vector<double> > crossSection;
    //! Table of average secondary electron energy
    std::vector<std::vector<double> > transferredEnergy;
    //! Table of average incident electron energy lost
    std::vector<std::vector<double> > lostEnergy;
    
    //! Temporary stuff before patches arrive
    Particles new_electrons;
    virtual void finish(Species *s1, Species *s2, Params&);
    SmileiMPI* smpi;
    
private:
    
    //! Atomic number
    int atomic_number;
    
    //! Table of first ionization energies of ions
    static const std::vector<std::vector<double> > ionizationEnergy;
    //! Table of binding energies of all electrons in neutral atoms
    static const std::vector<std::vector<double> > bindingEnergy;
    
    //! True if first group of species is the electron
    bool electronFirst;
    //! Ionizing electron density
    double ne;
    //! Ionizing "hybrid" density
    double nei;
    //! Ion density
    double ni;
    //! Coefficient for the ionization frequency
    double coeff;
    
    //! Method called by ::apply to calculate the ionization, being sure that electrons are the first species
    void calculate(double, double, Particles *pe, int ie, Particles *pi, int ii);
    
    //! Quantities used during computation
    int Zstar; // ion charge
    
};

//! Class to make empty ionization objects
class CollisionalNoIonization : public CollisionalIonization
{
public:
    CollisionalNoIonization() : CollisionalIonization(0,1.,NULL) {};
    ~CollisionalNoIonization(){};
    
    void prepare2(Particles*, int, Particles*, int, bool){};
    void prepare3(double, int){};
    void apply(double, double, Particles*, int, Particles*, int){};
    
    //! Temporary stuff before patches arrive
    void finish(Species*, Species*, Params&) {};
};

#endif
