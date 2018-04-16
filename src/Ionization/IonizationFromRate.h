#ifndef IONIZATIONFROMRATE_H
#define IONIZATIONFROMRATE_H

#include <cmath>

#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;

//! calculate the particle FromRate ionization
class IonizationFromRate : public Ionization
{

public:
    //! Constructor for IonizationFromRate: with no input argument
    IonizationFromRate(Params& params, Species * species);
    
    //! apply the FromRate Ionization model to the species 
    void operator() (Particles*, unsigned int, unsigned int, std::vector<double>*, ElectroMagn*, Projector*) override;
    
private:
    
    int itime;
    unsigned int maximum_charge_state_;
    PyObject* ionization_rate;
    
};


#endif
