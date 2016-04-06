#ifndef IONIZATION_H
#define IONIZATION_H

#include <map>

#include "Tools.h"
#include "Params.h"
#include "Field.h"
#include "Particles.h"


//! Class Ionization: generic class allowing to define Ionization physics
class Ionization
{

public:
    //! Constructor for Ionization
    Ionization(Params& params, Species * species);
    virtual ~Ionization();

    //! Overloading of () operator
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart) = 0;

    //! Overloading of () operator
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Jion) = 0;

    Particles new_electrons;

protected:
    std::vector<double> Potential;
    std::vector<double> Azimuthal_quantum_number;

    double eV_to_au;
    double EC_to_au;
    double au_to_w0;

    double referenceAngularFrequency_SI;
    double dt;
    unsigned int nDim_field;
    unsigned int nDim_particle;
    unsigned int atomic_number_;
    unsigned int ionized_species_mass;

private:


};

#endif
