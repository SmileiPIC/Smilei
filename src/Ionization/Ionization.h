#ifndef IONIZATION_H
#define IONIZATION_H

#include <functional>
#include <map>

#include "Field.h"
#include "Params.h"
#include "Particles.h"
#include "Patch.h"
#include "Projector.h"
#include "Tools.h"

using namespace std;

//! Class Ionization: generic class allowing to define Ionization physics
class Ionization
{
public:
    //! Constructor for Ionization
    Ionization(Params &params, Species *species);
    virtual ~Ionization();

    //! Overloading of () operator
    virtual void operator()(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int = 0) {};
    //! method for envelope ionization
    virtual void envelopeIonization( Particles *, unsigned int, unsigned int, std::vector<double> *, std::vector<double> *, std::vector<double> *, std::vector<double> *, Patch *, Projector *, int = 0, int = 0 ){};

    Particles new_electrons;
    
    //! Whether the initial charge (of the atom that was ionized) should be saved
    bool save_ion_charge_ = false;
    //! Temporarily contains the initial charge of the atom that was ionized
    std::vector<short> ion_charge_;

protected:
    double eV_to_au;
    double au_to_mec2;
    double EC_to_au;
    double au_to_w0;

    double reference_angular_frequency_SI;
    double dt;
    double invdt;
    unsigned int nDim_field;
    unsigned int nDim_particle;
    double ionized_species_invmass;

private:


};

#endif
