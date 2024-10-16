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
    virtual void envelopeIonization(Particles *, unsigned int, unsigned int, std::vector<double> *, std::vector<double> *,
                                    std::vector<double> *, std::vector<double> *, Patch *, Projector *, int = 0, int = 0) {};

    // method for tunnel ionization using tasks
    virtual void ionizationTunnelWithTasks(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *,
                                           Projector *, int, int, double *, double *, double *, int = 0) {};
    // join the lists of electrons created through ionization when tasks are used
    void joinNewElectrons(unsigned int Nbins);

    Particles new_electrons;
    Particles *new_electrons_per_bin;

    //! Whether the initial charge (of the atom that was ionized) should be saved
    bool save_ion_charge_ = false;
    //! Temporarily contains the initial charge of the atom that was ionized
    std::vector<short> ion_charge_;
    std::vector<std::vector<short> > ion_charge_per_bin_;

   protected:
    double eV_to_au;
    double au_to_mec2;
    double EC_to_au;
    double au_to_w0;
    double one_third;

    double reference_angular_frequency_SI;
    double dt;
    double invdt;
    unsigned int nDim_field;
    unsigned int nDim_particle;
    double ionized_species_invmass;

    unsigned int atomic_number_;
    std::vector<double> Potential, Azimuthal_quantum_number;
    std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel;

   private:
};

#endif
