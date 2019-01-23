#ifndef PROJECTOR2D2ORDERV_H
#define PROJECTOR2D2ORDERV_H

#include "Projector2D.h"


class Projector2D2OrderV : public Projector2D {
public:
    Projector2D2OrderV(Params&, Patch* patch);
    ~Projector2D2OrderV();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
void currents(double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int istart, unsigned int iend, std::vector<double> *invgf, std::vector<unsigned int> &b_dim, int* iold, double *deltaold, int ipart_ref = 0);
    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
 inline void currentsAndDensity(double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold, int nparts_in_buf);

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void densityFrozen(double* rhoj, Particles &particles, unsigned int ipart, unsigned int bin) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents(Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion) override final;

    //!Wrapper
    void currentsAndDensityWrapper(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int icell, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec, int ipart_ref = 0) override final;

    // Project susceptibility
    void susceptibility(ElectroMagn* EMfields, Particles &particles, double species_mass, SmileiMPI* smpi, int istart, int iend,  int ithread, int ibin, int ipart_ref = 0) override final;

private:
};

#endif

