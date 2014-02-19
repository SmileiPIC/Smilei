
#ifndef PROJECTOR1D2ORDER_H
#define PROJECTOR1D2ORDER_H

#include "Projector1D.h"

class Projector1D2Order : public Projector1D {
public:
    Projector1D2Order(PicParams*, SmileiMPI* smpi);
    ~Projector1D2Order();
    void operator() (ElectroMagn* champs, Particles &particles, int ipart, double gf);
    void operator() (Field* rho, Particles &particles, int ipart);
    void operator() (double* Jx, double* Jy, double* Jz, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim0);
    void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf);
    void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion);
private:
    double dx_ov_dt;
};

#endif

