#ifndef PROJECTOR2D4ORDER_H
#define PROJECTOR2D4ORDER_H

#include "Projector2D.h"


class Projector2D4Order : public Projector2D {
public:
    Projector2D4Order(PicParams*, SmileiMPI* smpi);
    ~Projector2D4Order();

    void operator() (ElectroMagn* champs, Particles &particles, int ipart, double gf);
    void operator() (double* Jx, double* Jy, double* Jz, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim0);
    void operator() (Field* rho, Particles &particles, int ipart);
    void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf);
    void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion);

private:
    double one_third;
};

#endif

