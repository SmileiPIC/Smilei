
#ifndef PROJECTOR1D4ORDER_H
#define PROJECTOR1D4ORDER_H

#include "Projector1D.h"

class Projector1D4Order : public Projector1D {
public:
    Projector1D4Order(PicParams*, SmileiMPI* smpi);
    ~Projector1D4Order();
    void operator() (ElectroMagn* champs, Particles &particles, int ipart, double gf);
    void operator() (Field* rho, Particles &particles, int ipart);
    void operator() (double* Jx, double* Jy, double* Jz, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim0);
    void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf);
    void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion);

private:
    double dx_ov_dt;
    double dble_1_ov_384 ;
    double dble_1_ov_48 ;
    double dble_1_ov_16 ;
    double dble_1_ov_12 ;
    double dble_1_ov_24 ;
    double dble_19_ov_96 ;
    double dble_11_ov_24 ;
    double dble_1_ov_4 ;
    double dble_1_ov_6 ;
    double dble_115_ov_192 ;
    double dble_5_ov_8 ;
};

#endif

