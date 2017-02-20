#include <iostream>
#include <cmath>

#include "Species.h"
#include "PusherBorisNR.h"
#include "Particles.h"

using namespace std;

PusherBorisNR::PusherBorisNR(Params& params, Species *species)
    : Pusher(params, species)
{
}

PusherBorisNR::~PusherBorisNR()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (Boris) scheme
***********************************************************************/

void PusherBorisNR::operator() (Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread)
{
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);

    double charge_over_mass_ ;
    double umx, umy, umz;
    double upx, upy, upz;
    double alpha;
    double Tx, Ty, Tz;
    double T2;
    double Sx, Sy, Sz;

    for (int ipart=istart ; ipart<iend; ipart++ ) {

        charge_over_mass_ = static_cast<double>(particles.charge(ipart))*one_over_mass_;
        alpha = charge_over_mass_*dts2;

        // uminus = v + q/m * dt/2 * E
        umx = particles.momentum(0, ipart) * one_over_mass_ + alpha * (*Epart)[ipart].x;
        umy = particles.momentum(1, ipart) * one_over_mass_ + alpha * (*Epart)[ipart].y;
        umz = particles.momentum(2, ipart) * one_over_mass_ + alpha * (*Epart)[ipart].z;


        // Rotation in the magnetic field

        Tx    = alpha * (*Bpart)[ipart].x;
        Ty    = alpha * (*Bpart)[ipart].y;
        Tz    = alpha * (*Bpart)[ipart].z;

        T2 = Tx*Tx + Ty*Ty + Tz*Tz;

        Sx = 2*Tx/(1.+T2);
        Sy = 2*Ty/(1.+T2);
        Sz = 2*Tz/(1.+T2);

        // uplus = uminus + uprims x S
        upx = umx + umy*Sz - umz*Sy;
        upy = umy + umz*Sx - umx*Sz;
        upz = umz + umx*Sy - umy*Sx;


        particles.momentum(0, ipart) = mass_ * (upx + alpha*(*Epart)[ipart].x);
        particles.momentum(1, ipart) = mass_ * (upy + alpha*(*Epart)[ipart].y);
        particles.momentum(2, ipart) = mass_ * (upz + alpha*(*Epart)[ipart].z);

        // Move the particle
        for ( int i = 0 ; i<nDim_ ; i++ )
            particles.position(i, ipart)     += dt*particles.momentum(i, ipart);
    }
}
