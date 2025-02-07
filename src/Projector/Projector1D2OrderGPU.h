#ifndef SMILEI_PROJECTOR_PROJECTOR1D2ORDERGPU_H
#define SMILEI_PROJECTOR_PROJECTOR1D2ORDERGPU_H

#include "Projector1D.h"


class Projector1D2OrderGPU : public Projector1D
{
public:
    Projector1D2OrderGPU( Params &parameters, Patch *a_patch );
    ~Projector1D2OrderGPU();

    /// For initialization and diags, doesn't use the standard scheme
    void basic( double      *rhoj,
                Particles   &particles,
                unsigned int ipart,
                unsigned int type,
                int bin_shift = 0 ) override;
    /// Projection wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields,
                                    Particles   &particles,
                                    SmileiMPI   *smpi,
                                    int          istart,
                                    int          iend,
                                    int          ithread,
                                    bool         diag_flag,
                                    bool         is_spectral,
                                    int          ispec,
                                    int          icell     = 0,
                                    int          ipart_ref = 0 ) override;

    void susceptibility( ElectroMagn *EMfields,
                         Particles   &particles,
                         double       species_mass,
                         SmileiMPI   *smpi,
                         int          istart,
                         int          iend,
                         int          ithread,
                         int          icell     = 0,
                         int          ipart_ref = 0 ) override;
                         
    void ionizationCurrents( Field      *Jx,
                             Field      *Jy,
                             Field      *Jz,
                             Particles  &particles,
                             int         ipart,
                             LocalFields Jion ) override;


protected:
    double dts2_;
    double dts4_;
    int    not_spectral_;
    unsigned int x_dimension_bin_count_;
};

#endif
