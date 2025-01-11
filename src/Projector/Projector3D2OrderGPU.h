#ifndef SMILEI_PROJECTOR_PROJECTOR3D2ORDERGPU_H
#define SMILEI_PROJECTOR_PROJECTOR3D2ORDERGPU_H

#include "Projector3D.h"

/// Particle to grid projector (~~dual to the grid to particle the interpolator
/// does)
///
/// NOTE: we could have inherited from Projector2D2Order but the interface is final for most of the member functions
///
class Projector3D2OrderGPU : public Projector3D
{
public:
    Projector3D2OrderGPU( Params &parameters, Patch *a_patch );
    ~Projector3D2OrderGPU();

    /// For initialization and diags, doesn't use the standard scheme
    ///
    void basic( double      *rhoj,
                Particles   &particles,
                unsigned int ipart,
                unsigned int type,
                int bin_shift = 0 ) override;

    /// Project global current densities (ionize)
    ///
    void ionizationCurrents( Field      *Jx,
                             Field      *Jy,
                             Field      *Jz,
                             Particles  &particles,
                             int         ipart,
                             LocalFields Jion ) override;

    /// Projection wrapper
    ///
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

    /// Project susceptibility, used as source term in envelope equation
    ///
    void susceptibility( ElectroMagn *EMfields,
                         Particles   &particles,
                         double       species_mass,
                         SmileiMPI   *smpi,
                         int          istart,
                         int          iend,
                         int          ithread,
                         int          icell     = 0,
                         int          ipart_ref = 0 ) override;

protected:
    double dt;
    double dts2;
    double dts4;
    int    not_spectral_;
    unsigned int x_dimension_bin_count_;
    unsigned int y_dimension_bin_count_;
    unsigned int z_dimension_bin_count_;
};

#endif
