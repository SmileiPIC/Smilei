#ifndef PROJECTOR_H
#define PROJECTOR_H

#include <complex>

#include "Params.h"
#include "Field.h"

class PicParams;
class Patch;

class ElectroMagn;
class ElectroMagnAM;
class Field;
class Particles;

//----------------------------------------------------------------------------------------------------------------------
//! class Projector: contains the virtual operators used during the current projection
//----------------------------------------------------------------------------------------------------------------------
class Projector
{

public:
    //! Creator for the Projector
    Projector( Params &, Patch * );
    virtual ~Projector() {};
    virtual void mv_win( unsigned int shift ) = 0;
    virtual void setMvWinLimits( unsigned int shift ) = 0;

    //! Project global current charge (EMfields->rho_ , J), for initialization and diags
    virtual void basic( double *, Particles &, unsigned int, unsigned int, int = 0 ) {};
    virtual void basicForComplex( std::complex<double> *, Particles &, unsigned int, unsigned int, int ) {};
    // basic rho projector for tasks
    virtual void basicForComplexOnBuffer( std::complex<double> *, Particles &, unsigned int, unsigned int, int, int, int ) {};

    //! Apply boundary conditions on axis for Rho and J in AM geometry
    virtual void axisBC(ElectroMagnAM *, bool) {};

    //! Apply boundary conditions on axis for Env_Chi in AM geometry
    virtual void axisBCEnvChi( double * ) {};

    //! Project global current densities if Ionization in Species::dynamics,
    virtual void ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion ) = 0;

    //! Project global current densities if Ionization in Species::dynamics,
    virtual void ionizationCurrentsForTasks( double *, double *, double *, Particles &, int, LocalFields, int ) {};

    //!Wrapper
    virtual void currentsAndDensityWrapper( ElectroMagn *, Particles &, SmileiMPI *, int, int, int, bool, bool, int, int = 0, int = 0 ) = 0;

    //!Wrapper for tasks
    virtual void currentsAndDensityWrapperOnBuffers( double *b_Jx, double *b_Jy, double *b_Jz, double *b_rho, int bin_width, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell = 0, int ipart_ref = 0 ) = 0;

    //!Wrapper for projection on AM buffers
    virtual void currentsAndDensityWrapperOnAMBuffers( ElectroMagn *, std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *, int, int, Particles &, SmileiMPI *, int, int, int, bool, int = 0 ) {};
    
    virtual void susceptibility( ElectroMagn *, Particles &, double , SmileiMPI *, int, int, int, int = 0, int = 0 )
    {
        ERROR( "Envelope not implemented with this geometry and this order" );
    };

    virtual void susceptibilityOnBuffer( ElectroMagn *, double *, int, int, Particles &, double, SmileiMPI *, int, int, int, int = 0, int = 0 )
    {
        ERROR( "Envelope with tasks not implemented with this geometry and this order" );
    };

    
protected:
    double inv_cell_volume;
};

#endif
