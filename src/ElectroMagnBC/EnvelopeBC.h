#ifndef ENVELOPEBC_H
#define ENVELOPEBC_H

#include <vector>

#include "Patch.h"
#include "ElectroMagn.h"

class Params;
class Patch;
class LaserEnvelope;
class Field;
class Solver;

class EnvelopeBC
{
public:
    EnvelopeBC( Params &params, Patch *patch, unsigned int i_boundary );
    virtual ~EnvelopeBC();
    void clean();

    virtual void apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;

    virtual cField* getA2Dnp1PML() { ERROR("Not using PML");return NULL;}
    virtual cField* getA2DnPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getA2Dnm1PML() { ERROR("Not using PML");return NULL;}

    virtual cField* getA3Dnp1PML() { ERROR("Not using PML");return NULL;}
    virtual cField* getA3DnPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getA3Dnm1PML() { ERROR("Not using PML");return NULL;}

    virtual cField* getG2Dnp1PML() { ERROR("Not using PML");return NULL;}
    virtual cField* getG2DnPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getG2Dnm1PML() { ERROR("Not using PML");return NULL;}

    virtual Field* getPhi2DPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getPhi3DPML() { ERROR("Not using PML");return NULL;}

    virtual cField* getu1np1xPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2np1xPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3np1xPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu1nm1xPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2nm1xPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3nm1xPML() { ERROR("Not using PML");return NULL;}

    virtual cField* getu1np1yPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2np1yPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3np1yPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu1nm1yPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2nm1yPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3nm1yPML() { ERROR("Not using PML");return NULL;}

    virtual cField* getu1np1zPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2np1zPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3np1zPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu1nm1zPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2nm1zPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3nm1zPML() { ERROR("Not using PML");return NULL;}

    virtual cField* getu1np1lPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2np1lPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3np1lPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu1nm1lPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2nm1lPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3nm1lPML() { ERROR("Not using PML");return NULL;}

    virtual cField* getu1np1rPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2np1rPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3np1rPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu1nm1rPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu2nm1rPML() { ERROR("Not using PML");return NULL;}
    virtual cField* getu3nm1rPML() { ERROR("Not using PML");return NULL;}

protected:

    //! time-step
    double dt;

    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;
    unsigned int nl_p;

    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;
    unsigned int nr_p;

    //! Number of nodes on the primal grid in the z-direction
    unsigned int nz_p;

    //! Spatial step dx for 2D3V cartesian simulations
    double dx;
    double dl;

    //! Spatial step dy for 2D3V cartesian simulations
    double dy;
    double dr;

    //! Spatial step dy for 3D3V cartesian simulations
    double dz

    //! Ratio of the time-step by the spatial-step dt/dx for 2D3V cartesian simulations
    double dt_ov_dx;

    //! Ratio of the time-step by the spatial-step dt/dy for 2D3V cartesian simulations
    double dt_ov_dy;

    //! Ratio of the time-step by the spatial-step dt/dz for 3D3V cartesian simulations
    double dt_ov_dz;

    //! Ratio of the spatial-step by the time-step dx/dt for 2D3V cartesian simulations
    double dx_ov_dt;

    //! Ratio of the spatial-step by the time-step dy/dt for 2D3V cartesian simulations
    double dy_ov_dt;

    //! Ratio of the spatial-step by the time-step dz/dt for 3D3V cartesian simulations
    double dz_ov_dt;

    // side of BC is applied 0:xmin 1:xmax 2:ymin 3:ymax 4:zmin 5:zmax
    unsigned int i_boundary_;

    Solver* pml_solver_envelope_ ;

};

#endif
