
#ifndef LaserENVELOPE_H
#define LaserENVELOPE_H

#include <vector>

#include "Params.h"
#include "Patch.h"

class Field;
class ElectroMagn;
class Profile;
class EnvelopeBC;
class SimWindow;

// Class for envelope
class LaserEnvelope {
public:
    LaserEnvelope( Params& params, Patch* patch, ElectroMagn* EMfields ); // Main constructor
    LaserEnvelope( LaserEnvelope *envelope, Patch* patch, ElectroMagn* EMfields, Params& params, unsigned int n_moved ); // Cloning constructor
    virtual void initEnvelope( Patch* patch , ElectroMagn* EMfields) = 0;
    virtual ~LaserEnvelope();
    virtual void compute(ElectroMagn* EMfields) = 0;
    virtual void compute_Phi(ElectroMagn* EMfields) = 0;
    virtual void compute_gradient_Phi(ElectroMagn* EMfields) = 0;
    void boundaryConditions(int itime, double time_dual, Patch* patch, Params &params, SimWindow* simWindow);
    
    Profile *profile_;
    const std::vector<double> cell_length;
    const double timestep;
    
    Field* A_;
    Field* A0_;

    Field* Phi_;
    Field* GradPhix_;
    Field* GradPhiy_;
    Field* GradPhiz_;
    Field* Phiold_;
    Field* GradPhixold_;
    Field* GradPhiyold_;
    Field* GradPhizold_;

    //! Vector of boundary-condition per side for the envelope field
    std::vector<EnvelopeBC*> EnvBoundCond;
    //EnvBoundCond = EnvelopeBC_Factory::create(params, patch);
};


// Class for envelope
class LaserEnvelope3D : public LaserEnvelope {
public:
    LaserEnvelope3D( Params& params, Patch* patch, ElectroMagn* EMfields );
    LaserEnvelope3D( LaserEnvelope *envelope, Patch* patch, ElectroMagn* EMfields, Params& params, unsigned int n_moved );
    void initEnvelope( Patch* patch,ElectroMagn* EMfields ) override final;
    ~LaserEnvelope3D();
     void compute(ElectroMagn* EMfields) override final;
     void compute_Phi(ElectroMagn* EMfields) override final;
     void compute_gradient_Phi(ElectroMagn* EMfields) override final;
};

#endif
