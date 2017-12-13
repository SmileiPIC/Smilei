
#ifndef LaserENVELOPE_H
#define LaserENVELOPE_H

#include <vector>

#include "Params.h"
#include "Patch.h"

class Field;
class ElectroMagn;
class Profile;


// Class for envelope
class LaserEnvelope {
public:
    LaserEnvelope( Params& params, Patch* patch, ElectroMagn* EMfields );
    LaserEnvelope( LaserEnvelope *envelope, Patch* patch, ElectroMagn* EMfields );
    virtual void initEnvelope( Patch* patch , ElectroMagn* EMfields) = 0;
    virtual ~LaserEnvelope();
    virtual void compute(ElectroMagn* EMfields) = 0;

    Profile *profile_;
    const std::vector<double> cell_length;
    const double timestep;
    
    Field* A_;
    Field* A0_;
  
};


// Class for envelope
class LaserEnvelope3D : public LaserEnvelope {
public:
    LaserEnvelope3D( Params& params, Patch* patch, ElectroMagn* EMfields );
    LaserEnvelope3D( LaserEnvelope *envelope, Patch* patch, ElectroMagn* EMfields );
    void initEnvelope( Patch* patch,ElectroMagn* EMfields ) override final;
    ~LaserEnvelope3D();
     void compute(ElectroMagn* EMfields) override final;
};

#endif
