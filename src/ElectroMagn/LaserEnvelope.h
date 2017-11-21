
#ifndef LaserENVELOPE_H
#define LaserENVELOPE_H

#include <vector>
class Params;
class Field;
class ElectroMagn;

// Class for envelope
class LaserEnvelope {
public:
    LaserEnvelope(Params& params);
    LaserEnvelope(LaserEnvelope *envelope);
    virtual ~LaserEnvelope();
    virtual void compute(ElectroMagn* EMfields) = 0;
    
    Field* A_;
    Field* A0_;
};

// Class for envelope
class LaserEnvelope3D : public LaserEnvelope {
public:
    LaserEnvelope3D(Params& params);
    LaserEnvelope3D(LaserEnvelope *envelope);
    ~LaserEnvelope3D();
     void compute(ElectroMagn* EMfields) override final;
};



#endif
