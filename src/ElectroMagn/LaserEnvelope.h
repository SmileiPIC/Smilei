
#ifndef LaserENVELOPE_H
#define LaserENVELOPE_H

#include <vector>
class Field;
class ElectroMagn;

// Class for envelope
class LaserEnvelope {
public:
    LaserEnvelope(std::vector<unsigned int> dimPrim);
    virtual ~LaserEnvelope();
    virtual void compute(ElectroMagn* EMfields) = 0;
protected:
    Field* A_;
    Field* A0_;
};

// Class for envelope
class LaserEnvelope3D : public LaserEnvelope {
public:
    LaserEnvelope3D(std::vector<unsigned int> dimPrim);
    ~LaserEnvelope3D();
     void compute(ElectroMagn* EMfields) override final;
};



#endif
