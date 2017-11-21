
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
    LaserEnvelope( Params& params, Patch* patch );
    LaserEnvelope( LaserEnvelope *envelope, Patch* patch );
    virtual void initEnvelope( Patch* patch ) = 0;
    virtual ~LaserEnvelope();
    virtual void compute(ElectroMagn* EMfields) = 0;

    Profile *profile_;
    const std::vector<double> cell_length;
        
    Field* A_;
    Field* A0_;
};


// Class for envelope
class LaserEnvelope3D : public LaserEnvelope {
public:
    LaserEnvelope3D( Params& params, Patch* patch );
    LaserEnvelope3D( LaserEnvelope *envelope, Patch* patch );
    void initEnvelope( Patch* patch ) override final;
    ~LaserEnvelope3D();
     void compute(ElectroMagn* EMfields) override final;
};

#endif
