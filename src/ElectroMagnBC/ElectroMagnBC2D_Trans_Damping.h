
#ifndef ELECTROMAGNBC2D_Trans_DAMPING_H
#define ELECTROMAGNBC2D_Trans_DAMPING_H

#include "ElectroMagnBC2D.h"

class Params;
class ElectroMagn;

class ElectroMagnBC2D_Trans_Damping : public ElectroMagnBC2D
{
public:
    ElectroMagnBC2D_Trans_Damping( Params &params, Patch *patch, unsigned int i_boundary );
    ~ElectroMagnBC2D_Trans_Damping() {};
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
private:
    
    // number of dumping layers
    unsigned int ny_l;
    // Damping coefficient
    double cdamp;
    // array of coefficient per layer
    std::vector<double> coeff;
    
    
};

#endif

