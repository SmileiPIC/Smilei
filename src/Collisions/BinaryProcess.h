#ifndef BINARYPROCESS_H
#define BINARYPROCESS_H

#include <vector>

#include "Tools.h"
#include "Random.h"
#include "Diagnostic.h"
#include "BinaryProcessData.h"

class BinaryProcess
{
public:
    BinaryProcess() {};
    
    virtual ~BinaryProcess() {};
    
    virtual void prepare() = 0;
    virtual void apply( Random *random, BinaryProcessData &D ) = 0;
    virtual void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool intra, std::vector<unsigned int> sg1, std::vector<unsigned int> sg2, int itime ) = 0;
    virtual std::string name() = 0;
};

#endif
