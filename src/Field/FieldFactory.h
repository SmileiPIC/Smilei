#ifndef FIELDFACTORY_H
#define FIELDFACTORY_H

#include <vector>
#include <string>

#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"

#include "Params.h"

class FieldFactory {
public:
    static Field* create( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in, Params& params ) {
        if (params.geometry == "1Dcartesian")
            return new Field1D( dims, mainDim, isPrimal, name_in );
        else if (params.geometry == "2Dcartesian")
            return new Field2D( dims, mainDim, isPrimal, name_in );
        else if (params.geometry == "3Dcartesian")
            return new Field3D( dims, mainDim, isPrimal, name_in );
        else
            return nullptr;
    }

};

#endif
