#ifndef FIELDFACTORY_H
#define FIELDFACTORY_H

#include <vector>
#include <string>

#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "cField2D.h"

#include "Params.h"

class FieldFactory {
public:
    static Field* create( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in, Params& params ) {
        if ( !params.is_pxr ) { // FDTD, staggered grids
            if (params.geometry == "1Dcartesian")
                return new Field1D( dims, mainDim, isPrimal, name_in );
            else if (params.geometry == "2Dcartesian")
                return new Field2D( dims, mainDim, isPrimal, name_in );
            else if (params.geometry == "3Dcartesian")
                return new Field3D( dims, mainDim, isPrimal, name_in );
            else
                return nullptr;
        }
        else { // if PXR, same size for all grids
            if (params.geometry == "1Dcartesian")
                return new Field1D( dims, name_in );
            else if (params.geometry == "2Dcartesian")
                return new Field2D( dims, name_in );
            else if (params.geometry == "3Dcartesian")
                return new Field3D( dims, name_in );
            else
                return nullptr;
        }
    }

    static Field* create( std::string name, std::vector<unsigned int> dims, Params& params ) {
        if (params.geometry == "1Dcartesian")
            return new Field1D( name, dims );
        else if (params.geometry == "2Dcartesian")
            return new Field2D( name, dims );
        else if (params.geometry == "3Dcartesian")
            return new Field3D( name, dims );
        else
            return nullptr;
    }

    static cField2D* createComplex( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in, Params& params ) {
        if ( !params.is_pxr ) { // FDTD, staggered grids
            if( params.geometry == "AMcylindrical" ) {
                //return new cField2D( dims, name_in );
                return new cField2D( dims, mainDim, isPrimal, name_in );
            }
            else
                return nullptr;
        }
        else { // if PXR, same size for all grids
            if( params.geometry == "AMcylindrical" )
                return new cField2D( dims, name_in );
            else
                return nullptr;
        }
    }

    static cField2D* createComplex( std::string name, std::vector<unsigned int> dims, Params& params ) {
        if( params.geometry == "AMcylindrical" )
            return new cField2D( name, dims );
        else
            return nullptr;
    }



};

#endif
