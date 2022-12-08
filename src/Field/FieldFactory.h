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
    static Field* create1D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in, Params& params, bool allocate=true ) {
        if( allocate ) {
            if ( !params.is_pxr ) { // FDTD, staggered grids
                return new Field1D( dims, mainDim, isPrimal, name_in );
            } else { // if PXR, same size for all grids
                return new Field1D( dims, name_in );
            }
        } else {
            return new Field1D( name_in, dims );
        }
    }
    static Field* create2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in, Params& params, bool allocate=true ) {
        if( allocate ) {
            if ( !params.is_pxr ) { // FDTD, staggered grids
                return new Field2D( dims, mainDim, isPrimal, name_in );
            } else { // if PXR, same size for all grids
                return new Field2D( dims, name_in );
            }
        } else {
            return new Field2D( name_in, dims );
        }
    }
    static Field* create3D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in, Params& params, bool allocate=true ) {
        if( allocate ) {
            if ( !params.is_pxr ) { // FDTD, staggered grids
                return new Field3D( dims, mainDim, isPrimal, name_in );
            } else { // if PXR, same size for all grids
                return new Field3D( dims, name_in );
            }
        } else {
            return new Field3D( name_in, dims );
        }
    }
    static cField2D* createAM( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in, Params& params, bool allocate=true ) {
        if( allocate ) {
            if ( !params.is_pxr ) { // FDTD, staggered grids
                return new cField2D( dims, mainDim, isPrimal, name_in );
            } else { // if PXR, same size for all grids
                return new cField2D( dims, name_in );
            }
        } else {
            return new cField2D( name_in, dims );
        }
    }
    
    static Field* create1D( std::vector<unsigned int> dims, std::string name_in, bool allocate=true ) {
        if( allocate ) {
            return new Field1D( dims, name_in );
        } else {
            return new Field1D( name_in, dims );
        }
    }
    static Field* create2D( std::vector<unsigned int> dims, std::string name_in, bool allocate=true ) {
        if( allocate ) {
            return new Field2D( dims, name_in );
        } else {
            return new Field2D( name_in, dims );
        }
    }
    static Field* create3D( std::vector<unsigned int> dims, std::string name_in, bool allocate=true ) {
        if( allocate ) {
            return new Field3D( dims, name_in );
        } else {
            return new Field3D( name_in, dims );
        }
    }
    static cField2D* createAM( std::vector<unsigned int> dims, std::string name_in, bool allocate=true ) {
        if( allocate ) {
            return new cField2D( dims, name_in );
        } else {
            return new cField2D( name_in, dims );
        }
    }



};

#endif
