#ifndef DensityFactory_H
#define DensityFactory_H

#include "DensityProfile.h"
#include "DensityProfile1D.h"
#include "DensityProfile2D.h"

#include "PicParams.h"
#include "SmileiMPI.h"

#include "Tools.h"

//! staic creator of plasma density profiles (one per species)
class DensityFactory {
public:
    static DensityProfile* create(PicParams& params, unsigned int speciesNumber) {
        DensityProfile* densityProfile = NULL;
        // ---------------
        // 1d3v simulation
        // ---------------
        if (params.geometry == "1d3v") {
            densityProfile = new DensityProfile1D(params.species_param[speciesNumber]);
        }
        // ---------------
        // 2d3v simulation
        // ---------------
        else if (params.geometry == "2d3v") {
            densityProfile = new DensityProfile2D(params.species_param[speciesNumber]);
        }
        else {
            ERROR( "Unsupported geometry : " << params.geometry);
        }
        return densityProfile;        
    }
};

#endif

