#ifndef TemperatureFactory_H
#define TemperatureFactory_H

#include "TemperatureProfile.h"
#include "TemperatureProfile1D.h"
#include "PicParams.h"
#include "SmileiMPI.h"

#include "Tools.h"
#include <cmath>

//! staic creator of plasma Temperature profiles (one per species)
class TemperatureFactory {
public:
    static TemperatureProfile* create(PicParams& params, unsigned int speciesNumber, int direction) {
       
        TemperatureProfile* TemperatureProfile = NULL;
        // ---------------
        // 1d3v simulation
        // ---------------
        if (params.geometry == "1d3v") {
            switch (direction) {
                case 0:
                    TemperatureProfile = new TemperatureProfile1D(params.species_param[speciesNumber].temp_x_profile);
                    break;
                case 1:
                    TemperatureProfile = new TemperatureProfile1D(params.species_param[speciesNumber].temp_y_profile);
                    break;
                case 2:
                    TemperatureProfile = new TemperatureProfile1D(params.species_param[speciesNumber].temp_z_profile);
                    break;
                default:
                    break;
            }
        }
        
        else {
            ERROR( "Unsupported geometry : " << params.geometry);
        }
        return TemperatureProfile;        
    }
};

#endif

