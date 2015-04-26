#ifndef VelocityFactory_H
#define VelocityFactory_H

#include "VelocityProfile.h"
#include "VelocityProfile1D.h"
#include "VelocityProfile2D.h"

#include "PicParams.h"
#include "SmileiMPI.h"

#include "Tools.h"
#include <cmath>

//! staic creator of plasma velocity profiles (one per species)
class VelocityFactory {
public:
    static VelocityProfile* create(PicParams& params, unsigned int speciesNumber, int direction) {
       
        VelocityProfile* VelocityProfile = NULL;
        // ---------------
        // 1d3v simulation
        // ---------------
        if (params.geometry == "1d3v") {
            switch (direction) {
                case 0:
                    VelocityProfile = new VelocityProfile1D(params.species_param[speciesNumber].mvel_x_profile);
                    break;
                case 1:
                    VelocityProfile = new VelocityProfile1D(params.species_param[speciesNumber].mvel_y_profile);
                    break;
                case 2:
                    VelocityProfile = new VelocityProfile1D(params.species_param[speciesNumber].mvel_z_profile);
                    break;
                default:
                    break;
            }
        }
        // ---------------
        // 2d3v simulation
        // ---------------
        else if (params.geometry == "2d3v") {
            switch (direction) {
                case 0:
                    VelocityProfile = new VelocityProfile2D(params.species_param[speciesNumber].mvel_x_profile);
                    break;
                case 1:
                    VelocityProfile = new VelocityProfile2D(params.species_param[speciesNumber].mvel_y_profile);
                    break;
                case 2:
                    VelocityProfile = new VelocityProfile2D(params.species_param[speciesNumber].mvel_z_profile);
                    break;
                default:
                    break;
            }
        }
        else {
            ERROR( "Unsupported geometry : " << params.geometry);
        }
        return VelocityProfile;        
    }
};

#endif

