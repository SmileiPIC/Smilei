#include "ExtFieldProfile2D.h"

using namespace std;

ExtFieldProfile2D::ExtFieldProfile2D(ExtFieldStructure &extfield_struct) : ExtFieldProfile(extfield_struct) {
    
    // ------------------------------------------------
    // Checking for errors & attributing default values
    // ------------------------------------------------
    
    // if no double_params defined: put to zero and send a warning
    if (my_struct.double_params.size()<1) {
        my_struct.double_params.resize(1);
        my_struct.double_params[0] = 0.0;
        WARNING("double_params not defined for external field: automatically put to 0");
    }
    
    if (my_struct.profile == "constant") {
        // Constant (uniform) field distribution through all the box
        // nothing to be done
        
    } else if (my_struct.profile == "harris") {
        // Harris field distribution (initialization for reconnection)
        if (my_struct.double_params.size()<1)
            ERROR("one double_params must be defined for Harris profile" );
        // requires at least to values for the length_params_y
        if (my_struct.length_params_y.size()<3)
            ERROR("three length_params_y must be defined for Harris profile" );
        
    } else {
        ERROR("unknown or empty profile :" << my_struct.profile );
    }
    
}


double ExtFieldProfile2D::operator() (vector<double> x_cell) {
    
    if (my_struct.profile == "constant") {
        // ------------------------------------------
        // Constant field profile through all the box
        // double_params[0] = field amplitude
        // ------------------------------------------
        return my_struct.double_params[0];
    }
    
    else if (my_struct.profile == "harris") {
        // ---------------------------------------------------------------
        // Harris initialization for reconnection (field correspond to Bx)
        // double_params[0]   = field amplitude
        // length_params_y[0] = characteristic width of the field
        // length_params_y[1] = position in y of the 1st field node
        // length_params_y[2] = position in y of the 2nd field node
        // ---------------------------------------------------------------
        double A0 = my_struct.double_params[0];
        double L  = my_struct.length_params_y[0];
        double y0 = my_struct.length_params_y[1];
        double y1 = my_struct.length_params_y[2];
        
        return A0 * ( tanh((x_cell[1]-y0)/L) - tanh((x_cell[1]-y1)/L) - 1.0 );
    }
    
    return 0;
}

