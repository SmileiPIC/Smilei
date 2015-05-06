#include "ExtFieldProfile2D.h"

using namespace std;

ExtFieldProfile2D::ExtFieldProfile2D(ExtFieldStructure &extfield_struct) : ExtFieldProfile(extfield_struct) {
    
    // ------------------------------------------------
    // Checking for errors & attributing default values
    // ------------------------------------------------
    
/* 354   // if no double_params defined: put to zero and send a warning
    if (my_struct.double_params.size()<1) {
        my_struct.double_params.resize(1);
        my_struct.double_params[0] = 0.0;
        WARNING("double_params not defined for external field: automatically put to 0");
    } 354 */
    
    if (my_struct.profile == "constant") {
        // Constant (uniform) field distribution through all the box
        // nothing to be done
        
    } else if (my_struct.profile == "cosine") {
        // Cosine field distribution
        if (my_struct.int_params.size()<1)
            ERROR("one int_params must be defined for Cosine profile" );
        if (my_struct.double_params.size()<2)
            ERROR("two double_params must be defined for Cosine profile" );
        if ( (my_struct.length_params_x.size()<1) && (my_struct.length_params_y.size()<1)
            && (my_struct.length_params_z.size()<1) )
            ERROR("one length has to be defined for Cosine profile" );
    
    } else if ( (my_struct.profile=="harrisbx") || (my_struct.profile=="harrisby") ) {
        // Harris field distribution (initialization for reconnection)
        if (my_struct.double_params.size()<2)
            ERROR("two double_params must be defined for Harris profile" );
        // requires at least to values for the length_params_x
        if (my_struct.length_params_y.size()<3)
            ERROR("three length_params_y must be defined for Harris profile" );
        // requires at least to values for the length_params_y
        if (my_struct.length_params_y.size()<3)
            ERROR("three length_params_y must be defined for Harris profile" );
        
    }     
    else if (my_struct.profile=="python") {
    }
    else {
        ERROR("unknown or empty profile: " << my_struct.profile );
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
    
    else if (my_struct.profile == "cosine") {
        // -----------------------------------------------------------------
        // Cosine field profile
        // int_params[0]    = direction in which the perturbation is imposed
        // double_params[0] = amplitude of the constant field background
        // double_params[1] = amplitude of the field cosine perturbation
        // length_params_x/y/z[0] = wavelength of the field perturbation
        // -----------------------------------------------------------------
        unsigned int idim = my_struct.int_params[0];
        double B0 = my_struct.double_params[0];
        double dB = my_struct.double_params[1];
        double x = x_cell[idim];
        double k=0;
        if (idim==0) {
            k = 2.0*M_PI / my_struct.length_params_x[0];
        } else if (idim==1) {
            k = 2.0*M_PI / my_struct.length_params_y[0];
        } else if (idim==2) {
            k = 2.0*M_PI / my_struct.length_params_z[0];
        }//endif

        return B0 + dB*cos(k*x);
    }
    
    else if (my_struct.profile == "harrisbx") {
        // ---------------------------------------------------------------
        // Harris initialization for reconnection (field correspond to Bx)
        // double_params[0]   = field amplitude
        // length_params_y[0] = characteristic width of the field
        // length_params_y[1] = position in y of the 1st field node
        // length_params_y[2] = position in y of the 2nd field node
        // ---------------------------------------------------------------
        double B0  = my_struct.double_params[0];
        double dB  = my_struct.double_params[1];
        double L   = my_struct.length_params_y[0];
        double y0  = my_struct.length_params_y[1];
        double y1  = my_struct.length_params_y[2];
        double sgl = my_struct.length_params_x[0];
        double x0  = my_struct.length_params_x[1];
        double x1  = my_struct.length_params_x[2];
        
        double Dx0 = (x_cell[0]-x0)/sgl;
        double Dx1 = (x_cell[0]-x1)/sgl;
        double Dy0 = (x_cell[1]-y0)/sgl;
        double Dy1 = (x_cell[1]-y1)/sgl;
        
        return B0 * ( tanh((x_cell[1]-y0)/L) - tanh((x_cell[1]-y1)/L) - 1.0 )
        -      2.0*dB * Dy0 * exp(-Dx0*Dx0-Dy0*Dy0)
        +      2.0*dB * Dy1 * exp(-Dx1*Dx1-Dy1*Dy1) ;
    }
    
    else if (my_struct.profile == "harrisby") {
        // ---------------------------------------------------------------
        // Harris initialization for reconnection (field correspond to Bx)
        // double_params[0]   = field amplitude
        // length_params_y[0] = characteristic width of the field
        // length_params_y[1] = position in y of the 1st field node
        // length_params_y[2] = position in y of the 2nd field node
        // ---------------------------------------------------------------
        double B0  = my_struct.double_params[0];
        double dB  = my_struct.double_params[1];
        double L   = my_struct.length_params_y[0];
        double y0  = my_struct.length_params_y[1];
        double y1  = my_struct.length_params_y[2];
        double sgl = my_struct.length_params_x[0];
        double x0  = my_struct.length_params_x[1];
        double x1  = my_struct.length_params_x[2];
        
        double Dx0 = (x_cell[0]-x0)/sgl;
        double Dx1 = (x_cell[0]-x1)/sgl;
        double Dy0 = (x_cell[1]-y0)/sgl;
        double Dy1 = (x_cell[1]-y1)/sgl;
        
        return 2.0*dB * Dx0 * exp(-Dx0*Dx0-Dy0*Dy0)
        -      2.0*dB * Dx1 * exp(-Dx1*Dx1-Dy1*Dy1) ;
        
    }    
    else if (my_struct.profile=="python") {
        PyObject *pyresult = PyObject_CallFunction(my_struct.py_profile, const_cast<char *>("dd"), x_cell[0], x_cell[1]);
        if (pyresult == NULL) {
            ERROR("can't evaluate python function");
        }
        double cppresult = PyFloat_AsDouble(pyresult);
        Py_XDECREF(pyresult);
        return cppresult;
    } else {
        return 0;
    }
    
}

